#------------------------------------------------------------------------------
# 04_mNFIPreprocess.R
#------------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(Rcpp)
library(Rsamtools)
library(data.table)
library(ggplot2)
library(rtracklayer)
library(SummarizedExperiment)
'%ni%' <- Negate("%in%")

#------------------
# functions
#------------------
bamToFragmentGR <- function(
  bamPATH = NULL,
  bamNAME = NULL,
  offsetPlus = 4,
  offsetMinus = -5,
  bamFlag = NULL
){
  if(is.null(bamPATH)){
    stop("Please set PATH to bam files")
  }
  if(is.null(bamNAME)){
    stop("No input bamNAME; please recheck your input")
  }
  if(is.null(bamFlag)){
    stop("Please set bamFlag using Rsamtools's scanBamFlag!")
  }
  
  #1. Read In Bam File
  sF <- scanBam(bamPATH, param = ScanBamParam(flag = bamFlag, what = c("rname","pos", "isize")))[[1]]
  
  #2. Make Fragment Table
  dt <- data.table(seqnames = sF$rname, start = sF$pos + offsetPlus, end = sF$pos + abs(sF$isize) - 1 + offsetMinus)
  
  #3. Make fragment Granges and remove unwanted chromosomes
  gr <- GRanges(seqnames = dt$seqnames, IRanges(start = dt$start, end = dt$end))
  idy = which(seqnames(gr) %in% seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  gr <- gr[-idy]
  gr <- dropSeqlevels(gr, seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  mcols(gr) <- DataFrame(sample = bamNAME)
  
  #4. output Granges List
  return(gr)
}

FragmentGRToBED <- function(
  gr = NULL,
  name = NULL,
  outputDir = NULL
){
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0(outputDir, "/", name, "-fragments.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}

RunMacs2 <- function(
  inputBamPATH = NULL,
  inputBamName = NULL,
  genome = NULL,
  outputDir = NULL,
  method = c("p", "q"),
  cutoff = 0.05
){
  if(genome %ni% c("hg19", "hg38", "mm10")){
    stop("Please set genome as hg19, hg38 and mm10!")
  }
  if(genome %in% c("hg19", "hg38")){
    gen <- "hs"
  }
  if(genome == "mm10"){
    gen <- "mm"
  }
  
  commandPeaks <- sprintf("macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BAM --call-summits --keep-dup auto", 
                          gen, inputBamName, inputBamPATH, outputDir)
  
  if (tolower(method) == "p") {
    commandPeaks <- sprintf("%s -p %s", commandPeaks, cutoff)
  } else {
    commandPeaks <- sprintf("%s -q %s", commandPeaks, cutoff)
  }
  
  message("Running Macs2...")
  message(commandPeaks)
  system(commandPeaks, intern = TRUE)
  
  return(NULL)
}

MakeSamplePeakSet <- function(gr, by = "score"){
  #nonOverlappingGRanges
  stopifnot(by %in% colnames(mcols(gr)))
  
  #function for picking up most significant peaks
  clusterGRanges <- function(gr, by = "score"){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = TRUE),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  #iteration of filtering overlapping peaks
  i <-  0
  gr_converge <- gr
  while(length(gr_converge) > 0){
    i <-  i + 1
    gr_selected <- clusterGRanges(gr = gr_converge, by = by)
    gr_converge <- subsetByOverlaps(gr_converge, gr_selected, invert=TRUE) #blacklist selected gr
    if(i == 1){ #if i=1 then set gr_all to clustered
      gr_all <- gr_selected
    }else{
      gr_all <- c(gr_all, gr_selected)
    }
  }
  gr_all <- sort(sortSeqlevels(gr_all))
  return(gr_all)
}

MakeCTSummarizedExperiment <- function(
  fragmentGRangesList = NULL,
  unionPeaks = NULL,
  blacklist = NULL,
  sampleName = NULL,
  prior.count = 1,
  by = "sample"
){
  fragments <- unlist(as(fragmentGRangesList, "GRangesList"))
  overlapDF <- DataFrame(findOverlaps(unionPeaks, fragments, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(fragments)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(unionPeaks), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  sparseM <- sparseM[, sampleName]
  rownames(sparseM) <- unionPeaks$name
  sparseM.cpm <- edgeR::cpm(sparseM, log = TRUE, prior.count = prior.count)
  rownames(sparseM.cpm) <- unionPeaks$name
  #SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = sparseM, log2cpm = sparseM.cpm),
                                                   rowRanges = unionPeaks, 
                                                   colData = DataFrame(sample = sampleName))
  return(se)
}

#------------------
# analysis
#------------------

#--------- Importing fragments ---------#
bamFiles  <- list.files("data/", pattern = ".mapped.bam")
names(bamFiles) <- stringr::str_split(bamFiles, pattern = "_align", simplify = T)[,1]
fragments <- parallel::mclapply(names(bamFiles), 
                                function(x){
                                  out <- bamToFragmentGR(bamPATH = paste0("data/", bamFiles[x]), 
                                                         bamNAME = x,
                                                         bamFlag = scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE))
                                  return(out)
                                }, mc.cores = 9) %>% `names<-`(., names(bamFiles))

#--------- Calculating and Visualizing Quality ---------#
#fragment size
fragmentsSizeCT <- lapply(fragments, function(x) width(x))
df <- lapply(names(fragmentsSizeCT), function(i) data.frame(l = fragmentsSizeCT[[i]], sample = i)) %>% do.call(rbind, .)
p1 <- ggplot(df, aes(x = l, color = sample)) + geom_line(stat = "density", size = 0.5) + ArchR::theme_ArchR() + 
      scale_color_manual(values = ArchR::ArchRPalettes$stallion[c(1:18)] %>% `names<-`(., names(fragmentsSizeCT))) +
      xlim(0, 600) + ylab("Density") + xlab("Size of fragments (bp)") 

pdf("output/Plots/04_FragmentWidth_NFI.pdf", width = 4, height = 4)
p1
dev.off()

#output
parallel::mclapply(names(fragments), 
                   function(x) FragmentGRToBED(gr = fragments[[x]], name = x, outputDir = "output/nfi_fragments_bed/"), 
                   mc.cores = 9)

#--------- Peak call ---------#
parallel::mclapply(names(bamFiles), function(x){
  RunMacs2(inputBamPATH = paste0("data/", bamFiles[x]),
           inputBamName = x,
           genome = "mm10",
           outputDir = "output/nfi_sample_peaks/",
           method = "p",
           cutoff = 0.01)
}, mc.cores = 9)

#--------- Making peak set per sample ---------#
BSgenome   <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome))) %>% 
              GenomeInfoDb::keepStandardChromosomes(., pruning.mode = "coarse")
blacklist  <- rtracklayer::import.bed("ref/mm10-blacklist.v2.bed")
gr_ls <- GenomicRanges::GRangesList(lapply(list.files("output/nfi_sample_peaks/", pattern = "summits.bed", full.names = T), function(x) import.bed(x)))
names(gr_ls) <- gsub("_summits.bed", "", list.files("output/nfi_sample_peaks/", pattern = "summits.bed"))

gr_ls_proc <- parallel::mclapply(gr_ls, function(x){
  gr <- resize(x, width = 501, fix = "center") %>%
    subsetByOverlaps(., chromSizes, type = "within") %>%
    subsetByOverlaps(., blacklist, invert=TRUE) %>%
    MakeSamplePeakSet(., by = "score")
  mcols(gr)$scorePerMillion <- mcols(gr)$score / (sum(mcols(gr)$score) / 1000000)
  return(gr)
}, mc.cores = 9) %>% GenomicRanges::GRangesList(.)

saveRDS(gr_ls_proc, "rds/nfi_gr_ls_proc.rds")

#--------- Making Consensus Peak Set ---------#
gr_cumulative <- MakeSamplePeakSet(unlist(gr_ls_proc), by = "scorePerMillion")
mcols(gr_cumulative)$sampleOverlap <- countOverlaps(gr_cumulative, gr_ls_proc)
reproduciblePeaks <- gr_cumulative[which(mcols(gr_cumulative)$scorePerMillion >= 10 &
                                           mcols(gr_cumulative)$sampleOverlap >= 2 &
                                           seqnames(gr_cumulative) %ni% c("chrY", "chrM"))]
names(reproduciblePeaks) <- paste0("mNFI_", c(1:length(reproduciblePeaks)))
mcols(reproduciblePeaks)$name <- names(reproduciblePeaks)

#--------- Constructing NFI SE ---------#
fragments_named <- lapply(names(fragments), function(x){
  fr <- fragments[[x]]
  mcols(fr)$sample <- x
  return(fr)
})
se <- MakeCTSummarizedExperiment(fragmentGRangesList = fragments_named,
                                 unionPeaks = reproduciblePeaks,
                                 blacklist = blacklist,
                                 sampleName = names(fragments),
                                 by = "sample",
                                 prior.count = 5)
se$sampleType <- stringr::str_split(se$sample, "_", simplify = T)[,1]
se$ab <- stringr::str_split(se$sample, "_", simplify = T)[,2]
se$rep <- stringr::str_split(se$sample, "_", simplify = T)[,3]

saveRDS(se, "rds/mNFI_se.rds")
