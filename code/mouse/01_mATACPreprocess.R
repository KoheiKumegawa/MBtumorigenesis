#------------------------------------------------------------------------------
# 01_mATACPreprocess.R
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

sourceCpp(code='
          #include <Rcpp.h>
          using namespace Rcpp;
          using namespace std;
          // [[Rcpp::export]]
          IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
          if(x1.size() != y1.size()){
          stop("width must equal size!");
          }
          IntegerVector x = clone(x1);
          IntegerVector y = clone(y1);
          int n = x.size();
          IntegerVector rx = seq(xmin,xmax);
          IntegerVector ry = seq(ymin,ymax);
          IntegerMatrix mat( ry.size() , rx.size() );
          int xi,yi;
          for(int i = 0; i < n; i++){
          xi = (x[i] - xmin);
          yi = (y[i] - ymin);
          if(yi >= 0 && yi < ry.size()){
          if(xi >= 0 && xi < rx.size()){
          mat( yi , xi ) = mat( yi , xi ) + 1; 
          }
          }
          }
          return mat;
          }')

TSSenrich <- function(
  fragments = NULL,
  TSSgranges = NULL,
  flank = 2000,
  norm = 100,
  smooth = 51,
  range = 50,
  by = NULL
){
  message(sprintf("Convert %s fragments to inserts", length(fragments)))
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), sample = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), sample = mcols(fragments)[,by]))
  
  message(sprintf("center the features"))
  center <- unique(resize(TSSgranges, width = 1, fix = "center", ignore.strand = FALSE))
  
  message(sprintf("get overlaps between the feature and insertions only up to %s bp", flank))
  overlap <- DataFrame(findOverlaps(query = center, subject = inserts, maxgap = flank, ignore.strand = TRUE))
  overlap$strand <- strand(center)[overlap[,1]]
  overlap$name <- mcols(inserts)[overlap[,2],by]
  overlap <- transform(overlap, id=match(name, unique(name)))
  ids <- length(unique(overlap$name))
  
  message(sprintf("calculate distance from TSS"))
  overlap$dist <- NA
  minus <- which(overlap$strand == "-")
  other <- which(overlap$strand != "-")
  overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(inserts[overlap[minus,2]])
  overlap$dist[other] <- start(inserts[overlap[other,2]]) - start(center[overlap[other,1]])
  
  message(sprintf("make insertion matrix"))
  profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
  colnames(profile_mat) <- unique(overlap$name)
  profile <- rowSums(profile_mat)
  
  message(sprintf("calculate TSS enrichments"))
  profile_mat_norm <- apply(profile_mat, 2, function(x) x/mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]))
  profile_norm <- profile/mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])
  profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
  profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)
  max_finite <- function(x){
    suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
  }
  
  #TSS enrichment score                                 
  e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])
  
  #Insersion profiles                                  
  rownames(profile_mat_norm_smooth) <- c(-2000:2000)
  plotDF <- reshape2::melt(profile_mat_norm_smooth)
  colnames(plotDF) <- c("distance", "sample", "normInsert")
  
  #fragment width                                   
  w <- as.numeric(width(fragments))
  
  out <- list(enrichScore = e, plotDF = plotDF, fragmentWidth = w)
  message(sprintf("ATAC QC metrics successfully calculated!"))
  return(out)
}

FragmentGRToInsert <- function(
  gr = NULL,
  name = NULL,
  outputDirBED = NULL
){
  inserts <- c(
    GRanges(seqnames = seqnames(gr), ranges = IRanges(start(gr), start(gr))),
    GRanges(seqnames = seqnames(gr), ranges = IRanges(end(gr), end(gr)))
  )
  #save bed file (2bp)
  d <- data.frame(seqnames = seqnames(inserts), start = start(inserts)-1, end = end(inserts))
  write.table(d, paste0(outputDirBED, "/", name, "_inserts.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}

RunMacs2 <- function(
  inputBedPATH = NULL,
  inputBedName = NULL,
  genome = NULL,
  outputDir = NULL,
  shift = NULL,
  extsize = NULL,
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
  
  commandPeaks <- sprintf(
    "macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all",
    gen, inputBedName, inputBedPATH, outputDir)
  
  if (!is.null(shift) & !is.null(extsize)) {
    commandPeaks <- sprintf("%s --shift %s --extsize %s", commandPeaks, shift, extsize)
  }
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

MakeATACSummarizedExperiment <- function(
  fragmentGRangesList = NULL,
  unionPeaks = NULL,
  blacklist = NULL,
  sampleName = NULL,
  prior.count = 5,
  by = "sample"
){
  fragments <- unlist(as(fragmentGRangesList, "GRangesList"))
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), sample = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), sample = mcols(fragments)[,by])
  )
  overlapDF <- DataFrame(findOverlaps(unionPeaks, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
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
  sparseM.norm <- preprocessCore::normalize.quantiles(sparseM.cpm)
  colnames(sparseM.norm) <- colnames(sparseM.cpm)
  rownames(sparseM.norm) <- unionPeaks$name
  #SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = sparseM, normcounts = sparseM.norm),
                                                   rowRanges = unionPeaks, 
                                                   colData = DataFrame(sample = sampleName))
  return(se)
}

#------------------
# analysis
#------------------

#--------- Importing fragments ---------#
bamFiles  <- list.files("data/", pattern = ".rmdup.bam")
names(bamFiles) <- stringr::str_split(bamFiles, pattern = "_ATAC", simplify = T)[,1]

fragments <- parallel::mclapply(names(bamFiles), 
                                function(x){
                                  out <- bamToFragmentGR(bamPATH = paste0("data/", bamFiles[x]), 
                                                         bamNAME = x, 
                                                         bamFlag = scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE))
                                  return(out)
                                }, mc.cores = 9) %>% `names<-`(., names(bamFiles))

#--------- Calculating and Visualizing Quality ---------#
tss <- TxDb.Mmusculus.UCSC.mm10.knownGene %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique
QualityMetricsATAC <- parallel::mclapply(names(fragments), FUN = function(x){
  message(sprintf("Data processing : %s", x))
  a <- fragments[[x]]
  b <- TSSenrich(fragments = a, TSSgranges = tss, by = "sample")
  return(b)
}, mc.cores = 9)
names(QualityMetricsATAC) <- names(fragments)

#tss score
tssscore <- lapply(QualityMetricsATAC, function(i) i$enrichScore) %>% unlist
tssscore <- data.frame(score = tssscore, tumor = factor(names(tssscore), levels = names(tssscore)))
p1 <- ggplot(tssscore, aes(x = tumor, y = score)) + geom_bar(stat = "identity") + ArchR::theme_ArchR() +
      geom_hline(yintercept = 5, color = "red", lty = "dashed") + labs(x = "Tumor", y = "TSS enrichment score") +
      theme(axis.text.x = element_text(angle = 90))
#normalized insertions
df <- lapply(names(QualityMetricsATAC), function(i) QualityMetricsATAC[[i]]$plotDF) %>% do.call(rbind, .)
p2 <- ggplot(df, aes(x = distance, y = normInsert, color = sample)) + geom_line(size = 0.5) +
  labs(x = "Distance From Center (bp)", y = "Normalized Insertion Profile") + ArchR::theme_ArchR() +
  scale_color_manual(values = ArchR::ArchRPalettes$stallion[c(1:9)] %>% `names<-`(., names(bamFiles))) +
  scale_y_continuous(limits = c(0, max(df$normInsert)*1.05), expand = c(0,0)) +
  scale_x_continuous(limits = c(min(df$distance), max(df$distance)), expand = c(0,0))
#fragment width
df <- lapply(names(QualityMetricsATAC), function(i) data.frame(l = QualityMetricsATAC[[i]]$fragmentWidth, sample = i)) %>% do.call(rbind, .)
p3 <- ggplot(df, aes(x = l, color = sample)) + geom_line(stat = "density", size = 0.5) + ArchR::theme_ArchR() + 
  scale_color_manual(values = ArchR::ArchRPalettes$stallion[c(1:9)] %>% `names<-`(., names(bamFiles))) +
  xlim(0, 600) + ylab("Density") + xlab("Size of fragments (bp)") 

pdf("output/Plots/01_TSSenrichment.pdf", width = 4, height = 3)
p1
dev.off()
pdf("output/Plots/01_NormalizedInsertions.pdf", width = 4, height = 4)
p2
dev.off()
pdf("output/Plots/01_FragmentWidth.pdf", width = 4, height = 4)
p3
dev.off()

#output
saveRDS(QualityMetricsATAC, "rds/QualityMetricsATAC.rds")
parallel::mclapply(names(fragments), 
                   function(x) FragmentGRToBED(gr = fragments[[x]], name = x, outputDir = "output/fragments_bed/"),
                   mc.cores = 9)

#--------- Fragments to Inserts ---------#
#read fragment bed as GRanges
# FragFiles   <- list.files("output/fragments_bed/") %>% `names<-`(., gsub("-fragments.bed", "", .))
# fragments <- parallel::mclapply(FragFiles, function(x){
#   fr <- fread(paste0("output/fragments_bed/", x)) %>% data.frame(.)
#   out <- GRanges(seqnames = fr[,1], IRanges(start = fr[,2]+1, end = fr[,3]))
#   return(out)
# }, mc.cores = 9)

#fragments to inserts
parallel::mclapply(names(fragments), function(x){
  FragmentGRToInsert(gr = fragments[[x]], 
                     name = x, 
                     outputDirBED = "output/inserts_bed/")
}, mc.cores = 9)

#--------- Peak call ---------#
insertFiles <- list.files("output/inserts_bed/") %>% `names<-`(., gsub("_inserts.bed", "", .))
parallel::mclapply(names(insertFiles), function(x){
  RunMacs2(inputBedPATH = paste0("output/inserts_bed/", insertFiles[x]),
           inputBedName = x,
           genome = "mm10",
           outputDir = "output/sample_peaks/",
           shift = -75,
           extsize = 150,
           method = "p",
           cutoff = 0.01)
}, mc.cores = 9)

#--------- Making peak set per sample ---------#
BSgenome   <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome))) %>% 
              GenomeInfoDb::keepStandardChromosomes(., pruning.mode = "coarse")
blacklist  <- rtracklayer::import.bed("ref/mm10-blacklist.v2.bed")
gr_ls <- GenomicRanges::GRangesList(lapply(list.files("output/sample_peaks/", pattern = "summits.bed", full.names = T), function(x) import.bed(x)))
names(gr_ls) <- gsub("_summits.bed", "", list.files("output/sample_peaks/", pattern = "summits.bed"))

gr_ls_proc <- parallel::mclapply(gr_ls, function(x){
  gr <- resize(x, width = 501, fix = "center") %>%
    subsetByOverlaps(., chromSizes, type = "within") %>%
    subsetByOverlaps(., blacklist, invert=TRUE) %>%
    MakeSamplePeakSet(., by = "score")
  mcols(gr)$scorePerMillion <- mcols(gr)$score / (sum(mcols(gr)$score) / 1000000)
  return(gr)
}, mc.cores = 9) %>% GenomicRanges::GRangesList(.)

saveRDS(gr_ls_proc, "rds/gr_ls_proc.rds")

#--------- Making Consensus Peak Set ---------#
gr_cumulative <- MakeSamplePeakSet(unlist(gr_ls_proc), by = "scorePerMillion")
mcols(gr_cumulative)$sampleOverlap <- countOverlaps(gr_cumulative, gr_ls_proc)
reproduciblePeaks <- gr_cumulative[which(mcols(gr_cumulative)$scorePerMillion >= 5 &
                                           mcols(gr_cumulative)$sampleOverlap >= 2 &
                                           seqnames(gr_cumulative) %ni% "chrY")]
names(reproduciblePeaks) <- paste0("mATAC_", c(1:length(reproduciblePeaks)))
mcols(reproduciblePeaks)$origName <- mcols(reproduciblePeaks)$name
mcols(reproduciblePeaks)$name <- names(reproduciblePeaks)

#--------- Constructing ATAC SE ---------#
fragments_named <- lapply(names(fragments), function(x){
  fr <- fragments[[x]]
  mcols(fr)$sample <- x
  return(fr)
})

se <- MakeATACSummarizedExperiment(fragmentGRangesList = fragments_named,
                                   unionPeaks = reproduciblePeaks,
                                   blacklist = blacklist,
                                   sampleName = names(fragments),
                                   by = "sample")
se$sampleType <- stringr::str_split(se$sample, "_", simplify = T)[,1]
saveRDS(se, "rds/mATAC_se.rds")
