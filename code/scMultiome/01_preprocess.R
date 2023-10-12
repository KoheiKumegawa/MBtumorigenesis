#----------------------------------------------------------------------------
# 01_preprocess
#----------------------------------------------------------------------------
library(ArchR)
library(Seurat)

#--------- ArchR Proprocess ---------#
addArchRThreads(threads = 16) 

#sample assignment
sampleName <- paste0("data/", c("RS03056", "RS03060", "032103"), "_fragments.tsv.gz")
names(sampleName) <- c("RS03056_H", "RS03060_T", "032103_GNP")
outFile <- as.character(sampleName)

#make arrow files
addArchRGenome("mm10")
ArrowFiles = character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = outFile,
                               sampleNames = names(sampleName), 
                               minTSS = 4, minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE)

#doublet detection
doubScores <- addDoubletScores(ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

#make ArchR project and filter doublets!
pre_arc <- ArchRProject(ArrowFiles, outputDirectory = "output", copyArrows = F) %>% filterDoublets(.)

#check TSS score and number of fragments per sample
df <- data.frame(nFrags = log10(pre_arc$nFrags), TSSe = pre_arc$TSSEnrichment, Sample = pre_arc$Sample)
p1 <- lapply(names(sampleName), function(x){
  tmp <- df[which(df$Sample == x),]
  p <- ggplot(tmp, aes(x = nFrags, y = TSSe)) + geom_hex(bins = 100) + theme_classic() + 
    scale_fill_viridis_c(trans = "log") + lims(x = c(0, NA), y = c(0, NA)) + 
    geom_vline(xintercept = c(log10(1000), log10(2500), log10(3000))) +
    labs(x = "Log10 Unique Fragments", y = "TSS Enrichment score") + ggtitle(x)
})
plotPDF(p1, name = "01_QC_preprocess.pdf", ArchRProj = pre_arc, addDOC = FALSE, width = 5, height = 5)

#quality filter
atac_valid <- pre_arc$cellNames[which(pre_arc$TSSEnrichment > 6 & pre_arc$nFrags > 2500)]

#--------- Seurat Proprocess ---------#
#function
makeSeuratObj <- function(
  Dir_10x = NULL,
  SampleName = NULL
){
  data <- Read10X(Dir_10x)
  out <- CreateSeuratObject(counts = data$`Gene Expression`)
  
  #add %mitochondrial RNA 
  out[["percent.mt"]] <- PercentageFeatureSet(out, pattern = "mt-") #mouse "mt-"
  #rename samples
  out$orig.ident <- SampleName
  out <- RenameCells(out, new.names = paste0(SampleName, "#", colnames(out)))
  
  return(out)
}

#sample assignment
sampleName <- paste0("data/", c("RS03056", "RS03060", "032103"), "_filtered_feature_bc_matrix")
names(sampleName) <- c("RS03056_H", "RS03060_T", "032103_GNP")

#make Seurat object and change cell names corresponding to scATAC data
pre_seu_ls <- mclapply(seq_along(sampleName), function(x){
  a <- as.character(sampleName[x])
  b <- names(sampleName)[x]
  out <- makeSeuratObj(Dir_10x = a, SampleName = b)
  return(out)
}, mc.cores = 3)

#merge Seurat object
pre_seu <- merge(merge(pre_seu_ls[[1]], pre_seu_ls[[2]]), pre_seu_ls[[3]])

#quality filter
rna_valid <- colnames(pre_seu)[which(pre_seu$nFeature_RNA > 500 & pre_seu$nFeature_RNA < 9000 & pre_seu$percent.mt < 20)]

#--------- multi-omics quality valid cells ---------#
valid_cells <- intersect(rna_valid, atac_valid) # 15,454 cells

#final filter
arc <- pre_arc[valid_cells]
seu <- pre_seu[,valid_cells]

#check----------------
min(arc$TSSEnrichment)
min(arc$nFrags)
min(seu$nFeature_RNA)
max(seu$nFeature_RNA)
max(seu$percent.mt)
#---------------------

#--------- save objects ---------#
saveRDS(pre_arc, "rds/pre_arc.rds")
saveRDS(pre_seu, "rds/pre_seu.rds")
saveRDS(arc, "rds/arc.rds")
saveRDS(seu, "rds/seu.rds")
