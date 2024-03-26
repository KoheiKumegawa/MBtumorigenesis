#----------------------------------------------------------------------------
# 05_Trajectory.R
#----------------------------------------------------------------------------
library(ArchR)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(viridis)
arc <- readRDS("rds/arc.rds")

addArchRGenome("mm10")
addArchRThreads(threads = 16) 

arc2 <- arc[which(arc$Clusters %in% paste0("CA_C", c(6,7,9,10,11)))]

#--------- clustering ---------#
#LSI-ATAC
arc2 <- addIterativeLSI(
  ArchRProj = arc2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

#LSI-RNA
arc2 <- addIterativeLSI(
  ArchRProj = arc2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

#Combined Dims
arc2 <- addCombinedDims(arc2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

#UMAPs
arc2 <- addUMAP(arc2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
arc2 <- addUMAP(arc2, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
arc2 <- addUMAP(arc2, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
arc2 <- addClusters(arc2, reducedDims = "LSI_RNA", name = "Clusters2", resolution = 0.8, force = TRUE)

#Plot Embedding
sample_colors <- ArchR::ArchRPalettes$kelly[c(1:3)] %>% `names<-`(., c("GNP", "PNC", "Tumor"))
p1 <- plotEmbedding(arc2, name = "Clusters2", embedding = "UMAP_ATAC", size = 1.5)
p2 <- plotEmbedding(arc2, name = "Clusters2", embedding = "UMAP_RNA", size = 1.5)
p3 <- plotEmbedding(arc2, name = "Clusters2", embedding = "UMAP_Combined", size = 1.5)
p4 <- plotEmbedding(arc2, name = "Sample2", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F, pal = sample_colors)
pdf("output/Plots/08_UMAP_ProlifDiff.pdf", width = 5, height = 5)
p1
p2
p3
p4
dev.off()

#motif score
arc2 <- addImputeWeights(arc2)
p5 <- plotEmbedding(arc2, 
                    colorBy = "GeneExpressionMatrix", 
                    name = c("Nfia", "Nfib"), 
                    embedding = "UMAP_RNA", 
                    imputeWeights = getImputeWeights(arc2), 
                    pal = ArchRPalettes$sambaNight,
                    log2Norm = T, colorLimit =c(0,7),
                    plotAs = "points", size = 0.5, rastr = T)
p6 <- plotEmbedding(arc2, 
                    colorBy = "GeneScoreMatrix", 
                    name = c("Nfia", "Nfib"), 
                    embedding = "UMAP_RNA", 
                    imputeWeights = getImputeWeights(arc2), 
                    plotAs = "points", size = 0.5, rastr = T)
p7 <- plotEmbedding(arc2, 
                    colorBy = "homerMatrix", 
                    name = c("z:NF1.halfsite.CTF_175","z:NF1.CTF_176"), 
                    embedding = "UMAP_RNA", 
                    imputeWeights = getImputeWeights(arc2), 
                    plotAs = "points", size = 0.5, rastr = T)
pdf("output/Plots/08_UMAP_NFI.pdf", width = 5, height = 5)
p5
p6
p7
dev.off()

p8 <- plotEmbedding(arc2, 
                    colorBy = "GeneExpressionMatrix", 
                    name = "Smo", 
                    embedding = "UMAP_RNA", log2Norm = T, 
                    imputeWeights = getImputeWeights(arc2), 
                    pal = ArchRPalettes$sambaNight, colorLimit =c(0,1.5),
                    plotAs = "points", size = 0.5, rastr = T)
p9 <- plotEmbedding(arc2, 
                    colorBy = "GeneExpressionMatrix", 
                    name = "Neurod1", 
                    embedding = "UMAP_RNA", log2Norm = T, 
                    imputeWeights = getImputeWeights(arc2), 
                    pal = ArchRPalettes$sambaNight, colorLimit =c(0,3),
                    plotAs = "points", size = 0.5, rastr = T)

pdf("output/Plots/08_UMAP_Smo_Neurod1_Exp.pdf", width = 5, height = 5)
p8
p9
dev.off()

markerGenes <- c("Gli1", "Gli2", "Ccnd2", "Barhl1", "Atoh1", #GNP or GNP-like tumor
                 "Cntn2", "Grin2b") #differentiated granule cells/tumor cells
p10 <- plotEmbedding(arc2, 
                    colorBy = "GeneExpressionMatrix", 
                    name = markerGenes, 
                    embedding = "UMAP_RNA", log2Norm = T, 
                    imputeWeights = NULL, 
                    pal = ArchRPalettes$sambaNight, 
                    plotAs = "points", size = 0.5, rastr = T)
pdf("output/Plots/08_UMAP_Markers_Exp.pdf", width = 5, height = 5)
p10
dev.off()

#--------- trajectry ---------#
trajectory <- paste0("C", c(4,2,7,9,10,13,12,8))
trajectory
arc2 <- addTrajectory(
  ArchRProj = arc2, 
  name = "Trajectory", 
  groupBy = "Clusters2",
  trajectory = trajectory, 
  embedding = "UMAP_RNA", 
  force = TRUE
)
p11 <- plotTrajectory(arc2, embedding = "UMAP_RNA", trajectory = "Trajectory", colorBy = "cellColData", name = "Trajectory", continuousSet = "horizonExtra")

# trajectory expression/genescore/motifs
p12 <- plotTrajectory(arc2, embedding = "UMAP_RNA", trajectory = "Trajectory", colorBy = "homerMatrix", name = "z:NF1.CTF_176", continuousSet = "horizonExtra")
p13 <- plotTrajectory(arc2, embedding = "UMAP_RNA", trajectory = "Trajectory", colorBy = "GeneExpressionMatrix", name = "Nfia", 
                      log2Norm = T, continuousSet = "horizonExtra", quantCut = c(0,1), imputeWeights = NULL)
p14 <- plotTrajectory(arc2, embedding = "UMAP_RNA", trajectory = "Trajectory", colorBy = "GeneExpressionMatrix", name = "Nfib", 
                      log2Norm = T, continuousSet = "horizonExtra", quantCut = c(0,1), imputeWeights = NULL)
p15 <- plotTrajectory(arc2, embedding = "UMAP_RNA", trajectory = "Trajectory", colorBy = "GeneExpressionMatrix", name = "Smo", 
                      log2Norm = T, continuousSet = "horizonExtra", quantCut = c(0,1), imputeWeights = NULL)
p16 <- plotTrajectory(arc2, embedding = "UMAP_RNA", trajectory = "Trajectory", colorBy = "GeneExpressionMatrix", name = "Neurod1", 
                      log2Norm = T, continuousSet = "horizonExtra", quantCut = c(0,1), imputeWeights = NULL)

pdf("output/Plots/08_TrajectryPlots.pdf", width = 5, height = 5)
p11
p12
p13
p14
p15
p16
dev.off()

saveRDS(arc2, "rds/08_arc2.rds")

#--------- NFI bound site dynamics ---------#
arc2 <- readRDS("rds/08_arc2.rds")
peakSet <- arc2@peakSet
names(peakSet) <- paste0("peak_", c(1:length(peakSet)))

PeakClusters2 <- getMarkerFeatures(
  ArchRProj = arc2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
## Caution! Line up chromosomal order!! (chr1 -> chr10 -> chr11... to chr1 -> chr2 -> chr3)
rowranges1 <- GRanges(seqnames = rowData(PeakClusters2)$seqnames, 
                      IRanges(start = rowData(PeakClusters2)$start, end = rowData(PeakClusters2)$end))
idx <- findOverlaps(peakSet, rowranges1) %>% subjectHits()
PeakClusters2_mod <- PeakClusters2[idx, ]
rownames(PeakClusters2_mod) <- names(peakSet)

#NFI binding sites from CUT&Tag data
mNFI_se <- readRDS("../mouse/rds/mNFI_se.rds")
idx2 <- findOverlaps(peakSet, rowRanges(mNFI_se)) %>% queryHits() %>% unique() #scATAC:125746, NFI_C&T:99387, overlaps:66556

PeakClusters2_mod <- PeakClusters2_mod[names(peakSet)[idx2], paste0("C", c(4,2,7,9,10,13,12,8))]
mtx_peaks <- assays(PeakClusters2_mod)$Mean
mtx_peaks <- mtx_peaks[which(rowSums(mtx_peaks) != 0),]
mtx_peaks_z <- t(scale(t(mtx_peaks)))

#DBRs
DBRs <- lapply(c("common", "DBR_GNP", "DBR_PNC", "DBR_PT", "DBR_Tumor"), function(i){
                 gr <- import.bed(paste0("../mouse/output/output_bed/DBR/",i,".bed"))
                 out <- subsetByOverlaps(peakSet, gr) %>% names()
                 return(out)}
               ) 
names(DBRs) <- c("common", "GNP", "PNC", "PT", "Tumor")

DBR_anno <- rep("NA", nrow(mtx_peaks_z))
names(DBR_anno) <- rownames(mtx_peaks_z)
for(i in names(DBRs)){DBR_anno[DBRs[[i]]] <- i}
DBR_anno <- DBR_anno[rownames(mtx_peaks_z)]

km <- kmeans(mtx_peaks_z, centers = 5)
col_fun1 <- colorRamp2(c(-2,-1,0,1,2), viridis(5, option = "D"))
fh = function(x) hclust(dist(x), method="ward.D2")
ra1 = rowAnnotation(DBR = DBR_anno, 
                    col = list(DBR = c("common" = "darkgray", "GNP" = "blue", "PNC" = "orange", "PT" = "darkgreen", "Tumor" = "red", "NA" = "white")))
ht1 <- Heatmap(mtx_peaks_z, name = "Accessibility z-score", row_split = km$cluster,
               cluster_columns = F, cluster_rows = F, show_row_dend = F, show_row_names = F, 
               col = col_fun1, left_annotation = ra1, use_raster = T)
p20 <- draw(ht1)
pdf("output/Plots/08_Heatmap_NFIBoundSite_access_z_kmeans.pdf", height = 6, width = 4)
p20
dev.off()

df <- table(data.frame(DBR_anno, km$cluster)) %>% as.data.frame.matrix()
write.csv(df, "output/Tables/08_kmpeaks_DBRs_confusionMat.csv")

df <- df[-3,]
df2 <- t(df) / as.numeric(table(km$cluster))
df2 <- as.matrix(df2*100)
df2 <- df2[c(5,2,4,3,1),]

#overlap heatmap
col_fun2 <- colorRamp2(c(0,2,4,6,10,12,14), c("white", rev(viridis::viridis(6, option = "G"))))
ht2 <- Heatmap(df2,
               name = "Percent of overlaps", col = col_fun2, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = T, show_row_dend = F, 
               row_title_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 10),
               column_title = "NFI DBR", column_title_side = "top", row_title = "Kmeans", use_raster = F)
p21 <- draw(ht2)
pdf("output/Plots/08_HeatmapPercent_kmeans_DBR.pdf", height = 3, width = 3.5)
p21
dev.off()

#output as bed
lapply(c(1:5), function(i){
  g <- peakSet[names(which(km$cluster == i))]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/08_NFIBoundSiteAR/BS_km", i, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
})

#overlaps with DORC
DORC_p2g <- readRDS("rds/04_DORC_p2g.rds")
count_overlaps <- lapply(c(1:5), function(i){
  gr <- import.bed(paste0("output/output_bed/08_NFIBoundSiteAR/BS_km",i,".bed"))
  out <- countOverlaps(DORC_p2g, gr)
  return(out)
})
count_overlaps <- do.call(cbind, count_overlaps)
colnames(count_overlaps) <- paste0("km_", c(1:5))
count_overlaps <- count_overlaps[,paste0("km_", c(5,2,4,3,1))]
count_overlaps <- count_overlaps[order(rownames(count_overlaps)),]

write.csv(count_overlaps, "output/Tables/08_NFIBoundSite_DORC_overlaps.csv")
# count_overlaps <- read.csv("output/Tables/08_NFIBoundSite_DORC_overlaps.csv", row.names = 1)

#No.overlaps with key gene DORC
df <- count_overlaps[c("Smo","Ptch1", "Ccnd2", "Mcm2", "Cntn2", "Neurod1"), ]
p22 <- lapply(c("Smo","Ptch1", "Ccnd2", "Mcm2", "Cntn2", "Neurod1"), function(i){
  df2 <- df[i,] %>% reshape2::melt(.)
  df2$km <- factor(rownames(df2), levels = paste0("km_", c(5,2,4,3,1)))
  out <- ggplot(df2, aes(x = km, y = value)) + geom_bar(stat = "identity") + theme_ArchR() + 
    labs(x = "k-means clusters", y = paste0("No.overlaps with ", i, " DORC"))
  return(out)
})
pdf("output/Plots/08_Barplot_kmeans_overlap_DORC.pdf", height = 3, width = 5)
p22
dev.off()

#Example visualization, Smo
smo_overlaps <- lapply(c(1:5), function(i){
  gr <- import.bed(paste0("output/output_bed/08_NFIBoundSiteAR/BS_km",i,".bed"))
  out <- subsetByOverlaps(gr, DORC_p2g$Smo)
  return(out)
}) %>% GRangesList(.)
names(smo_overlaps) <- paste0("km_", c(1:5))
smo_overlaps <- smo_overlaps[paste0("km_", c(5,2,4,3,1))]
lapply(smo_overlaps, length)

p2g_500k_gr <- getPeak2GeneLinks(arc2, corCutOff = 0, FDRCutOff = 0.1, resolution = 500, returnLoops = T)
start_gr <- resize(p2g_500k_gr$Peak2GeneLinks, width = 1, fix = "start")
end_gr <- resize(p2g_500k_gr$Peak2GeneLinks, width = 1, fix = "end")
smo_gr <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == "Smo")]
smo_gr <- resize(smo_gr, width = 500, fix = "start")
n1 <- findOverlaps(start_gr, smo_gr) %>% queryHits()
n2 <- findOverlaps(end_gr, smo_gr) %>% queryHits()
gr <- p2g_500k_gr$Peak2GeneLinks[unique(sort(c(n1,n2)))]
gr$FDR2 <- -log10(gr$FDR)
gr$value <- gr$FDR2

p23 <- lapply(seq_along(smo_overlaps), function(i)
  plotBrowserTrack(
    ArchRProj = arc2, 
    groupBy = "Clusters2", useGroups = paste0("C", c(4,2,7,9,10,13,12,8)),
    geneSymbol = "Smo", 
    region = GRanges(seqnames = "chr6", ranges = IRanges(start = 29230001, end = 30200000)),
    loops = gr, features = smo_overlaps[[i]]))
plotPDF(p23, name = "08_GenomeTrack_Smo_kmeans.pdf", width = 10, height = 8, addDOC = F)

#DBR
DBRs
smo_overlaps_DBR <- lapply(names(DBRs), function(i){
  gr <- peakSet[DBRs[[i]]]
  out <- subsetByOverlaps(gr, unlist(smo_overlaps))
  return(out)
}) %>% GRangesList(.)
names(smo_overlaps_DBR) <- names(DBRs)

p24 <- lapply(seq_along(smo_overlaps_DBR), function(i)
  plotBrowserTrack(
    ArchRProj = arc2, 
    groupBy = "Clusters2", useGroups = paste0("C", c(4,2,7,9,10,13,12,8)),
    geneSymbol = "Smo", 
    region = GRanges(seqnames = "chr6", ranges = IRanges(start = 29230001, end = 30200000)),
    loops = gr, features = smo_overlaps_DBR[[i]]))
plotPDF(p24, name = "08_GenomeTrack_Smo_kmeans_DBR.pdf", width = 10, height = 8, addDOC = F)
