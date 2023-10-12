#----------------------------------------------------------------------------
# 03_scATACanalysis.R
#----------------------------------------------------------------------------
library(ArchR)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

arc <- readRDS("rds/arc.rds")

#--------- scATAC ---------#
#clustering
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI") %>% 
  addClusters(., reducedDims = "IterativeLSI", resolution = 0.5, force = T) %>%
  addUMAP(., reducedDims = "IterativeLSI", force = T)

#rename cluster id
arc$Clusters <- factor(paste0("CA_", arc$Clusters), levels = paste0("CA_C", c(1:11)))
arc$GEx_Clusters <- seu$Clusters
arc$cellType <- seu$cellType
#color settings from ArchR palettes
sample_colors <- ArchR::ArchRPalettes$kelly[c(1:length(table(seu$orig.ident)))] %>% `names<-`(., names(table(seu$orig.ident)))
cluster_colors <-  ArchR::ArchRPalettes$bear[seq_along(levels(seu$Clusters))] %>% `names<-`(., levels(seu$Clusters))
cellType_colors <- c(ArchR::ArchRPalettes$stallion[c(1:8)]) %>% `names<-`(., levels(seu$cellType))

#UMAP plot
p1 <- plotEmbedding(arc, name = "Sample", plotAs = "points", size = 0.7, pal = sample_colors)
p2 <- plotEmbedding(arc, name = "Clusters", plotAs = "points", size = 0.7)
p3 <- plotEmbedding(arc, name = "GEx_Clusters", plotAs = "points", size = 0.7, pal = cluster_colors)
p4 <- plotEmbedding(arc, name = "cellType", plotAs = "points", size = 0.7, pal = cellType_colors)

pdf("output/Plots/03_SC_CA_UMAP.pdf", width = 5.5, height = 5)
p1
p2
p3
p4
dev.off()

#confusion matrix
cM <- table(arc$Clusters, arc$Sample)[paste0("CA_C", c(1:11)),]
write.csv(cM, "output/Tables/03_CA_SampleClusterMatrix.csv", quote = F)

#peak call and add peak matrix to arrow files
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters") %>%
       addReproduciblePeakSet(., groupBy = "Clusters", pathToMacs2 = pathToMacs2, method = "q", cutOff = 0.05) %>%
       addPeakMatrix(.)

saveRDS(arc, "rds/arc.rds")

#--------- genomeTrack ---------#
markerGenes <- c("Col3a1", #Vascular fibroblast
                 "C1qb",#Microglia
                 "Sox10", "Olig2", #Oligodendrocytes
                 "Aqp4", #Astrocytes
                 "Pax2", "Ascl1", #interneuron progenitors/stem
                 "Gli1", "Gli2", "Ccnd2", "Barhl1", "Atoh1", #GNP or GNP-like tumor
                 "Cntn2", "Grin2b") #differentiated granule cells/tumor cells

p5 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "cellType", 
  useGroups = c("ProliferativeCells", "DifferentiatedCells", "InterneuronProgenitors", "Astrocytes","Microglia","Oligodendrocytes","VascularFibroblasts","Endothelial"),
  geneSymbol = markerGenes, 
  pal = cellType_colors,
  upstream = 5000, downstream = 5000)

plotPDF(p5, name = "03_GenomeTrack_Markers.pdf", width = 5, height = 5, addDOC = F)
pdf("output/Plots/03_GenomeTrack_Markers.pdf", width = 5, height = 5)
for(i in c(1:length(p5))){grid::grid.draw(p5[[i]])}
dev.off()

#--------- cluster heatmap ---------#
cM2 <- confusionMatrix(arc$GEx_Clusters, arc$Clusters)[paste0("GEx_C", c(1:15)), paste0("CA_C", c(1:11))] %>% as.matrix()

mtx <- cM2 / rowSums(cM2) * 100
apply(mtx, 2, function(i) order(i,decreasing = T)[1])


col_fun1 <- colorRamp2(c(0,50,100), c("white", "yellow", "red"))
ht1 <- Heatmap(mtx[c(13,15,14,11,10,2,8,9,12,6,4,7,3,5,1),], 
               cluster_rows = F, cluster_columns = F,
               col = col_fun1)
p6 <- draw(ht1)

pdf("output/Plots/03_Heatmap_ATAC_RNA_Cluster.pdf", width = 4, height = 3)
p6
dev.off()

#--------- UMAP ---------#
arc <- addImputeWeights(arc)
p7 <- plotEmbedding(arc, 
                    colorBy = "GeneScoreMatrix", 
                    name = markerGenes, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 0.5, rastr = T)
pdf("output/Plots/03_UMAPoverlay_markers.pdf", width = 4, height = 4)
p7
dev.off()

### after merge expression data
p8 <- plotEmbedding(arc, 
                    colorBy = "GeneExpressionMatrix", 
                    name = markerGenes, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc),
                    plotAs = "points", size = 0.5, rastr = T)
pdf("output/Plots/03_UMAPoverlay_GEx_markers.pdf", width = 4, height = 4)
p8
dev.off()

