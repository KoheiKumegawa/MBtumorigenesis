#----------------------------------------------------------------------------
# 02_scRNAanalysis.R
#----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(scales)
'%ni%' = Negate('%in%')
seu <- readRDS("rds/seu.rds")

#----- clustering -----#
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu, ndims = 50)

seu <- FindNeighbors(seu, dims = 1:30)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.6)

#rename cluster id
seu$Clusters <- factor(paste0("GEx_C", as.numeric(seu$seurat_clusters)) , levels = paste0("GEx_C", c(1:max(as.numeric(seu$seurat_clusters)))))
#color settings from ArchR palettes
cluster_colors <-  ArchR::ArchRPalettes$bear[seq_along(levels(seu$Clusters))] %>% `names<-`(., levels(seu$Clusters))
sample_colors <- ArchR::ArchRPalettes$kelly[c(1:length(table(seu$orig.ident)))] %>% `names<-`(., names(table(seu$orig.ident)))
#UMAP plot
p1 <- DimPlot(seu, reduction = "umap", label = T, group.by = "Clusters", cols = cluster_colors, raster = T) + ggtitle("clusters") + ArchR::theme_ArchR(legendPosition = "right")
p2 <- DimPlot(seu, reduction = "umap", label = F, group.by = "orig.ident", cols = sample_colors, raster = T) + ggtitle("samples") + ArchR::theme_ArchR(legendPosition = "right")

pdf("output/Plots/02_GEx_UMAP.pdf", width = 5.5, height = 5)
p1
p2
dev.off()

#confusion matrix
cM <- table(seu$Clusters, seu$orig.ident)
write.csv(cM, "output/Tables/02_GEx_SampleClusterMatrix.csv", quote = F)

#--------- Marker genes ---------#
#marker genes
markerGenes <- c("Col3a1", #Vascular fibroblast
                 "C1qb",#Microglia
                 "Sox10", "Olig2", #Oligodendrocytes
                 "Aqp4", #Astrocytes
                 "Pax2", "Ascl1", #interneuron progenitors/stem
                 "Gli1", "Gli2", "Ccnd2", "Barhl1", "Atoh1", #GNP or GNP-like tumor
                 "Cntn2", "Grin2b") #differentiated granule cells/tumor cells

#"Nes", #CSC
#"Igf1", "Igf2" #GOI

p3 <- lapply(markerGenes, function(x) FeaturePlot(seu, features = x, order = T, slot = "data") + scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR())
p4 <- DotPlot(seu, features = markerGenes, dot.scale = 7, dot.min = 0.05, group.by = "Clusters") + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke= 1) + ArchR::theme_ArchR() +
  scale_colour_gradientn(colours = viridis(256, option = "A"))

pdf("output/Plots/02_GEx_MarkerExp_UMAPoverlay.pdf", width = 4.5, height = 5)
p3
dev.off()
pdf("output/Plots/02_GEx_MarkerExp_DotPlot.pdf", width = 8, height = 5)
p4
dev.off()

#--------- singleR ---------#
library(SingleR)
library(celldex)
mrd.se <- MouseRNAseqData()
sce <- as.SingleCellExperiment(DietSeurat(seu))
pred <- SingleR(test = sce, ref = mrd.se, assay.type.test=1, labels = mrd.se$label.main)

seu$singleRanno <- pred$pruned.labels
excl.anno <- which(table(seu$singleRanno) < 10) %>% names
seu$singleRanno[seu$singleRanno %in% excl.anno] <- NA

p5 <- DimPlot(seu, reduction = "umap", group.by = "singleRanno", label = T, cols = as.character(ArchR::ArchRPalettes$kelly)) + 
      ggtitle("SingleR annotation") + ArchR::theme_ArchR(legendPosition = "right")

cM2 <- table(seu$Clusters, seu$singleRanno)
df <- reshape2::melt(cM2)
p6 <- ggplot(df, aes(x = factor(Var1, levels = paste0("GEx_C", c(1:15))), y = value, fill = Var2)) + 
  geom_bar(stat = "identity", position = "fill") + ArchR::theme_ArchR() + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = as.character(ArchR::ArchRPalettes$kelly)) + labs(x = "Cluster", y = "%Sample")

pdf("output/Plots/02_GEx_UMAP_SingleRAnnotation.pdf", height = 5.5, width = 6)
p5
dev.off()
pdf("output/Plots/02_GEx_BarPlotSingleRAnnotation.pdf", height = 4, width = 8)
p6
dev.off()
write.csv(cM2, "output/Tables/02_GEx_singleRAnno_ClusterMatrix.csv", quote = F)

#--------- Cluster Annotation ---------#
clusterAnno <- c("GEx_C15" = "VascularFibroblasts",
                 "GEx_C14" = "Oligodendrocytes",
                 "GEx_C13" = "Microglia",
                 "GEx_C12" = "Endothelial",
                 "GEx_C11" = "Astrocytes",
                 "GEx_C10" = "InterneuronProgenitors",
                 "GEx_C9" = "DifferentiatedCells",
                 "GEx_C8" = "DifferentiatedCells",
                 "GEx_C7" = "ProliferativeCells",
                 "GEx_C6" = "ProliferativeCells",
                 "GEx_C5" = "ProliferativeCells",
                 "GEx_C4" = "ProliferativeCells",
                 "GEx_C3" = "ProliferativeCells",
                 "GEx_C2" = "DifferentiatedCells",
                 "GEx_C1" = "ProliferativeCells")
clusterAnno2 <- as.character(seu$Clusters)
for(i in names(clusterAnno)){
  clusterAnno2[which(clusterAnno2 == i)] <- clusterAnno[i]
}
seu$cellType <- clusterAnno2
seu$cellType <- factor(seu$cellType, levels = c("ProliferativeCells", "DifferentiatedCells", "InterneuronProgenitors", "Astrocytes","Microglia","Oligodendrocytes","VascularFibroblasts","Endothelial"))

cellType_colors <- c(ArchR::ArchRPalettes$stallion[c(1:8)]) %>% `names<-`(., levels(seu$cellType))
p7 <- DimPlot(seu, reduction = "umap", label = F, group.by = "cellType", cols = cellType_colors, raster = T) + ggtitle("celltypes") + ArchR::theme_ArchR(legendPosition = "right")

pdf("output/Plots/02_GEx_UMAP_celltype.pdf", width = 6, height = 5)
p7
dev.off()

saveRDS(seu, "rds/seu.rds")

#--------- DEGs ---------#
DEG <- mclapply(paste0("GEx_C", c(1:15)), function(i){
  out <- FindMarkers(seu, group.by = "Clusters", ident.1 = i, only.pos = T)
  out$symbol <- rownames(out)
  out$cluster <- i
  return(out)
}, mc.cores = 8) 
DEG <- do.call(rbind, DEG)
write.csv(DEG, "output/Tables/02_DEG_clusters.csv")

#DEG heatmap
DEG_top50 <- lapply(paste0("GEx_C", c(1:15)), function(i) DEG$symbol[which(DEG$cluster == i)][c(1:50)])
names(DEG_top50) <- paste0("GEx_C", c(1:15))
avg2 <- AverageExpression(seu, return.seurat = T, group.by = "Clusters")
mtx <- avg2@assays$RNA@scale.data[unlist(DEG_top50),]
                                         
lapply(names(DEG_top50), function(i){
  write.table(data.frame(DEG_top50[[i]]), paste0("output/Tables/02_Top50DEG/DEG_top50_", i, ".txt"), col.names = F, row.names = F, quote = F)
  return(NULL)
})

col_fun1 <- colorRamp2(c(-2,0,2), c("blue", "white", "red"))
fh = function(x) hclust(dist(x), method="ward.D2")
rowsplit1 <- factor(lapply(paste0("C", c(1:15)), function(i) rep(i,50)) %>% unlist, levels = paste0("C", c(1:15)))
ht1 <- Heatmap(mtx, name = "Relative expression", cluster_columns = F, cluster_rows = F, 
               show_row_dend = T, show_row_names = T, row_names_gp = gpar(fontsize = 3),
               column_names_gp = gpar(fontsize = 6), col = col_fun1, 
               row_split = rowsplit1, use_raster = T,
               row_title_gp = gpar(fontsize = 6))
p8 <- draw(ht1)
pdf("output/Plots/02_Heatmap_top50DEGs.pdf", height = 10, width = 4)
p8
dev.off()

#----- Prmt4 expression -----#
p9 <- FeaturePlot(seu, features = "Carm1", order = T, slot = "data", raster = T) + scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR()
p10 <- VlnPlot(seu, features = "Carm1", group.by = "orig.ident", cols = sample_colors) + ArchR::theme_ArchR()

pdf("output/Plots/02_UMAPoverlay_Carm1.pdf", height = 5, width = 4)
p9
dev.off()
pdf("output/Plots/02_VlnPlot_Carm1.pdf", height = 4, width = 3)
p10
dev.off()


df <- data.frame(SampleType = seu$orig.ident, Carm1 = seu@assays$RNA@data["Carm1",])
p11 <- ggplot(df, aes(x = SampleType, y = Carm1, fill = SampleType)) + 
       #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.01) + 
       geom_violin() +
       scale_fill_manual(values = c("032103_GNP" = "blue", "RS03056_H" = "orange", "RS03060_T" = "red")) + ArchR::theme_ArchR()
pdf("output/Plots/02_VlnPlot_Carm1_v2.pdf", height = 3, width = 5)
p11
dev.off()

df <- data.frame(SampleType = seu$orig.ident, Carm1 = seu@assays$RNA@data["Carm1",])
df2 <- df[which(df$Carm1 ==0),]
df1 <- df[which(df$Carm1 !=0),]

p12 <- ggplot(df, aes(x = SampleType, y = Carm1, color = SampleType)) +
  geom_dotplot(data = df1, mapping = aes(x = SampleType, y = Carm1, color = SampleType), binaxis = "y", stackdir = "center", binwidth = 0.01) + 
  geom_dotplot(data = df2, mapping = aes(x = SampleType, y = Carm1, color = SampleType), binaxis = "y", stackdir = "center", binwidth = 0.00005) + 
  scale_color_manual(values = c("032103_GNP" = "blue", "RS03056_H" = "orange", "RS03060_T" = "red")) + ArchR::theme_ArchR()

pdf("output/Plots/02_VlnPlot_Carm1_v3.pdf", height = 3, width = 3)
p12
dev.off()
