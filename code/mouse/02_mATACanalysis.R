#------------------------------------------------------------------------------
# 02_mATACanalysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
source("code/edgeR_PairwiseFunction.R")
se <- readRDS("rds/mATAC_se.rds")

#--------- Annotate Peaks ---------#
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(rowRanges(se), tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Mm.eg.db")
mcols(se) <- mcols(as.GRanges(peakAnno))

pdf("output/Plots/02_peakAnnoPie.pdf", width = 7, height = 4)
plotAnnoPie(peakAnno) 
dev.off()

#--------- PCA ---------#
pca1 <- prcomp(t(assays(se)$normcounts))
summary(pca1) # PC1: 0.4904, PC2: 0.2506
df <- pca1$x[, c(1:2)] %>% data.frame
df$sample <- se$sampleType

p1 <- ggplot(df, aes(x = PC1, y = PC2, color = sample)) + geom_point() + ArchR::theme_ArchR() + 
      scale_color_manual(values = c(GNP = "blue", PNC = "orange", Tumor = "red")) +
      labs(x = "PC1 (49.0% variance)", y = "PC2 (25.1% variance)")
pdf("output/Plots/02_PCA.pdf", width = 4, height = 4.5)
p1
dev.off()

#--------- Correlation analysis ---------#
corr_promoter <- cor(assays(se[mcols(se)$annotation == "Promoter", ])$normcounts)
corr_distal   <- cor(assays(se[mcols(se)$annotation != "Promoter", ])$normcounts)

fh = function(x) hclust(dist(x), method="ward.D2")
col_fun1 <- colorRamp2(c(0.4, 0.7, 1.0), c("blue", "white", "red"))
ht1 <- Heatmap(corr_promoter, name = "Pearson's correlation", col = col_fun1, cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, show_column_dend = F, show_row_names = T, show_row_dend = T, row_names_gp = gpar(fontsize = 10), 
               column_title = paste0("Promoter elements (N = ", length(which(mcols(se)$annotation == "Promoter")),")"),
               column_title_side = "top", row_title = "")
p2 <- draw(ht1)
ht2 <- Heatmap(corr_distal, name = "Pearson's correlation", col = col_fun1, cluster_rows = fh, cluster_columns = fh,
               show_column_names = F, show_column_dend = F, show_row_names = T, show_row_dend = T, row_names_gp = gpar(fontsize = 10), 
               column_title = paste0("Distal elements (N = ", length(which(mcols(se)$annotation != "Promoter")),")"),
               column_title_side = "top", row_title = "")
p3 <- draw(ht2)
pdf("output/Plots/02_corHeatmap.pdf", width = 5, height = 3)
p2
p3
dev.off()

#--------- Identify differential accessible regions ---------#
#commonly accessible
stat1 <- data.frame(row.names = mcols(se)$name, median = rowMedians(assays(se)$normcounts), variance = rowVars(assays(se)$normcounts))
p4 <- ggplot(stat1, aes(x = variance, y = median)) + geom_hex(bins = 100) + ArchR::theme_ArchR() + 
      scale_fill_viridis_c(trans = "log") + labs(x = "Variance", y = "Median accessibility") + xlim(0,4) +
      geom_hline(yintercept = 4.5, color = "red", lty = "dashed") + geom_vline(xintercept = 0.25, color = "red", lty = "dashed") 
DAR_common <- rownames(stat1[which(stat1$variance < 0.25 & stat1$median > 4.5), ])
pdf("output/Plots/02_CREvarScat.pdf", height = 4.5, width = 4)
p4
dev.off()

#differential analysis
diff_ls <- list(GNP = list("GNP", c("PNC", "Tumor")), PNC = list("PNC", c("GNP", "Tumor")), Tumor = list("Tumor", c("GNP", "PNC")))
DiffTest <- lapply(diff_ls, function(x) edgeR_pairwise(se, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))

DAR_GNP  <- rownames(DiffTest[[1]])[which(assay(DiffTest[[1]])[,"log2FoldChange"] > 1 & assay(DiffTest[[1]])[,"FDR"] < 0.01)]
DAR_PNC  <- rownames(DiffTest[[2]])[which(assay(DiffTest[[2]])[,"log2FoldChange"] > 1 & assay(DiffTest[[2]])[,"FDR"] < 0.01)]
DAR_PT   <- rownames(DiffTest[[1]])[which(assay(DiffTest[[1]])[,"log2FoldChange"] < -1 & assay(DiffTest[[1]])[,"FDR"] < 0.01)]
DAR_Tumor<- rownames(DiffTest[[3]])[which(assay(DiffTest[[3]])[,"log2FoldChange"] > 1 & assay(DiffTest[[3]])[,"FDR"] < 0.01)]

DAR_PNC <- DAR_PNC[DAR_PNC %ni% c(DAR_PT, DAR_GNP, DAR_Tumor)]
DAR_Tumor <- DAR_Tumor[DAR_Tumor %ni% c(DAR_PT, DAR_GNP, DAR_PNC)]

#check exclusiveness
# gplots::venn(list(C = DAR_common, GNP = DAR_GNP, PNC = DAR_PNC, Tumor = DAR_Tumor, PT = DAR_PT))
# rowRanges(se)
# d <- data.frame(seqnames = seqnames(rowRanges(se)), start = start(rowRanges(se))-1, end = end(rowRanges(se)), name = names(rowRanges(se)))
# write.table(d, "output/output_bed/ATAC_allpeaks.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#export
peaks <- rowRanges(se)
gr <- GRangesList(common = peaks[DAR_common], 
                  DAR_GNP = peaks[DAR_GNP], 
                  DAR_PNC = peaks[DAR_PNC],
                  DAR_PT = peaks[DAR_PT], 
                  DAR_Tumor = peaks[DAR_Tumor])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/DAR/", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#DAR heatmap
col_fun2 <- colorRamp2(c(0,1,2,3,4,5), viridis::viridis(6))
tgtDARs <- list(GNP = DAR_GNP, PNC = DAR_PNC, PT = DAR_PT, Tumor = DAR_Tumor)
mtx <- assays(se)$normcounts[unlist(tgtDARs), ]
ht3 <- Heatmap(mtx, 
               name = "Normalized ATAC counts", col = col_fun2, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = F, show_row_dend = F, 
               #row_names_gp = gpar(fontsize = 10), 
               row_title_gp = gpar(fontsize = 8),
               row_split = unlist(lapply(tgtDARs, length)) %>% rep(names(.), .),
               column_title = "", column_title_side = "top", row_title = "", use_raster = T)
mtx_z <- t(scale(t(mtx)))
col_fun3 <- colorRamp2(c(-2,-1,0,1,2), viridis::viridis(5))
ht4 <- Heatmap(mtx_z, 
               name = "z-score[Normalized ATAC counts]", col = col_fun3, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = F, show_row_dend = F, 
               #row_names_gp = gpar(fontsize = 10), 
               row_title_gp = gpar(fontsize = 8),
               row_split = unlist(lapply(tgtDARs, length)) %>% rep(names(.), .),
               column_title = "", column_title_side = "top", row_title = "", use_raster = T)
p5 <- draw(ht3)
p6 <- draw(ht4)

pdf("output/Plots/02_HeatmapDAR.pdf", width = 4, height = 5)
p5
p6
dev.off()
