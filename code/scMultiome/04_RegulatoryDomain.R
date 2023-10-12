#----------------------------------------------------------------------------
# 04_RegulatoryDomain.R
#----------------------------------------------------------------------------
library(ArchR)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
arc <- readRDS("rds/arc.rds")
seu <- readRDS("rds/seu.rds")

#--------- Identifying Regulatory Domains (RD) ---------#
#DORC: Domains of Reguratory chromatin
#Ma et al, Cell, 2020 (SHARE-seq paper)

#read RNA count matrix
inputRNA <- c("data/032103_filtered_feature_bc_matrix.h5",
              "data/RS03056_filtered_feature_bc_matrix.h5",
              "data/RS03060_filtered_feature_bc_matrix.h5")
seRNA <- import10xFeatureMatrix(input = inputRNA, names = c("032103_GNP", "RS03056_H", "RS03060_T"))

#add gene expression to ArchR project
arc <- addGeneExpressionMatrix(input = arc, seRNA = seRNA)

#--------- 1st round: p2g association within +-50kb ---------#
arc <- addPeak2GeneLinks(arc, 
                         reducedDims = "IterativeLSI", 
                         useMatrix = "GeneExpressionMatrix", 
                         maxDist = 50000)
p2g_50k <- getPeak2GeneLinks(arc, corCutOff = 0, FDRCutOff = 0.1, resolution = 500, returnLoops = F) #49357 peak-gene associations

#associated peaks per gene
p2g_50k_gene <- lapply(seq_along(metadata(p2g_50k)$geneSet$name),
                       function(x) metadata(p2g_50k)$peakSet[p2g_50k$idxATAC[which(p2g_50k$idxRNA == x)]])
names(p2g_50k_gene) <- metadata(p2g_50k)$geneSet$name

#number of associated peaks per gene
p2g_50k_gene_no <- lapply(p2g_50k_gene, length) %>% unlist(.)
p2g_50k_gene_no <- sort(p2g_50k_gene_no, decreasing = T) %>% data.frame(gene = names(.), peaks = .)
p2g_50k_gene_no$rank <- c(1:nrow(p2g_50k_gene_no))
p2g_50k_gene_no$rank2 <- c(nrow(p2g_50k_gene_no):1)

p2 <- ggplot(p2g_50k_gene_no, aes(x = rank2, y = peaks)) + geom_point_rast(color = "black", size = 1) + ArchR::theme_ArchR() + 
      labs(x = "Rank sorted genes", y = "Number of correlated peaks") + ggtitle("Peak-gene associations")
write.csv(p2g_50k_gene_no, "output/Tables/04_p2g_50k_gene_no.csv")

#histogram
df <- data.frame(genes = as.numeric(table(p2g_50k_gene_no$peaks)), peaks = names(table(p2g_50k_gene_no$peaks)))
df[12,1] <- sum(df$genes[c(12:nrow(df))])
df <- df[c(2:12),]
df$peaks <- c(c(1:10), "≥11")
p3 <- ggplot(df, aes(x = factor(peaks, levels = c(c(1:10), "≥11")), y = genes)) + 
  geom_bar(stat = "identity", fill = "darkgray") + ArchR::theme_ArchR() + ylim(0, 5200) + 
  labs(x = "Correlated peaks per gene", y = "Number of genes") + ggtitle("Peak-gene associations")

pdf("output/Plots/04_RD_50kb_RankPlot.pdf", width = 6, height = 5)
p2
dev.off()
pdf("output/Plots/04_RD_50kb_Histo.pdf", width = 5, height = 5)
p3
dev.off()

#--------- 2nd round: p2g association within +-500kb ---------#
#re-calculate association for only genes with 5 peaks within +-50kb
DORC_genes <- p2g_50k_gene_no$gene[which(p2g_50k_gene_no$peaks >= 10)]

arc <- addPeak2GeneLinks(arc,
                         reducedDims = "IterativeLSI", 
                         useMatrix = "GeneExpressionMatrix", 
                         maxDist = 500000)
p2g_500k <- getPeak2GeneLinks(arc, corCutOff = 0, FDRCutOff = 0.1, resolution = 500, returnLoops = F) #328681 peak-gene associations

DORC_p2g <- lapply(DORC_genes, function(x){
  idx <- p2g_500k$idxATAC[p2g_500k$idxRNA == which(metadata(p2g_500k)$geneSet$name == x)]
  gr <- metadata(p2g_500k)$peakSet[idx]
  return(gr)
})
DORC_p2g <- GRangesList(DORC_p2g) %>% `names<-`(., DORC_genes)

#----- save objects -----#
saveRDS(arc, "rds/arc.rds")
saveRDS(DORC_p2g, "rds/04_DORC_p2g.rds")

lapply(names(DORC_p2g), function(x){
  g <- DORC_p2g[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/DORC/", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#--------- DORC score ---------#
sample2 <- c("RS03056_H" = "PNC", "RS03060_T" = "Tumor", "032103_GNP" = "GNP")
arc$Sample2 <- as.character(sample2[arc$Sample])
arc$cellType2 <- arc$cellType
idx <- which(arc$cellType %in% c("ProliferativeCells", "DifferentiatedCells"))
arc$cellType2[idx] <- paste0(arc$Sample2[idx], "_", arc$cellType[idx])

PeakMatrix <- getMatrixFromProject(arc, useMatrix = "PeakMatrix")
PeakMatrix <- PeakMatrix[, colnames(seu)]
assays(PeakMatrix)$NormCount <- assays(PeakMatrix)$PeakMatrix / colSums(assays(PeakMatrix)$PeakMatrix)

DORC_mtx <- lapply(DORC_p2g, function(x){
  fo <- findOverlaps(rowRanges(PeakMatrix), x)
  out <- colSums(assays(PeakMatrix)$NormCount[queryHits(fo),])
  return(out)
}) %>% do.call(rbind, .)
rownames(DORC_mtx) <- names(DORC_p2g)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
p4 <- lapply(c("Cntn2","Olig2","Barhl1","Ccnd2","Gli1","Pax2","Atoh1","Mcm2","Ptch1","Neurod1","Smo","Bcor"), function(x){
  tmp <- DORC_mtx[x,]
  tmp <- range01(tmp)
  seu$DORCscore <- tmp
  out <- FeaturePlot(seu, features = "DORCscore", order = T) + scale_colour_gradientn(colours = ArchRPalettes$comet) + 
    ArchR::theme_ArchR() + ggtitle(x)
})
#p5 <- lapply(rownames(DORC_mtx), function(x) FeaturePlot(seu, features = x, order = T, slot = "data") + scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR())

pdf("output/Plots/04_DORCscore_UMAPoverlay.pdf", width = 4.5, height = 5)
p4
dev.off()
# pdf("output/Plots/04_DORC_GEx_UMAPoverlay.pdf", width = 4.5, height = 5)
# p5
# dev.off()

#per celltype
DORC_mtx_cluster <- lapply(unique(arc$cellType2), function(x){
  idy <- arc$cellNames[which(arc$cellType2 == x)]
  out <- rowMeans(DORC_mtx[,idy])
  return(out)
}) %>% do.call(rbind, .)
rownames(DORC_mtx_cluster) <- unique(arc$cellType2)

DORC_mtx_cluster_scaled <- apply(DORC_mtx_cluster, 2, range01)
DORC_mtx_cluster_scaled_var <- DORC_mtx_cluster_scaled[,order(colVars(DORC_mtx_cluster_scaled), decreasing = T)[c(1:100)]]

col_fun1 <- colorRamp2(c(0,0.25,0.5,0.75,1), c("white",ArchRPalettes$comet))
fh <- function(x) hclust(dist(x), method = "ward.D2")

ht1 <- Heatmap(DORC_mtx_cluster_scaled_var, 
               cluster_rows = fh, cluster_columns = fh,
               col = col_fun1, show_column_names = T, column_names_gp = gpar(fontsize = 6))
p6 <- draw(ht1)
pdf("output/Plots/04_Heatmap_DORCscore_v1.pdf", width = 11, height = 4)
p6
dev.off()

DORC_mtx_cluster_scaled_var2 <- DORC_mtx_cluster_scaled[,order(colVars(DORC_mtx_cluster_scaled), decreasing = T)[c(1:250)]]
ht2 <- Heatmap(DORC_mtx_cluster_scaled_var2, 
               cluster_rows = fh, cluster_columns = fh,
               col = col_fun1, show_column_names = T, column_names_gp = gpar(fontsize = 5))
p7 <- draw(ht2)
pdf("output/Plots/04_Heatmap_DORCscore_v2.pdf", width = 18, height = 3.5)
p7
dev.off()

idx <- c(paste0(c("GNP","PNC","Tumor"),"_ProliferativeCells"),paste0(c("GNP","PNC","Tumor"),"_DifferentiatedCells"),"InterneuronProgenitors","Astrocytes","Microglia","Oligodendrocytes","VascularFibroblasts", "Endothelial")
DORC_mtx_cluster_scaled_goi <- DORC_mtx_cluster_scaled[idx, c("Cntn2","Olig2","Barhl1","Ccnd2","Gli1","Pax2","Atoh1","Mcm2","Ptch1","Neurod1","Smo","Bcor")]
ht3 <- Heatmap(DORC_mtx_cluster_scaled_goi, 
               cluster_rows = F, cluster_columns = fh,
               col = col_fun1, show_column_names = T, column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8))
p8 <- draw(ht3)
pdf("output/Plots/04_Heatmap_DORCscore_v3.pdf", width = 6, height = 3.5)
p8
dev.off()

#--------- GenomeTrack plot ---------#
p2g_500k_gr <- getPeak2GeneLinks(arc, corCutOff = 0, FDRCutOff = 0.1, resolution = 500, returnLoops = T) #309839 peak-gene associations

start_gr <- resize(p2g_500k_gr$Peak2GeneLinks, width = 1, fix = "start")
end_gr <- resize(p2g_500k_gr$Peak2GeneLinks, width = 1, fix = "end")

# Atoh1
atoh1_gr <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == "Atoh1")]

n1 <- findOverlaps(start_gr, atoh1_gr) %>% queryHits()
n2 <- findOverlaps(end_gr, atoh1_gr) %>% queryHits()
gr <- p2g_500k_gr$Peak2GeneLinks[unique(sort(c(n1,n2)))]
gr$FDR2 <- -log10(gr$FDR)
gr$value <- gr$FDR2

trackcolors <- c("#FFB300","#803E75","#FF6800") %>% `names<-`(., paste0(c("GNP", "PNC", "Tumor"), "_ProliferativeCells"))

p9 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "cellType2", useGroups = paste0(c("GNP", "PNC", "Tumor"), "_ProliferativeCells"),
  geneSymbol = "Atoh1", 
  pal = trackcolors,
  region = GRanges(seqnames = "chr6", ranges = IRanges(start = 64229001, end = 65190000)),
  loops = gr)
grid::grid.newpage()
grid::grid.draw(p9)

pdf("output/Plots/04_GenomeTrack_Atoh1.pdf", width = 10, height = 8)
grid::grid.draw(p9)
dev.off()

# Smo
smo_gr <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == "Smo")]
smo_gr <- resize(smo_gr, width = 500, fix = "start")

n1 <- findOverlaps(start_gr, smo_gr) %>% queryHits()
n2 <- findOverlaps(end_gr, smo_gr) %>% queryHits()
gr <- p2g_500k_gr$Peak2GeneLinks[unique(sort(c(n1,n2)))]
gr$FDR2 <- -log10(gr$FDR)
gr$value <- gr$FDR2

trackcolors <- c("#FFB300","#803E75","#FF6800") %>% `names<-`(., paste0(c("GNP", "PNC", "Tumor"), "_ProliferativeCells"))

p10 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "cellType2", useGroups = paste0(c("GNP", "PNC", "Tumor"), "_ProliferativeCells"),
  geneSymbol = "Smo", 
  pal = trackcolors,
  region = GRanges(seqnames = "chr6", ranges = IRanges(start = 29235001, end = 30261000)),
  loops = gr)
grid::grid.newpage()
grid::grid.draw(p10)

pdf("output/Plots/04_GenomeTrack_Smo.pdf", width = 10, height = 8)
grid::grid.draw(p10)
dev.off()

p11 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "cellType2", useGroups = paste0(c("GNP", "PNC", "Tumor"), "_ProliferativeCells"),
  geneSymbol = "Smo", 
  pal = trackcolors,
  region = GRanges(seqnames = "chr6", ranges = IRanges(start = 29230001, end = 30200000)),
  loops = gr)
grid::grid.newpage()
grid::grid.draw(p11)

pdf("output/Plots/04_GenomeTrack_Smo_v2.pdf", width = 10, height = 8)
grid::grid.draw(p11)
dev.off()
