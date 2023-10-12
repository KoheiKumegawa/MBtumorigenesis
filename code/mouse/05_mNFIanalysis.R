#------------------------------------------------------------------------------
# 05_mNFIanalysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
source("code/edgeR_PairwiseFunction.R")
se <- readRDS("rds/mNFI_se.rds")

#--------- Annotate Peaks ---------#
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(rowRanges(se), tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Mm.eg.db")
mcols(se) <- mcols(as.GRanges(peakAnno))

pdf("output/Plots/05_peakAnnoPie_NFI.pdf", width = 7, height = 4)
plotAnnoPie(peakAnno)
dev.off()

#--------- overlaps with public data ---------#
d <- data.frame(seqnames = seqnames(rowRanges(se)), start = start(rowRanges(se))-1, end = end(rowRanges(se)))
write.table(d, "output/output_bed/DBR/NFI_consensusPeak.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#By bedtools fisher
#sort -k1,1 -k2,2n my_file.bed > my_file.sorted.bed
#sort -k1,1 hg38.chrom.sizes > hg38.chrom.sorted.sizes
c(phyper(14170 - 1, 99387, 3337679 - 99387, 16993, lower.tail=F, log.p = T)/2.303, #GSM4407221
  phyper(18897 - 1, 99387, 3314912 - 99387, 23541, lower.tail=F, log.p = T)/2.303, #GSM4407222
  phyper(28367 - 1, 99387, 3121709 - 99387, 36429, lower.tail=F, log.p = T)/2.303, #GSM4407223
  phyper(25975 - 1, 99387, 3130419 - 99387, 33283, lower.tail=F, log.p = T)/2.303) #GSM4407224

filenames <- list.files("data/GSE146793_RAW/", pattern = ".narrowPeak")
gr <- lapply(filenames, function(i) rtracklayer::import(paste0("data/GSE146793_RAW/",i)))

n_overlaps <- lapply(gr, function(x) length(which(countOverlaps(x, rowRanges(se)) != 0))) %>% unlist
df <- data.frame(data = filenames, Overlaps = n_overlaps, Length = unlist(lapply(gr, length)))
df$Unique <- df$Length - df$Overlaps
df <- df[,-3]
df <- reshape2::melt(df)

p1 <- ggplot(df, aes(x = data, y = value, fill = factor(variable, levels = c("Unique","Overlaps")))) + 
      geom_bar(stat = "identity") + scale_fill_manual(values = c("darkgray", "orange")) + ArchR::theme_ArchR() + 
      labs(x = "", y = "Number of peaks", fill = "Overlap") + theme(axis.text.x = element_text(angle = 90, size = 6))
p1_2 <- ggplot(df, aes(x = data, y = value, fill = factor(variable, levels = c("Unique","Overlaps")))) + 
  geom_bar(stat = "identity", position = "fill") + scale_fill_manual(values = c("darkgray", "orange")) + ArchR::theme_ArchR() + 
  labs(x = "", y = "fraction of peaks", fill = "Overlap") + theme(axis.text.x = element_text(angle = 90, size = 6))

pdf("output/Plots/05_overlaps_NFI_published.pdf", height = 8, width = 4)
p1
p1_2
dev.off()

#--------- correlation between NFIA and NFIB  ---------#
df <- data.frame(GNP_NFIA = rowMeans(assays(se)$log2cpm[, which(se$sampleType == "GNP" & se$ab == "NFIA")]),
                 GNP_NFIB = rowMeans(assays(se)$log2cpm[, which(se$sampleType == "GNP" & se$ab == "NFIB")]),
                 PNC_NFIA = rowMeans(assays(se)$log2cpm[, which(se$sampleType == "PNC" & se$ab == "NFIA")]),
                 PNC_NFIB = rowMeans(assays(se)$log2cpm[, which(se$sampleType == "PNC" & se$ab == "NFIB")]),
                 Tumor_NFIA = rowMeans(assays(se)$log2cpm[, which(se$sampleType == "Tumor" & se$ab == "NFIA")]),
                 Tumor_NFIB = rowMeans(assays(se)$log2cpm[, which(se$sampleType == "Tumor" & se$ab == "NFIB")]))
idx <- list(GNP = c(1,2), PNC = c(3,4), Tumor = c(5,6))
colors <- c(GNP = "blue", PNC = "orange", Tumor = "red")
p2 <- lapply(names(idx), function(x){
  mtx <- df[, idx[[x]]]
  colnames(mtx) <- c("NFIA", "NFIB")
  cor(mtx)
  out <- ggplot(mtx, aes(x = NFIA, y = NFIB)) + ggrastr::geom_point_rast(size = 0.5, color = colors[x]) + ArchR::theme_ArchR() + 
         ggtitle(x) + geom_smooth(method = "lm", se = F, lty = "dashed", color = "black")
  return(out)
})
pdf("output/Plots/05_NFIAB_Scat.pdf", height = 4, width = 4)
p2
dev.off()
cor.test(df$GNP_NFIA, df$GNP_NFIB) #0.9458515, p-value < 2.2e-16
cor.test(df$PNC_NFIA, df$PNC_NFIB) #0.9137582, p-value < 2.2e-16
cor.test(df$Tumor_NFIA, df$Tumor_NFIB) #0.9172935, p-value < 2.2e-16

#--------- Identify differential accessible regions ---------#
#commonly accessible
stat1 <- data.frame(row.names = mcols(se)$name, median = rowMedians(assays(se)$log2cpm), variance = rowVars(assays(se)$log2cpm))
p3 <- ggplot(stat1, aes(x = variance, y = median)) + geom_hex(bins = 100) + ArchR::theme_ArchR() + 
      scale_fill_viridis_c(trans = "log", option = "E") + labs(x = "Variance", y = "Median accessibility") + xlim(0,3) +
      geom_hline(yintercept = 5, color = "red", lty = "dashed") + geom_vline(xintercept = 0.1, color = "red", lty = "dashed") 
DBR_common <- rownames(stat1[which(stat1$variance < 0.1 & stat1$median > 5), ])
pdf("output/Plots/05_NFIvarScat.pdf", height = 4.5, width = 4)
p3
dev.off()

#differential analysis
diff_ls <- list(GNP = list("GNP", c("PNC", "Tumor")), PNC = list("PNC", c("GNP", "Tumor")), Tumor = list("Tumor", c("GNP", "PNC")))
DiffTest <- lapply(diff_ls, function(x) edgeR_pairwise(se, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))

DBR_GNP  <- rownames(DiffTest[[1]])[which(assay(DiffTest[[1]])[,"log2FoldChange"] > 1 & assay(DiffTest[[1]])[,"FDR"] < 0.05)]
DBR_PNC  <- rownames(DiffTest[[2]])[which(assay(DiffTest[[2]])[,"log2FoldChange"] > 1 & assay(DiffTest[[2]])[,"FDR"] < 0.05)]
DBR_PT   <- rownames(DiffTest[[1]])[which(assay(DiffTest[[1]])[,"log2FoldChange"] < -1 & assay(DiffTest[[1]])[,"FDR"] < 0.05)]
DBR_Tumor<- rownames(DiffTest[[3]])[which(assay(DiffTest[[3]])[,"log2FoldChange"] > 1 & assay(DiffTest[[3]])[,"FDR"] < 0.05)]

DBR_PNC <- DBR_PNC[DBR_PNC %ni% c(DBR_PT, DBR_GNP, DBR_Tumor)]
DBR_Tumor <- DBR_Tumor[DBR_Tumor %ni% c(DBR_PT, DBR_GNP, DBR_PNC)]

#check exclusiveness
gplots::venn(list(C = DBR_common, GNP = DBR_GNP, PNC = DBR_PNC, Tumor = DBR_Tumor, PT = DBR_PT))

#export
peaks <- rowRanges(se)
gr <- GRangesList(common = peaks[DBR_common], 
                  DBR_GNP = peaks[DBR_GNP], 
                  DBR_PNC = peaks[DBR_PNC],
                  DBR_PT = peaks[DBR_PT], 
                  DBR_Tumor = peaks[DBR_Tumor])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/DBR/", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#DBR heatmap
col_fun2 <- colorRamp2(c(0,1,2,3,4,5,6), viridis::viridis(7, option = "E"))
tgtDBRs <- list(GNP = DBR_GNP, PNC = DBR_PNC, PT = DBR_PT, Tumor = DBR_Tumor)
mtx <- assays(se)$log2cpm[unlist(tgtDBRs), ]
ht1 <- Heatmap(mtx, 
               name = "log2CPM", col = col_fun2, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = F, show_row_dend = F, 
               #row_names_gp = gpar(fontsize = 10), 
               row_title_gp = gpar(fontsize = 8),
               row_split = unlist(lapply(tgtDBRs, length)) %>% rep(names(.), .),
               column_title = "", column_title_side = "top", row_title = "", use_raster = T)
p4 <- draw(ht1)
mtx_z <- t(scale(t(mtx)))
col_fun3 <- colorRamp2(c(-2,-1,0,1,2), viridis::viridis(5, option = "E"))
ht2 <- Heatmap(mtx_z, 
               name = "log2CPM", col = col_fun3, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = F, show_row_dend = F, 
               #row_names_gp = gpar(fontsize = 10), 
               row_title_gp = gpar(fontsize = 8),
               row_split = unlist(lapply(tgtDBRs, length)) %>% rep(names(.), .),
               column_title = "", column_title_side = "top", row_title = "", use_raster = T)
p5 <- draw(ht2)

pdf("output/Plots/05_HeatmapDBR.pdf", width = 4, height = 5)
p4
p5
dev.off()

#--------- Overlap with ATAC ---------#
#number of overlaps
filenames <- grep("_sorted", list.files("output/output_bed/DAR/"), value = T, invert = T)
ATAC_gr <- lapply(filenames, function(i) rtracklayer::import.bed(paste0("output/output_bed/DAR/", i)))
filenames <- grep("_sorted", list.files("output/output_bed/DBR/"), value = T, invert = T)[-6]
NFI_gr <- lapply(filenames, function(i) rtracklayer::import.bed(paste0("output/output_bed/DBR/", i)))

names(ATAC_gr) <- c("common", "GNP", "PNC", "PT", "Tumor")
names(NFI_gr) <- c("common", "GNP", "PNC", "PT", "Tumor")

ATAC_NFI_overlaps <- lapply(c("common", "GNP", "PNC", "PT", "Tumor"), function(x){
  gr1 <- ATAC_gr[[x]]
  out2 <- lapply(c("common", "GNP", "PNC", "PT", "Tumor"), function(y){
    gr2 <- NFI_gr[[y]]
    fo <- findOverlaps(gr1, gr2)
    out1 <- unique(queryHits(fo)) %>% length
    return(out1)
  })
  out2 <- unlist(out2)
  out2 <- (out2/length(gr1))*100
  return(out2)
}) %>% do.call(rbind, .)
colnames(ATAC_NFI_overlaps) <- c("common", "GNP", "PNC", "PT", "Tumor")
rownames(ATAC_NFI_overlaps) <- c("common", "GNP", "PNC", "PT", "Tumor")

#overlap heatmap
col_fun4 <- colorRamp2(c(0,10,20,30), c("white", rev(viridis::viridis(3, option = "G"))))
ht3 <- Heatmap(ATAC_NFI_overlaps,
               name = "Percent of overlaps", col = col_fun4, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = T, show_row_dend = F, 
               row_title_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 10),
               column_title = "NFI bound", column_title_side = "top", row_title = "ATAC open", use_raster = T)
p6 <- draw(ht3)

#overlap significance
# sort -k1,1 -k2,2n output/output_bed/DAR/common.bed > output/output_bed/DAR/common_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DAR/DAR_GNP.bed > output/output_bed/DAR/DAR_GNP_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DAR/DAR_PNC.bed > output/output_bed/DAR/DAR_PNC_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DAR/DAR_PT.bed > output/output_bed/DAR/DAR_PT_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DAR/DAR_Tumor.bed > output/output_bed/DAR/DAR_Tumor_sorted.bed
# 
# sort -k1,1 -k2,2n output/output_bed/DBR/common.bed > output/output_bed/DBR/common_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DBR/DBR_GNP.bed > output/output_bed/DBR/DBR_GNP_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DBR/DBR_PNC.bed > output/output_bed/DBR/DBR_PNC_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DBR/DBR_PT.bed > output/output_bed/DBR/DBR_PT_sorted.bed
# sort -k1,1 -k2,2n output/output_bed/DBR/DBR_Tumor.bed > output/output_bed/DBR/DBR_Tumor_sorted.bed

lapply(c("common", "DAR_GNP", "DAR_PNC", "DAR_PT", "DAR_Tumor"), function(i){
  x <- paste0("output/output_bed/DAR/", i, "_sorted.bed")
  res <- lapply(c("common", "DBR_GNP", "DBR_PNC", "DBR_PT", "DBR_Tumor"), function(n){
    y <- paste0("output/output_bed/DBR/", i=n, "_sorted.bed")
    out <- sprintf("bedtools fisher -a %s -b %s -g ref/mm10.chrom.sizes_sorted.txt", x, y)
    return(out)
  })
  return(unlist(res))
})

# bedtools fisher -a output/output_bed/DAR/common_sorted.bed -b output/output_bed/DBR/common_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/common_sorted.bed -b output/output_bed/DBR/DBR_GNP_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/common_sorted.bed -b output/output_bed/DBR/DBR_PNC_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/common_sorted.bed -b output/output_bed/DBR/DBR_PT_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/common_sorted.bed -b output/output_bed/DBR/DBR_Tumor_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt

# bedtools fisher -a output/output_bed/DAR/DAR_GNP_sorted.bed -b output/output_bed/DBR/common_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt 
# bedtools fisher -a output/output_bed/DAR/DAR_GNP_sorted.bed -b output/output_bed/DBR/DBR_GNP_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_GNP_sorted.bed -b output/output_bed/DBR/DBR_PNC_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt 
# bedtools fisher -a output/output_bed/DAR/DAR_GNP_sorted.bed -b output/output_bed/DBR/DBR_PT_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_GNP_sorted.bed -b output/output_bed/DBR/DBR_Tumor_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt

# bedtools fisher -a output/output_bed/DAR/DAR_PNC_sorted.bed -b output/output_bed/DBR/common_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt 
# bedtools fisher -a output/output_bed/DAR/DAR_PNC_sorted.bed -b output/output_bed/DBR/DBR_GNP_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PNC_sorted.bed -b output/output_bed/DBR/DBR_PNC_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PNC_sorted.bed -b output/output_bed/DBR/DBR_PT_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PNC_sorted.bed -b output/output_bed/DBR/DBR_Tumor_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt

# bedtools fisher -a output/output_bed/DAR/DAR_PT_sorted.bed -b output/output_bed/DBR/common_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PT_sorted.bed -b output/output_bed/DBR/DBR_GNP_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PT_sorted.bed -b output/output_bed/DBR/DBR_PNC_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PT_sorted.bed -b output/output_bed/DBR/DBR_PT_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_PT_sorted.bed -b output/output_bed/DBR/DBR_Tumor_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt

# bedtools fisher -a output/output_bed/DAR/DAR_Tumor_sorted.bed -b output/output_bed/DBR/common_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_Tumor_sorted.bed -b output/output_bed/DBR/DBR_GNP_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_Tumor_sorted.bed -b output/output_bed/DBR/DBR_PNC_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_Tumor_sorted.bed -b output/output_bed/DBR/DBR_PT_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt
# bedtools fisher -a output/output_bed/DAR/DAR_Tumor_sorted.bed -b output/output_bed/DBR/DBR_Tumor_sorted.bed -g ref/mm10.chrom.sizes_sorted.txt


## DAR >> DBR
overlapSig <- 
  -c(phyper(1470 - 1, 11718, 2719991 - 11718, 1503, lower.tail=F, log.p = T)/2.303,
     phyper(52 - 1, 11718, 2719991 - 11718, 3464, lower.tail=F, log.p = T)/2.303,
     phyper(15 - 1, 11718, 2719991 - 11718, 688, lower.tail=F, log.p = T)/2.303,
     phyper(126 - 1, 11718, 2719991 - 11718, 3662, lower.tail=F, log.p = T)/2.303,
     phyper(85 - 1, 11718, 2719991 - 11718, 1235, lower.tail=F, log.p = T)/2.303,
     
     phyper(3 - 1, 4986, 2719991 - 4986, 1503, lower.tail=F, log.p = T)/2.303,
     phyper(1651 - 1, 4986, 2719991 - 4986, 3464, lower.tail=F, log.p = T)/2.303,
     phyper(13 - 1, 4986, 2719991 - 4986, 688, lower.tail=F, log.p = T)/2.303,
     phyper(0 - 1, 4986, 2719991 - 4986, 3662, lower.tail=F, log.p = T)/2.303,
     phyper(14 - 1, 4986, 2719991 - 4986, 1235, lower.tail=F, log.p = T)/2.303,
     
     phyper(1 - 1, 1085, 2719991 - 1085, 1503, lower.tail=F, log.p = T)/2.303,
     phyper(5 - 1, 1085, 2719991 - 1085, 3464, lower.tail=F, log.p = T)/2.303,
     phyper(75 - 1, 1085, 2719991 - 1085, 688, lower.tail=F, log.p = T)/2.303,
     phyper(34 - 1, 1085, 2719991 - 1085, 3661, lower.tail=F, log.p = T)/2.303,
     phyper(0 - 1, 1085, 2719991 - 1085, 1235, lower.tail=F, log.p = T)/2.303,
     
     phyper(2 - 1, 8947, 2719991 - 8947, 1503, lower.tail=F, log.p = T)/2.303,
     phyper(3 - 1, 8947, 2719991 - 8947, 3464, lower.tail=F, log.p = T)/2.303,
     phyper(174 - 1, 8947, 2719991 - 8947, 688, lower.tail=F, log.p = T)/2.303,
     phyper(2511 - 1, 8947, 2719991 - 8947, 3662, lower.tail=F, log.p = T)/2.303,
     phyper(319 - 1, 8947, 2719991 - 8947, 1235, lower.tail=F, log.p = T)/2.303,
     
     phyper(2 - 1, 1411, 2719991 - 1411, 1503, lower.tail=F, log.p = T)/2.303,
     phyper(4 - 1, 1411, 2719991 - 1411, 3464, lower.tail=F, log.p = T)/2.303,
     phyper(0 - 1, 1411, 2719991 - 1411, 688, lower.tail=F, log.p = T)/2.303,
     phyper(89 - 1, 1411, 2719991 - 1411, 3661, lower.tail=F, log.p = T)/2.303,
     phyper(254 - 1, 1411, 2719991 - 1411, 1235, lower.tail=F, log.p = T)/2.303)

overlapSig <- matrix(overlapSig, nrow = 5, ncol = 5)
colnames(overlapSig) <- c("common", "GNP", "PNC", "PT", "Tumor")
rownames(overlapSig) <- c("common", "GNP", "PNC", "PT", "Tumor")

#overlap significant heatmap
col_fun5 <- colorRamp2(c(0,200,400,600,800,1000), c("white", rev(viridis::viridis(5, option = "A"))))
ht4 <- Heatmap(overlapSig,
               name = "-log10(P-value)", col = col_fun5, 
               cluster_rows = F, cluster_columns = F,
               show_column_names = T, show_column_dend = F, show_row_names = T, show_row_dend = F, 
               row_title_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 10),
               column_title = "NFI bound", column_title_side = "top", row_title = "ATAC open", use_raster = T)
p7 <- draw(ht4)
pdf("output/Plots/05_HeatmapOverlapSignif_DBRvsDAR.pdf", width = 4, height = 3)
p6
p7
dev.off()

#overlap regions
overlaps_ATAC_NFI_common <- subsetByOverlaps(ATAC_gr$common, NFI_gr$common)
overlaps_ATAC_NFI_GNP <- subsetByOverlaps(ATAC_gr$GNP, NFI_gr$GNP)
overlaps_ATAC_NFI_PNC <- subsetByOverlaps(ATAC_gr$PNC, NFI_gr$PNC)
overlaps_ATAC_NFI_PT <- subsetByOverlaps(ATAC_gr$PT, NFI_gr$PT)
overlaps_ATAC_NFI_Tumor <- subsetByOverlaps(ATAC_gr$Tumor, NFI_gr$Tumor)
overlaps_ATAC_NFI_Common_PT_Tumor <- subsetByOverlaps(ATAC_gr$common, c(NFI_gr$PT, NFI_gr$Tumor))

gr <- GRangesList(overlaps_ATAC_NFI_common = overlaps_ATAC_NFI_common, 
                  overlaps_ATAC_NFI_GNP = overlaps_ATAC_NFI_GNP, 
                  overlaps_ATAC_NFI_PNC = overlaps_ATAC_NFI_PNC,
                  overlaps_ATAC_NFI_PT = overlaps_ATAC_NFI_PT, 
                  overlaps_ATAC_NFI_Tumor = overlaps_ATAC_NFI_Tumor,
                  overlaps_ATAC_NFI_Common_PT_Tumor = overlaps_ATAC_NFI_Common_PT_Tumor)
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/overlaps/", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#--------- 208 commonly open & PT or Tumor specific NFI bound regions ---------#
library(org.Mm.eg.db)
gene <- genes(txdb)
anno <- AnnotationDbi::select(org.Mm.eg.db, keys = gene$gene_id, keytype = "ENTREZID", columns = c("SYMBOL"))
gene$symbol <- anno$SYMBOL

overlaps_ATAC_NFI_Common_PT_Tumor
gr <- ArchR::extendGR(overlaps_ATAC_NFI_Common_PT_Tumor, upstream = 5000, downstream = 5000)
gr <- subsetByOverlaps(gene, gr)
gr$symbol %>% unique %>% sort
write.table(data.frame(gene = gr$symbol %>% unique %>% sort), "output/Tables/05_proxgenes_ATAC_common_NFI_PTorTumor.txt", row.names = F, col.names = F, quote = F)

#--------- 2471 PT specific open & NFI bound regions ---------#
atac_se <- readRDS("rds/mATAC_se.rds")
gr1 <- subsetByOverlaps(rowRanges(atac_se), overlaps_ATAC_NFI_PT)
gr2 <- subsetByOverlaps(rowRanges(atac_se), overlaps_ATAC_NFI_GNP)

peakAnno2 <- annotatePeak(gr1, tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Mm.eg.db")
peakAnno3 <- annotatePeak(gr2, tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Mm.eg.db")

mcols(gr1) <- mcols(as.GRanges(peakAnno2))
mcols(gr2) <- mcols(as.GRanges(peakAnno3))

overlaps_ATAC_NFI_PT_promGenes <- gr1$SYMBOL[which(gr1$annotation == "Promoter")] %>% sort() %>% unique()
overlaps_ATAC_NFI_GNP_promGenes <- gr2$SYMBOL[which(gr2$annotation == "Promoter")] %>% sort() %>% unique()
intersect(overlaps_ATAC_NFI_PT_promGenes, overlaps_ATAC_NFI_GNP_promGenes)

# gr$annotation2 <-  stringr::str_split(gr$annotation, pattern = " ", simplify = T)
# overlaps_ATAC_NFI_PT_promGenes <- gr$SYMBOL[which(gr$annotation2 == "Promoter")] %>% sort() %>% unique()
write.table(data.frame(gene = overlaps_ATAC_NFI_PT_promGenes), 
            "output/Tables/05_overlaps_ATAC_NFI_PT_promGenes.txt", row.names = F, col.names = F, quote = F)
write.table(data.frame(gene = overlaps_ATAC_NFI_GNP_promGenes), 
            "output/Tables/05_overlaps_ATAC_NFI_GNP_promGenes.txt", row.names = F, col.names = F, quote = F)

pdf("output/Plots/05_peakAnnoPie_overlaps_ATAC_NFI_PT.pdf", width = 7, height = 4)
plotAnnoPie(peakAnno2)
dev.off()
pdf("output/Plots/05_peakAnnoPie_overlaps_ATAC_NFI_GNP.pdf", width = 7, height = 4)
plotAnnoPie(peakAnno3)
dev.off()

# not promoter, +-5kb
gr1 <- ArchR::extendGR(overlaps_ATAC_NFI_GNP, upstream = 5000, downstream = 5000)
gr1 <- subsetByOverlaps(gene, gr1)
gr1$symbol %>% unique %>% sort

gr2 <- ArchR::extendGR(overlaps_ATAC_NFI_PT, upstream = 5000, downstream = 5000)
gr2 <- subsetByOverlaps(gene, gr2)
gr2$symbol %>% unique %>% sort

write.table(data.frame(gene = gr1$symbol %>% sort), "output/Tables/05_proxgenes_ATAC_NFI_GNP.txt", row.names = F, col.names = F, quote = F)
write.table(data.frame(gene = gr2$symbol %>% sort), "output/Tables/05_proxgenes_ATAC_NFI_PT.txt", row.names = F, col.names = F, quote = F)

#--------- consensus peaks homer visualize ---------#
df <- data.table::fread("output/homer_motifs/DBR/consensusPeak/knownResults.txt", header = F, skip = 1)[c(1:10), c(1,3)] %>% data.frame %>% `colnames<-`(., c("Motif", "P"))
df$mlog10P <- as.numeric(gsub("1e-", "", df$P))

p8 <- ggplot(df, aes(x = mlog10P, y = reorder(Motif, mlog10P))) + geom_bar(stat = "identity") + ArchR::theme_ArchR() + 
      labs(x = "Motif", y = "-log10(P-value)") + ggtitle("NFI consensus peaks")
pdf("output/Plots/05_HomerMotif_NFI_consensusPeaks.pdf", height = 3, width = 10)
p8
dev.off()

#--------- Annotated DAR & DBR ---------#
# DBR
NFI_gr_anno <- lapply(NFI_gr, function(x) annotatePeak(x, tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Mm.eg.db"))
NFI_gr2 <- lapply(names(NFI_gr), function(x){
  gr <- NFI_gr[[x]]
  mcols(gr) <- mcols(as.GRanges(NFI_gr_anno[[x]]))
  return(gr)
}) %>% `names<-`(.,names(NFI_gr))
df_anno <- lapply(names(NFI_gr2), function(x){
  gr <- NFI_gr2[[x]]
  out <- data.frame(DBR_class = x, seqnames = seqnames(gr), start = start(gr), end = end(gr))
  out <- cbind(out, mcols(gr)) 
  return(out)
}) %>% do.call(rbind, .)
write.csv(df_anno, "output/Tables/05_DBR_annotated.csv", row.names = F, quote = F)

df_gRegion <- lapply(names(NFI_gr_anno), function(x){
  out <- NFI_gr_anno[[x]]@annoStat %>% as.data.frame()
  out$DBR_class <- x
  return(out)
}) %>% do.call(rbind,.)

p9 <- ggplot(df_gRegion, aes(x = DBR_class, y = Frequency, fill = Feature)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("lightblue","deepskyblue3", "green4", "pink", "indianred2", "darkorange1", "plum3", "thistle", "saddlebrown") %>% `names<-`(., df_gRegion$Feature[c(1:9)])) + ArchR::theme_ArchR()
pdf("output/Plots/05_DBR_genomicAnnotation.pdf", height = 4, width = 4)
p9
dev.off()

# DAR
ATAC_gr_anno <- lapply(ATAC_gr, function(x) annotatePeak(x, tssRegion=c(-1000, 100), TxDb = txdb, annoDb = "org.Mm.eg.db"))
ATAC_gr2 <- lapply(names(ATAC_gr), function(x){
  gr <- ATAC_gr[[x]]
  mcols(gr) <- mcols(as.GRanges(ATAC_gr_anno[[x]]))
  return(gr)
}) %>% `names<-`(.,names(ATAC_gr))
df_anno_ATAC <- lapply(names(ATAC_gr2), function(x){
  gr <- ATAC_gr2[[x]]
  out <- data.frame(DAR_class = x, seqnames = seqnames(gr), start = start(gr), end = end(gr))
  out <- cbind(out, mcols(gr)) 
  return(out)
}) %>% do.call(rbind, .)
write.csv(df_anno_ATAC, "output/Tables/05_DAR_annotated.csv", row.names = F, quote = F)
