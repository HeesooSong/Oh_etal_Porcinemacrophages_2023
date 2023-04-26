#### DeSeq2 analysis of RNAseq data 
library(BiocManager)
# BiocManager::install("DESeq2")
library("DESeq2")
library("dplyr")
library("tibble")
library(ggplot2)
library(gplots)
library(ggrepel)
library(RColorBrewer)
library(reshape)
library(cowplot)
library(ggpubr)

########################################
### Loading the data 
########################################
# ## LCM trimmed  (uncomment this if you want to analyse the LCM data)
# dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/"
# vis_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/Visualizations_3groups/"
# f_name <- "Dayoung_LCM_2021_htseq_counts_Trimmed.tsv"
# ##

### FACS trimmed (uncomment this if you want to analyse the FACS data)
dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/"
vis_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/Visualizations_3groups/"
f_name <- "Dayoung_MF_2019_htseq_counts_Trimmed.tsv"
###

cts <- as.matrix(read.csv(paste(dir_path, f_name, sep=""),sep="\t",row.names="X"))
dim(cts)

# #### LCM: Change column (A4 <-> B4) (Uncomment for LCM)
# colnames(cts) <- c("A1", "A2", "A3", "B4", "B1", "B2", "B3", "A4", "C1", "C2", "C3", "C4")
# cts <- cts[,c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")]

# read the metadata 
coldata <- read.csv("C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/coldata.tsv", sep="\t", row.names=1)
coldata <- coldata[,c("condition","type")]


coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

## Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(ENS_NCBI)
ENS_NCBI[0:5,]
ENS_NCBI[is.na(ENS_NCBI)] <- ""
ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$geneName_pig <- ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$Ensembl_GeneID_pig
ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$geneName_human <- ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$Ensembl_GeneID_pig

## Change some gene name
ENS_NCBI[ENS_NCBI$geneName_pig == "ENSSSCG00000012699",c("geneName_pig", "geneName_human")] <- "FHL1"
ENS_NCBI[ENS_NCBI$geneName_pig == "ENSSSCG00000036155",c("geneName_pig", "geneName_human")] <- "FAT4"
ENS_NCBI[ENS_NCBI$geneName_pig == "ENSSSCG00000011355",c("geneName_pig", "geneName_human")] <- "R-SSC-5694530"
ENS_NCBI[ENS_NCBI$geneName_pig == "ENSSSCG00000031851",c("geneName_pig", "geneName_human")] <- "CBR2"

## Deseq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)

dds
dds <- DESeq(dds)
resultsNames(dds)
rld <- rlog( dds)

## PCA
vsd <- vst(dds, blind=FALSE)

# ##### Change condition name
levels(colData(vsd)$condition) <- c("NaSn-", "NaSn+", "LuSn+")
# colData(vsd)[c("A1", "A2", "A3", "A4"), "condition"] <- "NaSn-"
# colData(vsd)[c("B1", "B2", "B3", "B4"), "condition"] <- "NaSn+"
# colData(vsd)[c("C1", "C2", "C3", "C4"), "condition"] <- "LuSn+"

filepath <- paste(vis_path, "PCA.pdf", sep="")
pdf(file = filepath, width = 8, height = 6)
par(mar = c(2,2,2,2))
plotPCA(vsd, intgroup=c("condition")) +
  geom_point(size = 6) +
  xlim(-40, 65) +
  ylim(-30, 35) +
  scale_color_discrete(labels = c(bquote("NaSn"^"-"), bquote("NaSn"^"+"), bquote("LuSn"^"+"))) +
  theme_bw() +
  theme(axis.title=element_text(size=30),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key.size = unit(1.5, 'cm'),
        legend.position = c(0.87, 0.76),
        legend.box.background = element_rect(colour = "black", size = 0.5))
dev.off()

## Calculate normalisationfactors
dds$sizeFactor
read_counts <- cts
read_counts[0:5,]
normalized_read_counts <- read_counts %*% diag(dds$sizeFactor)
normalized_read_counts[0:5,]
colnames(normalized_read_counts) <- colnames(read_counts)
normalized_read_counts[0:5,]

# Read in the manually cleaned out full table with verified gene names
# This is used for volcano plots in this script
clean_table <- as.data.frame(read.csv(paste(dir_path, "files_3_groups/DeSeq2_Full_Table_nozeroexpression_final.csv", sep=""),sep=","))
colnames(clean_table)[1] <- "Ensembl_GeneID_pig"


#########################################################################
## B vs A 
#########################################################################
B_NasalSnPos_vs_A_NasalSnNeg <- results(dds, contrast=c("condition", "B_NasalSnPos", "A_NasalSnNeg"), independentFiltering = FALSE)
dim(B_NasalSnPos_vs_A_NasalSnNeg)
normalized_read_counts_added <- merge(B_NasalSnPos_vs_A_NasalSnNeg, normalized_read_counts, by=0)
write.table(normalized_read_counts_added, paste(dir_path, "files_3_groups/DeSeq2_B_NasalSnPos_vs_A_NasalSnNeg.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
dim(B_NasalSnPos_vs_A_NasalSnNeg)
summary(B_NasalSnPos_vs_A_NasalSnNeg)

B_NasalSnPos_vs_A_NasalSnNeg_filtered <- B_NasalSnPos_vs_A_NasalSnNeg[which(abs(B_NasalSnPos_vs_A_NasalSnNeg$log2FoldChange)>1 & B_NasalSnPos_vs_A_NasalSnNeg$padj < 0.05),]
dim(B_NasalSnPos_vs_A_NasalSnNeg_filtered)
normalized_read_counts_added_filtered <- merge(B_NasalSnPos_vs_A_NasalSnNeg_filtered, normalized_read_counts, by=0)
dim(normalized_read_counts_added_filtered)
write.table(normalized_read_counts_added_filtered, paste(dir_path, "files_3_groups/DeSeq2_B_NasalSnPos_vs_A_NasalSnNeg_filtered.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
B_NasalSnPos_vs_A_NasalSnNeg$expression <- ifelse(B_NasalSnPos_vs_A_NasalSnNeg$padj < 0.05 & abs(B_NasalSnPos_vs_A_NasalSnNeg$log2FoldChange) > 1, 
                                                  ifelse(B_NasalSnPos_vs_A_NasalSnNeg$log2FoldChange> 1 ,"DEGs in NaSn+","DEGs in NaSn-"),
                                                  'non-DEGs')

# Volcano plot
test <- merge(ENS_NCBI, data.frame(B_NasalSnPos_vs_A_NasalSnNeg), by.x = "Ensembl_GeneID_pig", by.y =0)
test <- test[complete.cases(test), ]
test <- test %>%
  rownames_to_column(var="geneID")

## clean out genes (given curated table)
DEgenes <- test[which(abs(test$log2FoldChange) > 1 & test$padj < 0.05), "Ensembl_GeneID_pig"]
clean_geneName <- clean_table[which(clean_table$Ensembl_GeneID_pig %in% DEgenes),c("Ensembl_GeneID_pig", "geneName_pig")]
wrong_geneName <- clean_geneName[which(clean_geneName$geneName_pig == ""),"Ensembl_GeneID_pig"]
wrong_geneName <- c(wrong_geneName, setdiff(DEgenes, clean_geneName$Ensembl_GeneID_pig))
test <- test[!(test$Ensembl_GeneID_pig %in% wrong_geneName),]


# rename genes (given curated table)
DEG <- c()
for (gene in clean_geneName$Ensembl_GeneID_pig){
  geneName_initial <- test$geneName_pig[which(test$Ensembl_GeneID_pig == gene)]
  geneName_curated <- clean_table$geneName_pig[which(clean_table$Ensembl_GeneID_pig == gene)]
  
  if (length(geneName_initial) > 0 & geneName_curated != ""){
    DEG <- c(DEG, geneName_curated)
    if (geneName_initial != geneName_curated){
      test$geneName_pig[which(test$Ensembl_GeneID_pig == gene)] <- geneName_curated
    }
  }
}

## select genes to label
DEG_df <- test[which(abs(test$log2FoldChange) > 1 & test$padj < 0.05), ]
threshold_A <- sort(-log10(DEG_df$padj) * DEG_df$log2FoldChange, decreasing = TRUE)[6]
threshold_B <- sort(-log10(DEG_df$padj) * DEG_df$log2FoldChange, decreasing = FALSE)[6]
volcano_label <- DEG_df[which((-log10(DEG_df$padj)*DEG_df$log2FoldChange) > threshold_A | (-log10(DEG_df$padj)*DEG_df$log2FoldChange) < threshold_B), "Ensembl_GeneID_pig"]

volcano_df <- test
volcano_df$geneName_pig <- ""
volcano_df$geneName_pig[which(volcano_df$Ensembl_GeneID_pig %in% volcano_label)] <- test$geneName_pig[which(volcano_df$Ensembl_GeneID_pig %in% volcano_label)]

volcano_BA <- ggplot(data = volcano_df, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=expression)) +
  geom_point(alpha=0.8, size=2) +
  scale_color_manual(values=c("DEGs in NaSn-"="#F8766D", "DEGs in NaSn+"="#00BA38", "DEGs in LuSn+"="#619CFF", "non-DEGs"="grey"),
                     labels = c(bquote("DEGs in NaSn"^"-"), bquote("DEGs in NaSn"^"+"), bquote("DEGs in LuSn"^"+"), "non-DEGs")) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title = expression("NaSn"^"-"*" vs NaSn"^"+"))  +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5, colour = "#D39200"), 
        axis.title=element_text(size=20),
        axis.text = element_text(size=20),
        legend.position="right", 
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key.size = unit(1, 'cm'))+
  geom_text_repel(
    aes(label=geneName_pig),
    size=7,
    show.legend = FALSE,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 5, max.iter = Inf, force = 100, force_pull = 10,
    box.padding = 0.5,point.padding = 0)


pdf(file = paste(vis_path, "volcanoplot_AvsB.pdf", sep=""),  width = 8, height = 6)
volcano_BA
dev.off()
# 
# # Heatmap
# res = B_NasalSnPos_vs_A_NasalSnNeg[order(B_NasalSnPos_vs_A_NasalSnNeg$padj, na.last=NA), ]
# alpha = 0.05
# sigtab = res[(res$padj < alpha), ]
# sigtab = sigtab[(abs(sigtab$log2FoldChange)>1),]
# test <- merge(ENS_NCBI, as.data.frame(sigtab), by.x = "Ensembl_GeneID_pig", by.y =0)
# rownames(test) <- test$Ensembl_GeneID_pig
# test = test[!(nchar(test$geneName_human) == 18),]
# test = test[order(-test$log2FoldChange),]
# summary(sigtab)
# rownames(sigtab)
# samples <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")
# par(mar=c(1,1,1,1))
# pdf(file = paste(vis_path, "heatmap_AvsB.pdf", sep=""), width = 20, height = 15)
# heatmap.2( assay(rld)[ rownames(test), samples], scale="row", 
#            trace="none", dendrogram="column", 
#            margins = c(10, 12),
#            Rowv = FALSE,
#            #labRow= test[rownames(sigtab),]$geneName_pig,
#            labRow= test[rownames(test),]$geneName_human,
#            #labRow=as.expression(lapply(rownames(sigtab), function(a) bquote(italic(.(a))))),
#            keysize = 0.5,
#            key=TRUE,
#            density.info=c("none"),
#            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(20),
#            labCol = samples,
#            main = "A vs B")
# dev.off()

#########################################################################
## C vs A
#########################################################################
C_LungSnPos_vs_A_NasalSnNeg <- results(dds, contrast=c("condition", "C_LungSnPos", "A_NasalSnNeg"), independentFiltering = FALSE)
dim(C_LungSnPos_vs_A_NasalSnNeg)
normalized_read_counts_added <- merge(C_LungSnPos_vs_A_NasalSnNeg, normalized_read_counts, by=0)
write.table(normalized_read_counts_added, paste(dir_path, "files_3_groups/DeSeq2_C_LungSnPos_vs_A_NasalSnNeg.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
dim(C_LungSnPos_vs_A_NasalSnNeg)
summary(C_LungSnPos_vs_A_NasalSnNeg)
C_LungSnPos_vs_A_NasalSnNeg_filtered <- C_LungSnPos_vs_A_NasalSnNeg[which(abs(C_LungSnPos_vs_A_NasalSnNeg$log2FoldChange)>1 & C_LungSnPos_vs_A_NasalSnNeg$padj < 0.05),]
dim(C_LungSnPos_vs_A_NasalSnNeg_filtered)
normalized_read_counts_added_filtered <- merge(C_LungSnPos_vs_A_NasalSnNeg_filtered, normalized_read_counts, by=0)
dim(normalized_read_counts_added_filtered)
write.table(normalized_read_counts_added_filtered, paste(dir_path, "files_3_groups/DeSeq2_C_LungSnPos_vs_A_NasalSnNeg_filtered.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
C_LungSnPos_vs_A_NasalSnNeg$expression <- ifelse(C_LungSnPos_vs_A_NasalSnNeg$padj < 0.05 & abs(C_LungSnPos_vs_A_NasalSnNeg$log2FoldChange) > 1, 
                                                 ifelse(C_LungSnPos_vs_A_NasalSnNeg$log2FoldChange> 1 ,'DEGs in LuSn+','DEGs in NaSn-'),
                                                 'non-DEGs')

# Volcano Plot
test <- merge(ENS_NCBI, as.data.frame(C_LungSnPos_vs_A_NasalSnNeg), by.x = "Ensembl_GeneID_pig", by.y =0)
test <- test[complete.cases(test), ]
test <- test %>% rownames_to_column(var="geneID")

## clean out genes (given curated table)
DEgenes <- test[which(abs(test$log2FoldChange) > 1 & test$padj < 0.05), "Ensembl_GeneID_pig"]
clean_geneName <- clean_table[which(clean_table$Ensembl_GeneID_pig %in% DEgenes),c("Ensembl_GeneID_pig", "geneName_pig")]
wrong_geneName <- clean_geneName[which(clean_geneName$geneName_pig == ""),"Ensembl_GeneID_pig"]
wrong_geneName <- c(wrong_geneName, setdiff(DEgenes, clean_geneName$Ensembl_GeneID_pig))
test <- test[!(test$Ensembl_GeneID_pig %in% wrong_geneName),]


# rename genes (given curated table)
DEG <- c()
for (gene in clean_geneName$Ensembl_GeneID_pig){
  geneName_initial <- test$geneName_pig[which(test$Ensembl_GeneID_pig == gene)]
  geneName_curated <- clean_table$geneName_pig[which(clean_table$Ensembl_GeneID_pig == gene)]
  
  if (length(geneName_initial) > 0 & geneName_curated != ""){
    DEG <- c(DEG, geneName_curated)
    if (geneName_initial != geneName_curated){
      test$geneName_pig[which(test$Ensembl_GeneID_pig == gene)] <- geneName_curated
    }
  }
}

## select genes to label
DEG_df <- test[which(abs(test$log2FoldChange) > 1 & test$padj < 0.05), ]
threshold_A <- sort(-log10(DEG_df$padj) * DEG_df$log2FoldChange, decreasing = TRUE)[6]
threshold_B <- sort(-log10(DEG_df$padj) * DEG_df$log2FoldChange, decreasing = FALSE)[6]
volcano_label <- DEG_df[which((-log10(DEG_df$padj)*DEG_df$log2FoldChange) > threshold_A | (-log10(DEG_df$padj)*DEG_df$log2FoldChange) < threshold_B), "Ensembl_GeneID_pig"]

volcano_df <- test
volcano_df$geneName_pig <- ""
volcano_df$geneName_pig[which(volcano_df$Ensembl_GeneID_pig %in% volcano_label)] <- test$geneName_pig[which(volcano_df$Ensembl_GeneID_pig %in% volcano_label)]

volcano_CA <- ggplot(data = volcano_df, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=expression)) +
  geom_point(alpha=0.8, size=2) +
  scale_color_manual(values=c("DEGs in NaSn-"="#F8766D", "DEGs in NaSn+"="#00BA38", "DEGs in LuSn+"="#619CFF", "non-DEGs"="grey"),
                     labels = c(bquote("DEGs in NaSn"^"-"), bquote("DEGs in NaSn"^"+"), bquote("DEGs in LuSn"^"+"), "non-DEGs")) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title = expression("LuSn"^"+"*" vs NaSn"^"-"))  +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5, colour = "#DB72FB"), 
        axis.title=element_text(size=20),
        axis.text=element_text(size=20),
        legend.position="right", 
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key.size = unit(1, 'cm'))+
  geom_text_repel(
    aes(label=geneName_pig),
    size=7,
    show.legend = FALSE,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 5, max.iter = Inf, force = 100, force_pull = 10,
    box.padding = 0.5,point.padding = 0.2)

pdf(file = paste(vis_path, "volcanoplot_AvsC.pdf", sep=""), width = 8, height = 6)
volcano_CA
dev.off()

# # Heatmap
# res = C_LungSnPos_vs_A_NasalSnNeg[order(C_LungSnPos_vs_A_NasalSnNeg$padj, na.last=NA), ]
# alpha = 0.05
# sigtab = res[(res$padj < alpha), ]
# sigtab = sigtab[(abs(sigtab$log2FoldChange)>1),]
# sigtab = sigtab[order(-abs(sigtab$log2FoldChange)),]
# # sigtab = head(sigtab, 20)
# sigtab = sigtab[order(-sigtab$log2FoldChange),]
# test <- merge(ENS_NCBI, as.data.frame(sigtab), by.x = "Ensembl_GeneID_pig", by.y =0)
# rownames(test) <- test$Ensembl_GeneID_pig
# summary(sigtab)
# rownames(sigtab)
# samples <- c("A1", "A2", "A3", "A4", "C1", "C2", "C3", "C4")
# par(mar=c(1,1,1,1))
# pdf(file = paste(vis_path, "heatmap_CvsA.pdf", sep=""), width = 20, height = 15)
# heatmap.2( assay(rld)[ rownames(sigtab),samples], scale="row", 
#            trace="none", dendrogram="column", 
#            margins = c(10, 12),
#            Rowv = FALSE,
#            labRow= FALSE,
#            # labRow= test[rownames(sigtab),]$geneName_pig,
#            # labRow= test[rownames(sigtab),]$geneName_human,
#            # labRow=as.expression(lapply(rownames(sigtab), function(a) bquote(italic(.(a))))),
#            keysize = 0.5,
#            key=TRUE,
#            density.info=c("none"),
#            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(20),
#            labCol = samples,
#            main = "C vs A")
# dev.off()

#########################################################################
# C vs B
#########################################################################
C_LungSnPos_vs_B_NasalSnPos <- results(dds, contrast=c("condition","C_LungSnPos","B_NasalSnPos"), independentFiltering = FALSE)
write.table(C_LungSnPos_vs_B_NasalSnPos, paste(dir_path, "files_3_groups/DeSeq2_C_LungSnPos_vs_B_NasalSnPos.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

head(C_LungSnPos_vs_B_NasalSnPos)
C_LungSnPos_vs_B_NasalSnPos_filtered <- C_LungSnPos_vs_B_NasalSnPos[which(abs(C_LungSnPos_vs_B_NasalSnPos$log2FoldChange)>1 & C_LungSnPos_vs_B_NasalSnPos$padj < 0.05),]
dim(C_LungSnPos_vs_B_NasalSnPos_filtered)
write.table(C_LungSnPos_vs_B_NasalSnPos_filtered, paste(dir_path, "files_3_groups/DeSeq2_C_LungSnPos_vs_B_NasalSnPos_filtered.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

C_LungSnPos_vs_B_NasalSnPos$expression <- ifelse(C_LungSnPos_vs_B_NasalSnPos$padj < 0.05 & abs(C_LungSnPos_vs_B_NasalSnPos$log2FoldChange) > 1, 
                                                 ifelse(C_LungSnPos_vs_B_NasalSnPos$log2FoldChange> 1 ,'DEGs in LuSn+','DEGs in NaSn+'),
                                                 'non-DEGs')

# Volcano plot
test <- merge(ENS_NCBI, as.data.frame(C_LungSnPos_vs_B_NasalSnPos), by.x = "Ensembl_GeneID_pig", by.y =0)
test <- test[complete.cases(test), ]
test <- test %>%
  rownames_to_column(var="geneID")

## clean out genes (given curated table)
DEgenes <- test[which(abs(test$log2FoldChange) > 1 & test$padj < 0.05), "Ensembl_GeneID_pig"]
clean_geneName <- clean_table[which(clean_table$Ensembl_GeneID_pig %in% DEgenes),c("Ensembl_GeneID_pig", "geneName_pig")]
wrong_geneName <- clean_geneName[which(clean_geneName$geneName_pig == ""),"Ensembl_GeneID_pig"]
wrong_geneName <- c(wrong_geneName, setdiff(DEgenes, clean_geneName$Ensembl_GeneID_pig))
test <- test[!(test$Ensembl_GeneID_pig %in% wrong_geneName),]


# rename genes (given curated table)
DEG <- c()
for (gene in clean_geneName$Ensembl_GeneID_pig){
  geneName_initial <- test$geneName_pig[which(test$Ensembl_GeneID_pig == gene)]
  geneName_curated <- clean_table$geneName_pig[which(clean_table$Ensembl_GeneID_pig == gene)]
  
  if (length(geneName_initial) > 0 & geneName_curated != ""){
    DEG <- c(DEG, geneName_curated)
    if (geneName_initial != geneName_curated){
      test$geneName_pig[which(test$Ensembl_GeneID_pig == gene)] <- geneName_curated
    }
  }
}

## select genes to label
DEG_df <- test[which(abs(test$log2FoldChange) > 1 & test$padj < 0.05), ]
threshold_A <- sort(-log10(DEG_df$padj) * DEG_df$log2FoldChange, decreasing = TRUE)[6]
threshold_B <- sort(-log10(DEG_df$padj) * DEG_df$log2FoldChange, decreasing = FALSE)[6]
volcano_label <- DEG_df[which((-log10(DEG_df$padj)*DEG_df$log2FoldChange) > threshold_A | (-log10(DEG_df$padj)*DEG_df$log2FoldChange) < threshold_B), "Ensembl_GeneID_pig"]

volcano_df <- test
volcano_df$geneName_pig <- ""
volcano_df$geneName_pig[which(volcano_df$Ensembl_GeneID_pig %in% volcano_label)] <- test$geneName_pig[which(volcano_df$Ensembl_GeneID_pig %in% volcano_label)]


volcano_CB <- ggplot(data = volcano_df, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=expression)) +
  geom_point(alpha=0.8, size=2) +
  scale_color_manual(values=c("DEGs in NaSn-"="#F8766D", "DEGs in NaSn+"="#00BA38", "DEGs in LuSn+"="#619CFF", "non-DEGs"="grey"),
                     labels = c(bquote("DEGs in NaSn"^"-"), bquote("DEGs in NaSn"^"+"), bquote("DEGs in LuSn"^"+"), "non-DEGs"))+
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title = expression("LuSn"^"+"*" vs NaSn"^"+"))  +
  theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(plot.title = element_text(size = 30, hjust = 0.5, colour = "#00798C"), 
        axis.title=element_text(size=20),
        axis.text = element_text(size=20),
        legend.position="right", 
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key.size = unit(1, 'cm'))+
  geom_text_repel(
    aes(label=geneName_pig),
    size=7,
    show.legend = FALSE,
    min.segment.length = 0,
    max.time = 5, force = 100, force_pull = 10,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0)
pdf(file = paste(vis_path, "volcanoplot_BvsC.pdf", sep=""), width = 8, height = 6)
volcano_CB
dev.off()


# # Heatmap
# res = C_LungSnPos_vs_B_NasalSnPos[order(C_LungSnPos_vs_B_NasalSnPos$padj, na.last=NA), ]
# alpha = 0.05
# sigtab = res[(res$padj < alpha), ]
# sigtab = sigtab[(abs(sigtab$log2FoldChange)>1),]
# sigtab = sigtab[order(-abs(sigtab$log2FoldChange)),]
# #sigtab = head(sigtab, 20)
# sigtab = sigtab[order(-sigtab$log2FoldChange),]
# test <- merge(ENS_NCBI, as.data.frame(sigtab), by.x = "Ensembl_GeneID_pig", by.y =0)
# rownames(test) <- test$Ensembl_GeneID_pig
# summary(sigtab)
# rownames(sigtab)
# dim(sigtab)
# samples = c("B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")
# pdf(file = paste(vis_path, "heatmap_CvsB.pdf", sep=""), width = 20, height = 15)
# colnames(assay(rld))
# heatmap.2( assay(rld)[ rownames(sigtab),samples], scale="row", 
#            trace="none", dendrogram="column", 
#            margins = c(10, 12),
#            Rowv = FALSE,
#            # Colv = FALSE,
#            labRow = FALSE,
#            #labRow= test[rownames(sigtab),]$geneName_human,
#            #labCol = colnames(assay(rld)),
#            #labRow=as.expression(lapply(rownames(sigtab), function(a) bquote(italic(.(a))))),
#            keysize = 0.5,
#            key=TRUE,
#            density.info=c("none"),
#            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(20),
#            labCol = samples,
#            main = "B vs C")
# dev.off()

#########################################################################
# Combine volcano plots
#########################################################################

pdf(file = paste(vis_path, "volcanoplot_combined.pdf", sep=""), width = 19, height = 6)
ggarrange(volcano_BA, volcano_CA, volcano_CB, ncol = 3, common.legend = TRUE, legend = "right")
dev.off()


# #########################################################################
# # Dot Plot
# #########################################################################
# 
# # Read in marker gene data
# marker_genes_df <- read.csv("C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/cell_type_signature_gene_list_new2.csv", sep = ",")
# colnames(marker_genes_df)[1] <- "Isolation"
# levels(factor(marker_genes_df$Isolation))
# 
# ##Be aware to change LCM / FACS!!
# marker_genes_df <- marker_genes_df[marker_genes_df$Isolation == "FACS",]
# marker_genes_df %>% names()
# 
# marker_genes <- marker_genes_df[,"Ensembl_GeneID_pig"]
# 
# # Prepare data frame for ggplot
# ## select marker genes from normalized count data 
# dotplot_df <- assay(rld)[rownames(assay(rld)) %in% marker_genes,] %>% data.frame()
# 
# ## calculate average expression of each group
# dotplot_df$Avg_norm_A <- rowMeans(dotplot_df[,c("A1", "A2", "A3", "A4")])
# dotplot_df$Avg_norm_B <- rowMeans(dotplot_df[,c("B1", "B2", "B3", "B4")])
# dotplot_df$Avg_norm_C <- rowMeans(dotplot_df[,c("C1", "C2", "C3", "C4")])
# dotplot_df <- dotplot_df[,c("Avg_norm_A", "Avg_norm_B", "Avg_norm_C")]
# 
# ## add metadata column
# #dotplot_df <- merge(dotplot_df, ENS_NCBI[,c("Ensembl_GeneID_pig","geneName_pig")], by.x = 0, by.y = "Ensembl_GeneID_pig", all.x = TRUE) %>% column_to_rownames(var = "Row.names")
# dotplot_df <- merge(dotplot_df, marker_genes_df[,c("Ensembl_GeneID_pig", "geneName_pig", "Cell.type", "Group")], by.x = 0, by.y = "Ensembl_GeneID_pig", all.x = TRUE) %>% column_to_rownames(var = "Row.names")
# 
# 
# ## melt expression data for ggplot
# dotplot_df <- melt(dotplot_df)
# dotplot_df$Cell_type <- factor(dotplot_df$Cell.type)
# dotplot_df <- dotplot_df[order(dotplot_df$Cell.type),]
# Group_lvl <- c("A", "B", "C", "A, B, C")
# dotplot_df$Group <- factor(dotplot_df$Group, levels = Group_lvl)
# dotplot_df <- dotplot_df[order(dotplot_df$Group),]
# dotplot_df$geneName_pig <- factor(dotplot_df$geneName_pig, levels = rev(unique(dotplot_df$geneName_pig)))
# 
# ## Split data frame into macrophage and other cell types
# dotplot_df_macrophages <- dotplot_df[which(dotplot_df$Group == "A, B, C"),]
# dotplot_df_etc <- dotplot_df[!(dotplot_df$Group == "A, B, C"),]
# 
# # ggplot
# 
# scale.func <- switch(
#   EXPR = "radius",
#   'radius' = scale_radius,
# )
# 
# cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# 
# ##### ETC. ########################
# pdf(file = paste(vis_path, "DotPlot_etc.pdf", sep=""), width = 9, height = 12)
# ggplot(data = dotplot_df_etc, mapping = aes_string(x = "variable", y = 'geneName_pig')) +
#   geom_point(mapping = aes_string(size = 'value', color = "Cell.type")) +
#   scale_radius(range = c(0, 10), breaks = c(2.5, 5.0, 7.5, 10.0)) +
#   labs(
#     x = "Marker Genes",
#     y = 'Cell type',
#     title = "Bulk"
#   ) +
#   geom_hline(yintercept=c(5.5, 10.5),lwd=0.5,colour="black") +
#   theme_bw()+
#   scale_color_manual(values=cbPalette) +
#   theme(axis.title.x = element_blank(), 
#         axis.title.y = element_blank(), 
#         axis.text.x = element_text(angle = 90, hjust=1, size = 15), 
#         axis.text.y = element_text(size = 20),
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size=20),
#         legend.key.size = unit(0.5, 'cm')) +
#   guides(colour = guide_legend(override.aes = list(size=10), order = 1), size = guide_legend(title = 'Average Normalized Counts'))
# dev.off()
# 
# ###### Macrophages ####################
# pdf(file = paste(vis_path, "DotPlot_macrophages.pdf", sep=""), width = 8.5, height = 5)
# ggplot(data = dotplot_df_macrophages, mapping = aes_string(x = "variable", y = 'geneName_pig')) +
#   geom_point(mapping = aes_string(size = 'value', color = "Cell.type")) +
#   # Check values in dotplot_df_macrophages and control ranges and breaks
#   scale_radius(range = c(10, 15), breaks = c(7.5, 10)) +
#   guides(size = guide_legend(title = 'Average Normalized Counts')) +
#   labs(
#     x = "Marker Genes",
#     y = 'Cell type',
#     title = "Bulk: Macrophage markers"
#   ) +
#   theme_bw()+
#   scale_color_manual(values=cbPalette) +
#   theme(axis.title.x = element_blank(), 
#         axis.title.y = element_blank(), 
#         axis.text.x = element_blank(), 
#         axis.text.y = element_text(size = 20),
#         axis.ticks.x = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size=20),
#         legend.key.size = unit(0.5, 'cm'))
# dev.off()
# 
# ######## CD163 & Siglec1 ###############
# dotplot_df <- assay(rld)[rownames(assay(rld)) %in% c("ENSSSCG00000033146", "ENSSSCG00000007146"),] %>% data.frame()
# 
# metadata <- read.csv(paste0(dir_path, "files_3_groups/DeSeq2_Full_Table_expressed_genes.csv"), sep = ",")
# colnames(metadata)[1] <- "NCBI_GeneID_pig"
# 
# ## calculate average expression of each group
# dotplot_df$Avg_norm_A <- rowMeans(dotplot_df[,c("A1", "A2", "A3", "A4")])
# dotplot_df$Avg_norm_B <- rowMeans(dotplot_df[,c("B1", "B2", "B3", "B4")])
# dotplot_df$Avg_norm_C <- rowMeans(dotplot_df[,c("C1", "C2", "C3", "C4")])
# dotplot_df <- dotplot_df[,c("Avg_norm_A", "Avg_norm_B", "Avg_norm_C")]
# 
# ## add metadata column
# dotplot_df <- merge(dotplot_df, metadata[,c("Ensembl_GeneID_pig", "geneName_human")], by.x = 0, by.y = "Ensembl_GeneID_pig", all.x = TRUE) %>% column_to_rownames(var = "Row.names")
# dotplot_df <- melt(dotplot_df)
# 
# 
# # ggplot
# 
# scale.func <- switch(
#   EXPR = "radius",
#   'radius' = scale_radius,
# )
# 
# cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# 
# pdf(file = paste(vis_path, "DotPlot_CD163_Sn.pdf", sep=""), width = 9, height = 3)
# ggplot(data = dotplot_df, mapping = aes_string(x = "variable", y = 'geneName_human')) +
#   geom_point(mapping = aes_string(size = 'value')) +
#   scale_radius(range = c(2, 10), breaks = c(2.5, 5.0, 7.5, 10.0)) +
#   labs(
#     x = "Marker Genes",
#     y = 'Cell type',
#     title = "LCM"
#   ) +
#   geom_hline(yintercept=c(5.5, 10.5),lwd=0.5,colour="black") +
#   theme_bw()+
#   scale_color_manual(values=cbPalette) +
#   theme(axis.title.x = element_blank(), 
#         axis.title.y = element_blank(), 
#         axis.text.x = element_text(angle = 90, hjust=1, size = 15), 
#         axis.text.y = element_text(size = 20),
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size=20),
#         legend.key.size = unit(0.5, 'cm')) +
#   guides(colour = guide_legend(override.aes = list(size=10), order = 1), size = guide_legend(title = 'Average Normalized Counts'))
# dev.off()
# 
# 
# ######## New macrophage marker dot plot#####################
# # Don't need to run the code above
# vis_path_LCM <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/Visualizations_3groups/"
# vis_path_MF <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/Visualizations_3groups/"
# 
# LCM_mac <- read.csv(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/Visualizations_3groups/Macrophagemarkers_for_dotplot.csv")
# MF_mac <- read.csv(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/Visualizations_3groups/Macrophagemarkers_for_dotplot.csv")
# 
# colnames(LCM_mac) <- c('geneName_human', "NaSn-", "NaSn+", "LuSn+")
# colnames(MF_mac) <- c('geneName_human', "NaSn-", "NaSn+", "LuSn+")
# 
# LCM_mac_melt <- LCM_mac %>% melt()
# MF_mac_melt <- MF_mac %>% melt()
# 
# LCM_mac_melt[,"value"] <- log2(LCM_mac_melt[,"value"])
# MF_mac_melt[,"value"] <- log2(MF_mac_melt[,"value"])
# 
# # ggplot
# 
# scale.func <- switch(
#   EXPR = "radius",
#   'radius' = scale_radius,
# )
# 
# cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# Dotplot_mf_new <- function(vis_path, mac_melt, breaks, title){
#   
#   pdf(file = paste(vis_path, "DotPlot_MacrophageMarkers_new.pdf", sep=""), width = 9, height = 3)
#   dotplot <- ggplot(data = mac_melt, mapping = aes_string(x = "variable", y = 'geneName_human')) +
#     geom_point(mapping = aes_string(size = 'value')) +
#     scale_radius(range = c(7, 12), breaks = breaks) +
#     labs(
#       x = "Marker Genes",
#       y = 'Cell type',
#       title = title
#     ) +
#     theme_bw()+
#     scale_color_manual(values=cbPalette) +
#     theme(axis.title.x = element_blank(), 
#           axis.title.y = element_blank(), 
#           axis.text.x = element_text(size = 20), 
#           axis.text.y = element_text(size = 20),
#           plot.title = element_text(hjust = 0.5),
#           legend.title = element_text(size = 20),
#           legend.text = element_text(size=20),
#           legend.key.size = unit(0.5, 'cm')) +
#     guides(colour = guide_legend(override.aes = list(size=10), order = 1), size = guide_legend(title = 'Average Normalized Counts'))
#   print(dotplot)
#   dev.off()
# }
# 
# Dotplot_mf_new(vis_path_LCM, LCM_mac_melt,c(9.0, 10.0), "LCM")
# Dotplot_mf_new(vis_path_MF, MF_mac_melt, c(15, 16), "FACS")
