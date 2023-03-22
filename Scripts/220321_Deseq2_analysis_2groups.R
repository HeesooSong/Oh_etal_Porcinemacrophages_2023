################################### Lung vs nose
library(BiocManager)
# BiocManager::install("DESeq2")
library("DESeq2")
library("dplyr")
library("tibble")
library(ggplot2)
library(gplots)
library(ggrepel)
library(RColorBrewer)

########################################
### Using deseq2  
########################################
# ## LCM trimmed (uncomment this if you want to analyse the LCM data)
# dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/"
# vis_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/Visualizations_2groups/"
# f_name <- "Dayoung_LCM_2021_htseq_counts_Trimmed.tsv"
# ##

### MF trimmed (uncomment this if you want to analyse the MF data)
dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/"
vis_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/Visualizations_2groups/"
f_name <- "Dayoung_MF_2019_htseq_counts_Trimmed.tsv"
###

cts <- as.matrix(read.csv(paste(dir_path, f_name, sep=""),sep="\t",row.names="X"))
dim(cts)
# 
# #### IMPORTANT!! LCM: Change column (A4 <-> B4) (Uncomment for LCM)############################
# colnames(cts) <- c("A1", "A2", "A3", "B4", "B1", "B2", "B3", "A4", "C1", "C2", "C3", "C4")
# cts <- cts[,c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")]


#########################################################################
## A+B vs C (Nasal vs. Lung)
#########################################################################

# read the metadata 
coldata <- read.csv("C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/coldata_nosevslung.tsv", sep="\t", row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(ENS_NCBI)
ENS_NCBI[0:5,]
ENS_NCBI[is.na(ENS_NCBI)] <- ""
ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$geneName_pig <- ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$Ensembl_GeneID_pig
ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$geneName_human <- ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$Ensembl_GeneID_pig

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
filepath <- paste(vis_path, "PCA.pdf", sep="")
pdf(file = filepath, width = 8, height = 6)
par(mar = c(2,2,2,2))
plotPCA(vsd, intgroup=c("condition"))
dev.off()

## Calculate normalisationfactors
dds$sizeFactor
read_counts <- cts
read_counts[0:5,]
normalized_read_counts <- read_counts %*% diag(dds$sizeFactor)
normalized_read_counts[0:5,]
colnames(normalized_read_counts) <- colnames(read_counts)
normalized_read_counts[0:5,]


B_Lung_vs_A_Nasal <- results(dds, contrast=c("condition", "B_Lung", "A_Nasal"), independentFiltering = FALSE)
dim(B_Lung_vs_A_Nasal)
normalized_read_counts_added <- merge(B_Lung_vs_A_Nasal, normalized_read_counts, by=0)
write.table(normalized_read_counts_added, paste(dir_path, "files_2_groups/DeSeq2_B_Lung_vs_A_Nasal.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
summary(B_Lung_vs_A_Nasal)

B_Lung_vs_A_Nasal_filtered <- B_Lung_vs_A_Nasal[which(abs(B_Lung_vs_A_Nasal$log2FoldChange)>1 & B_Lung_vs_A_Nasal$padj < 0.05),]
dim(B_Lung_vs_A_Nasal_filtered)
normalized_read_counts_added_filtered <- merge(B_Lung_vs_A_Nasal_filtered, normalized_read_counts, by=0)
dim(normalized_read_counts_added_filtered)
write.table(normalized_read_counts_added_filtered, paste(dir_path, "files_2_groups/DeSeq2_B_Lung_vs_A_Nasal_filtered.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
B_Lung_vs_A_Nasal$expression <- ifelse(B_Lung_vs_A_Nasal$padj < 0.05 & abs(B_Lung_vs_A_Nasal$log2FoldChange) >= 1, 
                                                  ifelse(B_Lung_vs_A_Nasal$log2FoldChange> 1 ,'Up','Down'),
                                                  'Stable')

# Volcano Plot
test <- merge(ENS_NCBI, as.data.frame(B_Lung_vs_A_Nasal), by.x = "Ensembl_GeneID_pig", by.y =0)
test <- test[complete.cases(test), ]
test <- test %>%
  rownames_to_column(var="geneID")


options(ggrepel.max.overlaps = 20)
p <- ggplot(data = test, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=expression)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  geom_text_repel(
    data=test %>% filter(abs(log2FoldChange)>2 & padj<0.0000001), # Filter data first
    aes(label=geneName_human),
    size=2)

pdf(file = paste(vis_path, "volcanoplot_AvsB.pdf", sep=""), width = 8, height = 6)
p
dev.off()

# Heatmap
res = B_Lung_vs_A_Nasal[order(B_Lung_vs_A_Nasal$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = sigtab[(abs(sigtab$log2FoldChange)>1),]
sigtab = sigtab[order(-abs(sigtab$log2FoldChange)),]
# sigtab = head(sigtab, 20)
sigtab = sigtab[order(-sigtab$log2FoldChange),]
test <- merge(ENS_NCBI, as.data.frame(sigtab), by.x = "Ensembl_GeneID_pig", by.y =0)
rownames(test) <- test$Ensembl_GeneID_pig
summary(sigtab)
rownames(sigtab)
samples <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")
samples_pig <- c("NaSn-.Pig1", "NaSn-.Pig2", "NaSn-.Pig3", "NaSn-.Pig4",
                 "NaSn+.Pig1", "NaSn+.Pig2", "NaSn+.Pig3", "NaSn+.Pig4",
                 "LuSn+.Pig1", "LuSn+.Pig2", "LuSn+.Pig3", "LuSn+.Pig4")
# samples_pig <- c(expression("NaSn"^"low"*"_Pig1"), expression("NaSn"^"low"*"_Pig2"), expression("NaSn"^"low"*"_Pig3"), expression("NaSn"^"low"*"_Pig4"), 
#                  expression("NaSn"^"mid"*"_Pig1"), expression("NaSn"^"mid"*"_Pig2"), expression("NaSn"^"mid"*"_Pig3"), expression("NaSn"^"mid"*"_Pig4"),
#                  expression("LuSn"^"high"*"_Pig1"), expression("LuSn"^"high"*"_Pig2"), expression("LuSn"^"high"*"_Pig3"), expression("LuSn"^"high"*"_Pig4"))
#title <- expression("(NaSn"^"low"*" + NaSn"^"mid"*") vs LuSn"^"high")
title <- "(NaSn- + NaSn+) vs LuSn+"

##WARNING: takes more than 15 minutes to compute and plot
pdf(file = paste(vis_path, "heatmap_ABvsC.pdf", sep=""), width = 20, height = 15)
par(mar=c(4,1,1,1), cex.main = 3)
heatmap.2( assay(rld)[ rownames(sigtab),samples], scale="row", 
           trace="none", dendrogram="none", 
           margins = c(15, 5),
           Rowv = FALSE,
           labRow= FALSE,
           # labRow= test[rownames(sigtab),]$geneName_human,
           #labRow=as.expression(lapply(rownames(sigtab), function(a) bquote(italic(.(a))))),
           keysize = 0.5,
           key=TRUE,
           density.info=c("none"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(20),
           cexCol = 2,
           lhei=c(0.1, 1),
           labCol = samples_pig,
           main = title)
dev.off()


#########################################################################
## A vs B+C 
#########################################################################

# read the metadata 
coldata <- read.csv("C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/coldata_AvsBC.tsv", sep="\t", row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(ENS_NCBI)
ENS_NCBI[0:5,]
ENS_NCBI[is.na(ENS_NCBI)] <- ""
ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$geneName_pig <- ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$Ensembl_GeneID_pig
ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$geneName_human <- ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$Ensembl_GeneID_pig

## Deseq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)

dds
dds <- DESeq(dds)
resultsNames(dds)
rld <- rlog( dds)


## Calculate normalisationfactors
dds$sizeFactor
read_counts <- cts
read_counts[0:5,]
normalized_read_counts <- read_counts %*% diag(dds$sizeFactor)
normalized_read_counts[0:5,]
colnames(normalized_read_counts) <- colnames(read_counts)
normalized_read_counts[0:5,]


BC_vs_A <- results(dds, contrast=c("condition", "BC", "A"), independentFiltering = FALSE)
dim(BC_vs_A)
normalized_read_counts_added <- merge(BC_vs_A, normalized_read_counts, by=0)
write.table(normalized_read_counts_added, paste(dir_path, "files_2_groups/DeSeq2_BC_vs_A.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
summary(BC_vs_A)

BC_vs_A_filtered <- BC_vs_A[which(abs(BC_vs_A$log2FoldChange)>1 & BC_vs_A$padj < 0.05),]
dim(BC_vs_A_filtered)
normalized_read_counts_added_filtered <- merge(BC_vs_A_filtered, normalized_read_counts, by=0)
dim(normalized_read_counts_added_filtered)
write.table(normalized_read_counts_added_filtered, paste(dir_path, "files_2_groups/DeSeq2_BC_vs_A_filtered.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

BC_vs_A$expression <- ifelse(BC_vs_A$padj < 0.05 & abs(BC_vs_A$log2FoldChange) >= 1, 
                                       ifelse(BC_vs_A$log2FoldChange> 1 ,'Up','Down'),
                                      'Stable')

test <- merge(ENS_NCBI, as.data.frame(BC_vs_A), by.x = "Ensembl_GeneID_pig", by.y =0)
test <- test[complete.cases(test), ]
test <- test %>%
  rownames_to_column(var="geneID")


# Heatmap
res = BC_vs_A[order(BC_vs_A$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = sigtab[(abs(sigtab$log2FoldChange)>1),]
sigtab = sigtab[order(-abs(sigtab$log2FoldChange)),]
# sigtab = head(sigtab, 20)
sigtab = sigtab[order(-sigtab$log2FoldChange),]
test <- merge(ENS_NCBI, as.data.frame(sigtab), by.x = "Ensembl_GeneID_pig", by.y =0)
rownames(test) <- test$Ensembl_GeneID_pig
test = test[!(nchar(test$geneName_human) == 18),]
test = test[order(-test$log2FoldChange),]
summary(sigtab)
rownames(sigtab)
samples <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")
#title <- expression("NaSn"^"low"*" vs (NaSn"^"mid"*" + LuSn"^"high"*")")
title <- "NaSn- vs (NaSn+ vs LuSn+)"
pdf(file = paste(vis_path, "heatmap_AvsBC.pdf", sep=""), width = 20, height = 15)
par(mar=c(4,1,1,1), cex.main = 3)
heatmap.2( assay(rld)[ rownames(sigtab),samples], scale="row", 
           trace="none", dendrogram="none", 
           margins = c(15, 5),
           Rowv = FALSE,
           labRow= FALSE,
           # labRow= test[rownames(sigtab),]$geneName_human,
           #labRow=as.expression(lapply(rownames(sigtab), function(a) bquote(italic(.(a))))),
           keysize = 0.5,
           key=TRUE,
           density.info=c("none"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(20),
           cexCol = 2,
           lhei=c(0.1, 1),
           labCol = samples_pig,
           main = title)
dev.off()


#########################################################################
## A+C vs B 
#########################################################################

# read the metadata 
coldata <- read.csv("C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/coldata_ACvsB.tsv", sep="\t", row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(ENS_NCBI)
ENS_NCBI[0:5,]
ENS_NCBI[is.na(ENS_NCBI)] <- ""
ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$geneName_pig <- ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$Ensembl_GeneID_pig
ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$geneName_human <- ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$Ensembl_GeneID_pig

## Deseq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)

dds
dds <- DESeq(dds)
resultsNames(dds)
rld <- rlog( dds)

## Calculate normalisationfactors
dds$sizeFactor
read_counts <- cts
read_counts[0:5,]
normalized_read_counts <- read_counts %*% diag(dds$sizeFactor)
normalized_read_counts[0:5,]
colnames(normalized_read_counts) <- colnames(read_counts)
normalized_read_counts[0:5,]


AC_vs_B <- results(dds, contrast=c("condition", "AC", "B"), independentFiltering = FALSE)
dim(BC_vs_A)
normalized_read_counts_added <- merge(AC_vs_B, normalized_read_counts, by=0)
write.table(normalized_read_counts_added, paste(dir_path, "files_2_groups/DeSeq2_AC_vs_B.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
summary(AC_vs_B)

AC_vs_B_filtered <- AC_vs_B[which(abs(AC_vs_B$log2FoldChange)>1 & AC_vs_B$padj < 0.05),]
dim(AC_vs_B_filtered)
normalized_read_counts_added_filtered <- merge(AC_vs_B_filtered, normalized_read_counts, by=0)
dim(normalized_read_counts_added_filtered)
write.table(normalized_read_counts_added_filtered, paste(dir_path, "files_2_groups/DeSeq2_AC_vs_B_filtered.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

AC_vs_B$expression <- ifelse(AC_vs_B$padj < 0.05 & abs(AC_vs_B$log2FoldChange) >= 1, 
                             ifelse(AC_vs_B$log2FoldChange> 1 ,'Up','Down'),
                             'Stable')

test <- merge(ENS_NCBI, as.data.frame(AC_vs_B), by.x = "Ensembl_GeneID_pig", by.y =0)
test <- test[complete.cases(test), ]
test <- test %>%
  rownames_to_column(var="geneID")


# Heatmap
res = AC_vs_B[order(AC_vs_B$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = sigtab[(abs(sigtab$log2FoldChange)>1),]
sigtab = sigtab[order(-abs(sigtab$log2FoldChange)),]
# sigtab = head(sigtab, 20)
sigtab = sigtab[order(-sigtab$log2FoldChange),]
test <- merge(ENS_NCBI, as.data.frame(sigtab), by.x = "Ensembl_GeneID_pig", by.y =0)
rownames(test) <- test$Ensembl_GeneID_pig
test = test[!(nchar(test$geneName_human) == 18),]
test = test[order(-test$log2FoldChange),]
summary(sigtab)
rownames(sigtab)
samples <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")
#title <- expression("(NaSn"^"low"*"+LuSn"^"high"*") vs NaSn"^"mid")
title <- "(NaSn- + LuSn+) vs NaSn+"
pdf(file = paste(vis_path, "heatmap_ACvsB.pdf", sep=""), width = 20, height = 15)
par(mar=c(4,1,1,1), cex.main = 3)
heatmap.2( assay(rld)[ rownames(sigtab),samples], scale="row", 
           trace="none", dendrogram="none", 
           margins = c(15, 5),
           Rowv = FALSE,
           labRow= FALSE,
           # labRow= test[rownames(sigtab),]$geneName_human,
           #labRow=as.expression(lapply(rownames(sigtab), function(a) bquote(italic(.(a))))),
           keysize = 0.5,
           key=TRUE,
           density.info=c("none"),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(20),
           cexCol = 2,
           lhei=c(0.1, 1),
           labCol = samples_pig,
           main = title)
dev.off()


