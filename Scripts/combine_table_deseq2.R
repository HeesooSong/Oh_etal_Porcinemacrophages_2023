# Combine datatables for DeSeq2

##############################################################################
# Comment out the wanted directory path
# LCM
dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/files_3_groups/"
# ## MF
# dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/files_3_groups/"
##############################################################################

ba_fname <- "DeSeq2_B_NasalSnPos_vs_A_NasalSnNeg.tsv"
ca_fname <- "DeSeq2_C_LungSnPos_vs_A_NasalSnNeg.tsv"
cb_fname <- "DeSeq2_C_LungSnPos_vs_B_NasalSnPos.tsv"

# B/A
BvsA <- read.table(file = paste(dir_path, ba_fname, sep=""),
                   sep = "\t",
                   header = TRUE)
BvsA <- subset(BvsA, select = -c(stat, lfcSE))
rownames(BvsA) <- BvsA[,1]
BvsA <- BvsA[,-1]
colnames(BvsA) <- paste("BvsA", colnames(BvsA), sep = "_")

BvsA["ENSSSCG00000047610",]
# C/A
CvsA <- read.table(file = paste(dir_path, ca_fname, sep=""),
                   sep = "\t",
                   header = TRUE)
CvsA <- subset(CvsA, select = -c(stat, lfcSE, baseMean))
rownames(CvsA) <- CvsA[,1]
CvsA <- CvsA[,-1]
colnames(CvsA) <- paste("CvsA", colnames(CvsA), sep = "_")
CvsA["ENSSSCG00000047610",]

# C/B
CvsB <- read.table(file = paste(dir_path, cb_fname, sep=""),
                   sep = "\t",
                   header = TRUE)
CvsB <- subset(CvsB, select = -c(stat, lfcSE, baseMean))
colnames(CvsB) <- paste("CvsB", colnames(CvsB), sep = "_")
CvsB["ENSSSCG00000047610",]

BvsA[0:5,]
CvsA[0:5,]
CvsB[0:5,]

full_table <- merge(BvsA, CvsA, by=0)
full_table[0:5,]
rownames(full_table) <- full_table[,1]
full_table <- full_table[,-1]
full_table["ENSSSCG00000047610",]
full_table <- merge(full_table, CvsB, by=0)
full_table[0:50,]
dim(full_table)


### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
ENS_NCBI[0:5]
#Change NA to ""
ENS_NCBI[is.na(ENS_NCBI)] <- ""
colnames(full_table)
colnames(ENS_NCBI)
ENS_NCBI[ which(ENS_NCBI$Ensembl_GeneID_pig == "ENSSSCG00000047610"), ]


full_table_ENS_NCBI_added <- merge(ENS_NCBI, full_table, by.x = "Ensembl_GeneID_pig", by.y ="Row.names")
full_table_ENS_NCBI_added[0:5,] %>% colnames()
full_table_ENS_NCBI_added_reordered <- full_table_ENS_NCBI_added[, c("NCBI_GeneID_pig", "NCBI_GeneID_human", "Ensembl_GeneID_pig", "Ensembl_GeneID_human", "geneName_pig", "geneName_human", "BvsA_log2FoldChange", "CvsA_log2FoldChange", "CvsB_log2FoldChange", "BvsA_pvalue", "CvsA_pvalue", "CvsB_pvalue", "BvsA_padj", "CvsA_padj", "CvsB_padj", "BvsA_A1", "BvsA_A2", "BvsA_A3", "BvsA_A4", "BvsA_B1", "BvsA_B2", "BvsA_B3", "BvsA_B4", "BvsA_C1", "BvsA_C2", "BvsA_C3", "BvsA_C4")]
#full_table_ENS_NCBI_added_reordered <- full_table_ENS_NCBI_added[, c("NCBI_GeneID_pig", "NCBI_GeneID_human", "Ensembl_GeneID_pig", "Ensembl_GeneID_human", "geneName_pig", "geneName_human", "BvsA_log2FoldChange", "CvsA_log2FoldChange", "CvsB_log2FoldChange", "BvsA_pvalue", "CvsA_pvalue", "CvsB_pvalue", "BvsA_padj", "CvsA_padj", "CvsB_padj", "BvsA_A1", "BvsA_A2", "BvsA_A3", "BvsA_B1", "BvsA_B2", "BvsA_B3", "BvsA_C1", "BvsA_C2", "BvsA_C3")]
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_A1"] <- "A1"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_A2"] <- "A2"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_A3"] <- "A3"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_A4"] <- "A4"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_B1"] <- "B1"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_B2"] <- "B2"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_B3"] <- "B3"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_B4"] <- "B4"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_C1"] <- "C1"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_C2"] <- "C2"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_C3"] <- "C3"
colnames(full_table_ENS_NCBI_added_reordered)[colnames(full_table_ENS_NCBI_added_reordered) == "BvsA_C4"] <- "C4"
full_table_ENS_NCBI_added_reordered[0:100,]

dim(full_table_ENS_NCBI_added_reordered)
## Remove all rows that have zeros in all counts but not the ones that have extreme outliers!! 
full_table_ENS_NCBI_added_reordered_onlyexpressed <- full_table_ENS_NCBI_added_reordered[ rowSums(full_table_ENS_NCBI_added_reordered[,16:ncol(full_table_ENS_NCBI_added_reordered)]) > 0, ]
dim(full_table_ENS_NCBI_added_reordered_onlyexpressed)

write.table(full_table_ENS_NCBI_added_reordered_onlyexpressed, paste(dir_path, "DeSeq2_Full_Table_nozeroexpression.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(full_table_ENS_NCBI_added_reordered, paste(dir_path, "DeSeq2_Full_Table.tsv", sep=""), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)




