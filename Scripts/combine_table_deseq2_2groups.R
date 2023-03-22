## for lung vs nasal
# Combine datatables for DeSeq2
##############################################################################
## Comment out the wanted directory path
## LCM
# dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/files_2_groups/"
## MF
dir_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/files_2_groups/"
##############################################################################

#############################################################################
####A+B vs C (Nasal vs. Lung)
#############################################################################

ba_fname = "DeSeq2_B_Lung_vs_A_Nasal.tsv"
# B/A
BvsA <- read.table(file = paste(dir_path, ba_fname, sep=""),
                   sep = "\t",
                   header = TRUE)
BvsA <- subset(BvsA, select = -c(stat, lfcSE))
rownames(BvsA) <- BvsA[,1]
BvsA <- BvsA[,-1]
colnames(BvsA) <- paste("BvsA", colnames(BvsA), sep = "_")

BvsA[0:5,]

### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(BvsA)
rownames(BvsA)
colnames(ENS_NCBI)
#Change NA to ""
ENS_NCBI[is.na(ENS_NCBI)] <- ""


full_table_ENS_NCBI_added <- merge(ENS_NCBI, BvsA, by.x = "Ensembl_GeneID_pig", by.y =0)
full_table_ENS_NCBI_added[0:50,]
full_table_ENS_NCBI_added_reordered <- full_table_ENS_NCBI_added[, c("NCBI_GeneID_pig", "NCBI_GeneID_human", "Ensembl_GeneID_pig", "Ensembl_GeneID_human", "geneName_pig", "geneName_human", "BvsA_log2FoldChange", "BvsA_pvalue", "BvsA_padj", "BvsA_A1", "BvsA_A2", "BvsA_A3", "BvsA_A4", "BvsA_B1", "BvsA_B2", "BvsA_B3", "BvsA_B4", "BvsA_C1", "BvsA_C2", "BvsA_C3", "BvsA_C4")]
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
full_table_ENS_NCBI_added_reordered[0:50,]

dim(full_table_ENS_NCBI_added_reordered)
full_table_ENS_NCBI_added_reordered_onlyexpressed <- full_table_ENS_NCBI_added_reordered[ rowSums(full_table_ENS_NCBI_added_reordered[,10:ncol(full_table_ENS_NCBI_added_reordered)]) > 0, ]
dim(full_table_ENS_NCBI_added_reordered_onlyexpressed)

write.table(full_table_ENS_NCBI_added_reordered_onlyexpressed, paste(dir_path, "DeSeq2_Full_Table_lungvsnose_nozeroexpression.csv", sep=""), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(full_table_ENS_NCBI_added_reordered, paste(dir_path, "DeSeq2_Full_Table_lungvsnose.csv", sep=""), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



#############################################################################
####A vs B+C 
#############################################################################

ba_fname = "DeSeq2_BC_vs_A.tsv"
# B/A
BvsA <- read.table(file = paste(dir_path, ba_fname, sep=""),
                   sep = "\t",
                   header = TRUE)
BvsA <- subset(BvsA, select = -c(stat, lfcSE))
rownames(BvsA) <- BvsA[,1]
BvsA <- BvsA[,-1]
colnames(BvsA) <- paste("BvsA", colnames(BvsA), sep = "_")

BvsA[0:5,]

### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(BvsA)
rownames(BvsA)
colnames(ENS_NCBI)
#Change NA to ""
ENS_NCBI[is.na(ENS_NCBI)] <- ""


full_table_ENS_NCBI_added <- merge(ENS_NCBI, BvsA, by.x = "Ensembl_GeneID_pig", by.y =0)
full_table_ENS_NCBI_added[0:50,]
full_table_ENS_NCBI_added_reordered <- full_table_ENS_NCBI_added[, c("NCBI_GeneID_pig", "NCBI_GeneID_human", "Ensembl_GeneID_pig", "Ensembl_GeneID_human", "geneName_pig", "geneName_human", "BvsA_log2FoldChange", "BvsA_pvalue", "BvsA_padj", "BvsA_A1", "BvsA_A2", "BvsA_A3", "BvsA_A4", "BvsA_B1", "BvsA_B2", "BvsA_B3", "BvsA_B4", "BvsA_C1", "BvsA_C2", "BvsA_C3", "BvsA_C4")]
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
full_table_ENS_NCBI_added_reordered[0:50,]

dim(full_table_ENS_NCBI_added_reordered)
full_table_ENS_NCBI_added_reordered_onlyexpressed <- full_table_ENS_NCBI_added_reordered[ rowSums(full_table_ENS_NCBI_added_reordered[,10:ncol(full_table_ENS_NCBI_added_reordered)]) > 0, ]
dim(full_table_ENS_NCBI_added_reordered_onlyexpressed)

write.table(full_table_ENS_NCBI_added_reordered_onlyexpressed, paste(dir_path, "DeSeq2_Full_Table_BCvsA_nozeroexpression.csv", sep=""), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(full_table_ENS_NCBI_added_reordered, paste(dir_path, "DeSeq2_Full_Table_BCvsA.csv", sep=""), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



#############################################################################
####A+C vs B 
#############################################################################

ba_fname = "DeSeq2_AC_vs_B.tsv"
# B/A
BvsA <- read.table(file = paste(dir_path, ba_fname, sep=""),
                   sep = "\t",
                   header = TRUE)
BvsA <- subset(BvsA, select = -c(stat, lfcSE))
rownames(BvsA) <- BvsA[,1]
BvsA <- BvsA[,-1]
colnames(BvsA) <- paste("BvsA", colnames(BvsA), sep = "_")

BvsA[0:5,]

### Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(BvsA)
rownames(BvsA)
colnames(ENS_NCBI)
#Change NA to ""
ENS_NCBI[is.na(ENS_NCBI)] <- ""


full_table_ENS_NCBI_added <- merge(ENS_NCBI, BvsA, by.x = "Ensembl_GeneID_pig", by.y =0)
full_table_ENS_NCBI_added[0:50,]
full_table_ENS_NCBI_added_reordered <- full_table_ENS_NCBI_added[, c("NCBI_GeneID_pig", "NCBI_GeneID_human", "Ensembl_GeneID_pig", "Ensembl_GeneID_human", "geneName_pig", "geneName_human", "BvsA_log2FoldChange", "BvsA_pvalue", "BvsA_padj", "BvsA_A1", "BvsA_A2", "BvsA_A3", "BvsA_A4", "BvsA_B1", "BvsA_B2", "BvsA_B3", "BvsA_B4", "BvsA_C1", "BvsA_C2", "BvsA_C3", "BvsA_C4")]
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
full_table_ENS_NCBI_added_reordered[0:50,]

dim(full_table_ENS_NCBI_added_reordered)
full_table_ENS_NCBI_added_reordered_onlyexpressed <- full_table_ENS_NCBI_added_reordered[ rowSums(full_table_ENS_NCBI_added_reordered[,10:ncol(full_table_ENS_NCBI_added_reordered)]) > 0, ]
dim(full_table_ENS_NCBI_added_reordered_onlyexpressed)

write.table(full_table_ENS_NCBI_added_reordered_onlyexpressed, paste(dir_path, "DeSeq2_Full_Table_ACvsB_nozeroexpression.csv", sep=""), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(full_table_ENS_NCBI_added_reordered, paste(dir_path, "DeSeq2_Full_Table_ACvsB.csv", sep=""), append = FALSE, sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
