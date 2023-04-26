#install.packages("extrafont")
library(VennDiagram)
library(RColorBrewer)
library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
library(extrafont)
font_import()
loadfonts(device = "win")
#fonts()

myCol <- c("#FC8D62", "#66C2A5", "#8DA0CB")
myCol_DEG <- c("#D39200", "#DB72FB", "#00798C")

venndiagram_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Venndiagrams_LCMswitchAB/"

## Read in the Human and Pig table
ENS_NCBI <- read.table(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Ensembl_NCBI_files/New_full_correctorder.tsv",
                       sep = "\t",
                       header = TRUE, row.names = NULL)
colnames(ENS_NCBI)
ENS_NCBI[0:5,]
ENS_NCBI[is.na(ENS_NCBI)] <- ""
ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$geneName_pig <- ENS_NCBI[which(ENS_NCBI$geneName_pig==""),]$Ensembl_GeneID_pig
ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$geneName_human <- ENS_NCBI[which(ENS_NCBI$geneName_human==""),]$Ensembl_GeneID_pig



###############################################################################
### LCM Presence vs absence
dir_path_LCM <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/files_3_groups/"

# Read in the manually cleaned out full table with verified gene names
clean_table_LCM <- as.data.frame(read.csv(paste(dir_path_LCM, "DeSeq2_Full_Table_nozeroexpression_final.csv", sep=""),sep=","))
colnames(clean_table_LCM)[1] <- "Ensembl_GeneID_pig"
###############################################################################
df_present_genes_nose1_LCM <- read.csv(file = paste(dir_path_LCM, "present_in_nose1.csv", sep=""), header = FALSE)
df_present_genes_nose2_LCM <- read.csv(file = paste(dir_path_LCM, "present_in_nose2.csv", sep=""), header = FALSE)
df_present_genes_lung_LCM <- read.csv(file = paste(dir_path_LCM, "present_in_lung.csv", sep=""), header = FALSE)

present_genes_nose1_LCM <- df_present_genes_nose1_LCM$V1
present_genes_nose2_LCM <- df_present_genes_nose2_LCM$V1
present_genes_lung_LCM <- df_present_genes_lung_LCM$V1

venn.diagram(
  x = list(present_genes_nose1_LCM, present_genes_nose2_LCM, present_genes_lung_LCM),
  category.names = c(expression("NaSn"^"-"), expression("NaSn"^"+"), expression("LuSn"^"+")),
  filename = paste(venndiagram_path, "LCM_presenceabsence_3groups.tif", sep=""),
  imagetype = "tiff",
  output=TRUE,
  print.mode = c("raw", "percent"),
  margin = 0.07,
  
  main.fontfamily = "Arial",
  main.cex = 1.5,
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1.7,
  fontface = "bold",
  fontfamily = "Arial",
  cat.col = myCol,
  cat.fontfamily = "Arial",
  cat.cex = 2.5,
  cat.just=list(c(0.3,0) , c(0.7,0) , c(0.5,0.2))
)

###############################################################################
### MF Presence vs absence
dir_path_MF <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/files_3_groups/"

# Read in the manually cleaned out full table with verified gene names
clean_table_MF <- as.data.frame(read.csv(paste(dir_path_MF, "DeSeq2_Full_Table_nozeroexpression_final.csv", sep=""),sep=","))
colnames(clean_table_MF)[1] <- "Ensembl_GeneID_pig"
###############################################################################
## Presence vs absence
df_present_genes_nose1_MF <- read.csv(file = paste(dir_path_MF, "present_in_nose1.csv", sep=""), header = FALSE)
df_present_genes_nose2_MF <- read.csv(file = paste(dir_path_MF, "present_in_nose2.csv", sep=""), header = FALSE)
df_present_genes_lung_MF <- read.csv(file = paste(dir_path_MF, "present_in_lung.csv", sep=""), header = FALSE)

present_genes_nose1_MF <- df_present_genes_nose1_MF$V1
present_genes_nose2_MF <- df_present_genes_nose2_MF$V1
present_genes_lung_MF <- df_present_genes_lung_MF$V1

venn.diagram(
  x = list(present_genes_nose1_MF, present_genes_nose2_MF, present_genes_lung_MF),
  category.names = c(expression("NaSn"^"-"), expression("NaSn"^"+"), expression("LuSn"^"+")),
  filename = paste(venndiagram_path, "MF_presenceabsence_3groups.tif", sep=""),
  imagetype = "tiff",
  output=TRUE,
  print.mode = c("raw", "percent"),
  margin = 0.07,
  
  main.fontfamily = "Arial",
  main.cex = 1.5,
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1.7,
  fontface = "bold",
  fontfamily = "Arial",
  cat.col = myCol,
  cat.fontfamily = "Arial",
  cat.cex = 2.5,
  cat.just=list(c(0.3,0) , c(0.7,0) , c(0.5,0.2))
)


###############################################################################
### LCM DEG
###############################################################################
## Venn diagram LCM 3 groups 
df_DEG_BvsA_LCM <- read.csv(file = paste(dir_path_LCM, "DEG_BvsA.csv", sep=""), header = FALSE)
df_DEG_CvsA_LCM <- read.csv(file = paste(dir_path_LCM, "DEG_CvsA.csv", sep=""), header = FALSE)
df_DEG_CvsB_LCM <- read.csv(file = paste(dir_path_LCM, "DEG_CvsB.csv", sep=""), header = FALSE)

#add gene symbol
# df_DEG_BvsA_LCM <- merge(df_DEG_BvsA_LCM, ENS_NCBI[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
# df_DEG_CvsA_LCM <- merge(df_DEG_CvsA_LCM, ENS_NCBI[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
# df_DEG_CvsB_LCM <- merge(df_DEG_CvsB_LCM, ENS_NCBI[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
df_DEG_BvsA_LCM <- merge(df_DEG_BvsA_LCM, clean_table_LCM[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
df_DEG_CvsA_LCM <- merge(df_DEG_CvsA_LCM, clean_table_LCM[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
df_DEG_CvsB_LCM <- merge(df_DEG_CvsB_LCM, clean_table_LCM[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)

# Remove NA in geneName
df_DEG_BvsA_LCM <- na.omit(df_DEG_BvsA_LCM)
df_DEG_CvsA_LCM <- na.omit(df_DEG_CvsA_LCM)
df_DEG_CvsB_LCM <- na.omit(df_DEG_CvsB_LCM)

# Remove "" in geneName
df_DEG_BvsA_LCM <- df_DEG_BvsA_LCM[which(df_DEG_BvsA_LCM$geneName_pig != ""),]
df_DEG_CvsA_LCM <- df_DEG_CvsA_LCM[which(df_DEG_CvsA_LCM$geneName_pig != ""),]
df_DEG_CvsB_LCM <- df_DEG_CvsB_LCM[which(df_DEG_CvsB_LCM$geneName_pig != ""),]


DEG_BvsA_LCM <- df_DEG_BvsA_LCM$geneName_pig
DEG_CvsA_LCM <- df_DEG_CvsA_LCM$geneName_pig
DEG_CvsB_LCM <- df_DEG_CvsB_LCM$geneName_pig

venn_DEG <- venn.diagram(
  x = list(DEG_BvsA_LCM, DEG_CvsA_LCM, DEG_CvsB_LCM),
  category.names = c(expression("NaSn"^"-"*" vs NaSn"^"+"), expression("LuSn"^"+"*" vs NaSn"^"-"), expression("LuSn"^"+"*" vs NaSn"^"+")),
  filename = paste(venndiagram_path, "LCM_DEG_3groups.tif", sep=""),
  imagetype = "tiff",
  output=TRUE,
  margin = 0.07,

  main = "LCM: DEG",
  main.fontfamily = "Arial",
  main.cex = 1.5,
  lwd = 2,
  lty = 'blank',
  fill = myCol_DEG,
  
  # Numbers
  cex = 2.5,
  fontface = "bold",
  fontfamily = "Arial",
  cat.col = myCol_DEG,
  cat.fontfamily = "Arial",
  cat.cex = 2,
  cat.just=list(c(0.3,-0.3) , c(0.7,-0.3) , c(0.5,0.25))
  
)

# legend <- legendGrob(labels=c("A","B"), pch=rep(19,length(c("A","B"))),
#                      gp=gpar(col=myCol_DEG, fill="gray"),
#                      byrow=TRUE)
# g <- gTree(children = gList(venn_DEG))
# gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))


# Export intersection information
LCM_DEG <- calculate.overlap(
  x = list(
    "NaSn-vsNaSn+" = DEG_BvsA_LCM,
    "LuSn+vsNaSn-" = DEG_CvsA_LCM,
    "LuSn+vsNaSn+" = DEG_CvsB_LCM
  )
)

LCM_DEG_matrix <- sapply(LCM_DEG, '[', seq(max(sapply(LCM_DEG, length))))
LCM_DEG_matrix[is.na(LCM_DEG_matrix)] <- ""
write.csv(LCM_DEG_matrix, file = paste0(venndiagram_path, "LCM_DEG.csv"), row.names = FALSE)

###############################################################################
### MF DEG
###############################################################################
df_DEG_BvsA_MF <- read.csv(file = paste(dir_path_MF, "DEG_BvsA.csv", sep=""), header = FALSE)
df_DEG_CvsA_MF <- read.csv(file = paste(dir_path_MF, "DEG_CvsA.csv", sep=""), header = FALSE)
df_DEG_CvsB_MF <- read.csv(file = paste(dir_path_MF, "DEG_CvsB.csv", sep=""), header = FALSE)

#add gene symbol
# df_DEG_BvsA_MF <- merge(df_DEG_BvsA_MF, ENS_NCBI[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
# df_DEG_CvsA_MF <- merge(df_DEG_CvsA_MF, ENS_NCBI[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
# df_DEG_CvsB_MF <- merge(df_DEG_CvsB_MF, ENS_NCBI[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)

#### Filter method 1: Exactly the same filtering method with DE analysis
# clean_geneName_BvsA_MF <- clean_table_MF[which(clean_table_MF$Ensembl_GeneID_pig %in% df_DEG_BvsA_MF$V1),c("Ensembl_GeneID_pig", "geneName_pig")]
# clean_geneName_CvsA_MF <- clean_table_MF[which(clean_table_MF$Ensembl_GeneID_pig %in% df_DEG_CvsA_MF$V1),c("Ensembl_GeneID_pig", "geneName_pig")]
# clean_geneName_CvsB_MF <- clean_table_MF[which(clean_table_MF$Ensembl_GeneID_pig %in% df_DEG_CvsB_MF$V1),c("Ensembl_GeneID_pig", "geneName_pig")]
# 
# wrong_geneName_BvsA_MF <- clean_geneName_BvsA_MF[which(clean_geneName_BvsA_MF$geneName_pig == ""),"Ensembl_GeneID_pig"]
# wrong_geneName_CvsA_MF <- clean_geneName_CvsA_MF[which(clean_geneName_CvsA_MF$geneName_pig == ""),"Ensembl_GeneID_pig"]
# wrong_geneName_CvsB_MF <- clean_geneName_CvsB_MF[which(clean_geneName_CvsB_MF$geneName_pig == ""),"Ensembl_GeneID_pig"]
# 
# wrong_geneName_BvsA_MF <- c(wrong_geneName_BvsA_MF, setdiff(df_DEG_BvsA_MF$V1, clean_geneName_BvsA_MF$Ensembl_GeneID_pig))
# wrong_geneName_CvsA_MF <- c(wrong_geneName_CvsA_MF, setdiff(df_DEG_CvsA_MF$V1, clean_geneName_CvsA_MF$Ensembl_GeneID_pig))
# wrong_geneName_CvsB_MF <- c(wrong_geneName_CvsB_MF, setdiff(df_DEG_CvsB_MF$V1, clean_geneName_CvsB_MF$Ensembl_GeneID_pig))
# 
# df_DEG_BvsA_MF <- df_DEG_BvsA_MF[!(df_DEG_BvsA_MF$V1 %in% wrong_geneName_BvsA_MF),]
# df_DEG_CvsA_MF <- df_DEG_CvsA_MF[!(df_DEG_CvsA_MF$V1 %in% wrong_geneName_CvsA_MF),]
# df_DEG_CvsB_MF <- df_DEG_CvsB_MF[!(df_DEG_CvsB_MF$V1 %in% wrong_geneName_CvsB_MF),]

#### Filter method 2: Since we are only dealing with DE genes here, use curated table from the first hand. Same result gained.
df_DEG_BvsA_MF <- merge(df_DEG_BvsA_MF, clean_table_MF[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
df_DEG_CvsA_MF <- merge(df_DEG_CvsA_MF, clean_table_MF[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)
df_DEG_CvsB_MF <- merge(df_DEG_CvsB_MF, clean_table_MF[,c("Ensembl_GeneID_pig", "geneName_pig")], by.x = "V1", by.y = "Ensembl_GeneID_pig", all.x = TRUE)

# Remove NA in geneName
df_DEG_BvsA_MF <- na.omit(df_DEG_BvsA_MF)
df_DEG_CvsA_MF <- na.omit(df_DEG_CvsA_MF)
df_DEG_CvsB_MF <- na.omit(df_DEG_CvsB_MF)

# Remove "" in geneName
df_DEG_BvsA_MF <- df_DEG_BvsA_MF[which(df_DEG_BvsA_MF$geneName_pig != ""),]
df_DEG_CvsA_MF <- df_DEG_CvsA_MF[which(df_DEG_CvsA_MF$geneName_pig != ""),]
df_DEG_CvsB_MF <- df_DEG_CvsB_MF[which(df_DEG_CvsB_MF$geneName_pig != ""),]



DEG_BvsA_MF <- df_DEG_BvsA_MF$geneName_pig
DEG_CvsA_MF <- df_DEG_CvsA_MF$geneName_pig
DEG_CvsB_MF <- df_DEG_CvsB_MF$geneName_pig

venn_DEG <- venn.diagram(
  x = list(DEG_BvsA_MF, DEG_CvsA_MF, DEG_CvsB_MF),
  category.names = c(expression("NaSn"^"-"*" vs NaSn"^"+"), expression("LuSn"^"+"*" vs NaSn"^"-"), expression("LuSn"^"+"*" vs NaSn"^"+")),
  filename = paste(venndiagram_path, "MF_DEG_3groups.tif", sep=""),
  imagetype = "tiff",
  output=TRUE,
  margin = 0.07,

  main = "Bulk: DEG",
  main.fontfamily = "Arial",
  main.cex = 1.5,
  lwd = 2,
  lty = 'blank',
  fill = myCol_DEG,
  
  # Numbers
  cex = 2.5,
  fontface = "bold",
  fontfamily = "Arial",
  cat.col = myCol_DEG,
  cat.fontfamily = "Arial",
  cat.cex = 2,
  cat.just=list(c(0.3,-0.3) , c(0.7,-0.3) , c(0.5,0.25))
  
)

# legend <- legendGrob(labels=c("NaSn-vsNaSn+","LuSn+vsNaSn-", "LuSn+vsNaSn+"), pch=rep(19,length(c("NaSn-vsNaSn+","LuSn+vsNaSn-", "LuSn+vsNaSn+"))),
#                      gp=gpar(col=myCol_DEG, fill="gray"),
#                      byrow=TRUE)
# g <- gTree(children = gList(venn_DEG))
# gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))


# Export intersection information
MF_DEG <- calculate.overlap(
  x = list(
    "NaSn-vsNaSn+" = DEG_BvsA_MF,
    "LuSn+vsNaSn-" = DEG_CvsA_MF,
    "LuSn+vsNaSn+" = DEG_CvsB_MF
  )
)

MF_DEG_matrix <- sapply(MF_DEG, '[', seq(max(sapply(MF_DEG, length))))
MF_DEG_matrix[is.na(MF_DEG_matrix)] <- ""
write.csv(MF_DEG_matrix, file = paste0(venndiagram_path, "MF_DEG.csv"), row.names = FALSE)
