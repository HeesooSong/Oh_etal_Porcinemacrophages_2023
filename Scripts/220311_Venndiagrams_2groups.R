library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venndiagram_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Venndiagrams_LCMWithoutPig4/"
###############################################################################
### LCM trimmed  (uncomment this if you want to analyse the LCM data)
dir_path_LCM <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithoutPig4/files_2_groups/"
###############################################################################

df_present_genes_nose_LCM <- read.csv(file = paste(dir_path_LCM, "present_nose.csv", sep=""), header = FALSE)
df_present_genes_lung_LCM <- read.csv(file = paste(dir_path_LCM, "present_lung.csv", sep=""), header = FALSE)

present_genes_nose_LCM <- df_present_genes_nose_LCM$V1
present_genes_lung_LCM <- df_present_genes_lung_LCM$V1

venn.diagram(
  x = list(present_genes_nose_LCM, present_genes_lung_LCM),
  category.names = c("LCM_nasal" , "LCM_lung"),
  filename = paste(venndiagram_path, "LCM_presenceabsence_2groups.png", sep=""),
  output=TRUE,
  
  lwd = 2,
  lty = 'blank',
  fill = myCol[0:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = myCol[0:2]
  
)

###############################################################################
### MF
dir_path_MF <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/files_2_groups/"
###############################################################################

df_present_genes_nose_MF <- read.csv(file = paste(dir_path_MF, "present_nose.csv", sep=""), header = FALSE)
df_present_genes_lung_MF <- read.csv(file = paste(dir_path_MF, "present_lung.csv", sep=""), header = FALSE)

present_genes_nose_MF <- df_present_genes_nose_MF$V1
present_genes_lung_MF <- df_present_genes_lung_MF$V1

venn.diagram(
  x = list(present_genes_nose_MF, present_genes_lung_MF),
  category.names = c("MF_nasal" , "MF_lung"),
  filename = paste(venndiagram_path, "MF_presenceabsence_2groups.png", sep=""),
  output=TRUE,
  
  lwd = 2,
  lty = 'blank',
  fill = myCol[0:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = myCol[0:2]
  
)

###############################################################################
### MF and LCM combined - presence/absence
###############################################################################
myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(present_genes_nose_LCM, present_genes_lung_LCM, present_genes_nose_MF, present_genes_lung_MF),
  category.names = c("LCM_nasal" , "LCM_lung", "MF_nasal", "MF_lung"),
  filename = paste(venndiagram_path, "LCM_MF_LungvsNose_presenceabsence_2groups.png", sep=""),
  output=TRUE,
  
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = myCol
  
)

###############################################################################
### MF and LCM combined - DEG 
###############################################################################
df_DEG_BvsA_LCM <- read.csv(file = paste(dir_path_LCM, "DEG_NoseVsLung.csv", sep=""), header = FALSE)
df_DEG_BvsA_MF <- read.csv(file = paste(dir_path_MF, "DEG_NoseVsLung.csv", sep=""), header = FALSE)

DEG_BvsA_LCM <- df_DEG_BvsA_LCM$V1
DEG_BvsA_MF <- df_DEG_BvsA_MF$V1

venn.diagram(
  x = list(DEG_BvsA_LCM, DEG_BvsA_MF),
  category.names = c("LCM_LungvsNose", "MF_LungvsNose"),
  filename = paste(venndiagram_path, "DEG_MF_LCM_2groups.png", sep=""),
  output=TRUE,
  
  lwd = 2,
  lty = 'blank',
  fill = myCol[0:2],
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = myCol[0:2]
)



