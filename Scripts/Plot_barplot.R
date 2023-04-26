##### Make library size barplots

## Dependencies
# Package names
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c("limma","edgeR"))
library("limma")
library("edgeR")

# Package names
packages <- c("ggplot2", "RColorBrewer", "Rtsne", "dplyr", "stringr", "dplyr", "gplots")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

########################################### MF #################################
dir_path <- "/home/bygul/5_Projects/220914_barplots/"
f_name <- "Dayoung_MF_2019_htseq_counts_Trimmed.tsv"
countData <- read.table(file = paste(dir_path, f_name, sep=""),
                        sep = "\t",
                        header = TRUE, row.names = NULL)
rownames(countData) <- countData$X
countData <- countData[,-1]
DGEobject <- DGEList(counts = countData)
metaData <- read.table(file = "/home/bygul/5_Projects/220914_barplots/coldata.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = TRUE)
# Add info to DGEobject$samples
DGEobject$samples <- cbind(DGEobject$samples, metaData[,c("sample","condition")])
### CHECKING LIBRARY SIZES
DGEobject$samples$lib.size

# BARPLOT LIB SIZES (all samples)
condit.col <- DGEobject$samples$condition
levels(condit.col)<- list("#f7746c"="A_NasalSnNeg",
                          "#00bb38"="B_NasalSnPos",
                          "#609dff"="C_LungSnPos")
condit.cols <- as.character(condit.col)
dev.off()
opar = par(oma = c(8,6,4,8))
# par(mar = c(8,6,4,8)) # more margin: bottom, left, top, right
bp <- barplot(DGEobject$samples$lib.size*1e-6,
              legend = TRUE, 
              beside = TRUE,
              axisnames = FALSE,
              cex.lab=1.5,
              col = condit.cols,
              ylab = "Library size (millions)", ylim=c(0,10), border=condit.cols, cex.axis = 1.5)
#axis(1, labels = DGEobject$samples$sample, 
#     at = bp, cex.axis = 2, srt=45)
labs <- DGEobject$samples$sample
text(cex=1.5, x=bp+0.4, y=-0.3, labs, xpd=TRUE, srt=45, pos=2)
par(opar)
opar = par(oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend(x = "right", legend = levels(DGEobject$samples$condition), fill = levels(condit.col), bty = "n", y.intersp = 2)
dev.off()

########################################### LCM #################################
dir_path <- "/home/bygul/5_Projects/220914_barplots/"
f_name <- "Dayoung_LCM_2021_htseq_counts_Trimmed.tsv"

countData <- read.table(file = paste(dir_path, f_name, sep=""),
                        sep = "\t",
                        header = TRUE, row.names = NULL)
rownames(countData) <- countData$X
countData[1:5,]
dim(countData)
countData <- countData[,-1]
countData[1:5,]
DGEobject <- DGEList(counts = countData)
dim(DGEobject)
metaData <- read.table(file = "/home/bygul/5_Projects/220914_barplots/coldata.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = TRUE)
# Add info to DGEobject$samples
DGEobject$samples <- cbind(DGEobject$samples, metaData[,c("sample","condition")])

### CHECKING LIBRARY SIZES 
DGEobject$samples$lib.size

# BARPLOT LIB SIZES (all samples)
condit.col <- DGEobject$samples$condition
levels(condit.col)<- list("#f7746c"="A_NasalSnNeg",
                          "#00bb38"="B_NasalSnPos",
                          "#609dff"="C_LungSnPos")
condit.col <- as.character(condit.col)
##
par(mar = c(8,6,4,4)) # more margin: bottom, left, top, right
bp <- barplot(DGEobject$samples$lib.size*1e-6,
              axisnames = FALSE,
              main = NULL,
              cex.lab=1.5,
              col = condit.col,
              ylab = "Library size (millions)", ylim=c(0,10), border=condit.col, cex.axis = 1.5)
#axis(1, labels = DGEobject$samples$sample, 
#     at = bp, cex.axis = 2, srt=45)
labs <- DGEobject$samples$sample
text(cex=1.5, x=bp+0.5, y=-0.2, labs, xpd=TRUE, srt=45, pos=2)
dev.off()

##

par(mar = c(8,6,4,4)) # more margin: bottom, left, top, right
bp <- barplot(DGEobject$samples$lib.size*1e-6,
              axisnames = FALSE,
              main = NULL,
              cex.lab=1.5,
              col = condit.col,
              ylab = "Library size (millions)", ylim=c(0,3), border=condit.col, cex.axis = 1.5)
labs <- DGEobject$samples$sample

#axis(1, labels = DGEobject$samples$sample, 
#     at = bp, las = 2, cex.axis = 2)
text(cex=1.5, x=bp+0.4, y=-0.1, labs, xpd=TRUE, srt=45, pos=2)
dev.off()
