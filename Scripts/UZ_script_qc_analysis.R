
cts_files <- list.files()[grep(list.files(), pattern = "*trimmed*")]

cts <- read.table(file="RNA018271_trimmed_htseq_counts.txt", stringsAsFactors = FALSE)
cts <- cts[-c(25881:25885),]
rownames(cts) <- cts$V1

cts <- as.data.frame(cts)

for(i in 1:10){

cts_tmp <- read.table(file=cts_files[i], stringsAsFactors = FALSE)
cts_tmp <- cts_tmp[-c(25881:25885),]

cts[,i] <- cts_tmp$V2

}
colMaxs(as.matrix(cts))

colnames(cts) <- c("RNA018264","RNA018265","RNA018266","RNA018267","RNA018268",
                   "RNA018269","RNA018270","RNA018271","RNA018272","RNA018273")

colnames(cts) <- c("10_cells","50_cells","100_cells","300_cells","500_cells",
                   "10_cells_car","50_cells_car","100_cells_car","300_cells_car","500_cells_car")

colnames(cts) <- c("10_c","50_c","100_c","300_c","500_c",
                   "10_c_car","50_c_car","100_c_car","300_c_car","500_c_car")

write.table(cts,file="count_table.txt",quote=FALSE,sep="\t",col.names = TRUE, row.names = TRUE)
barplot(colSums(cts),las=3, ylab="total counts (raw)")

plot(log(cts[,5],10), log(cts[,10],10), xlab="500_cells", ylab="500_cells_carrier")
plot(log(cts[,4],10), log(cts[,9],10), xlab="300_cells", ylab="300_cells_carrier")
plot(log(cts[,3],10), log(cts[,8],10), xlab="100_cells", ylab="100_cells_carrier")
plot(log(cts[,2],10), log(cts[,7],10), xlab="50_cells", ylab="50_cells_carrier")
plot(log(cts[,1],10), log(cts[,6],10), xlab="10_cells", ylab="10_cells_carrier")
plot(log(cts[,9],10), log(cts[,10],10), xlab="300_cells_carrier", ylab="500_cells_carrier")
plot(log(cts[,4],10), log(cts[,5],10), xlab="300_cells", ylab="500_cells", xlim=c(0,5), ylim=c(0,5))
plot(log(cts[,9],10), log(cts[,10],10), xlab="300_cells_carrier", ylab="500_cells_carrier", xlim=c(0,5), ylim=c(0,5))


sum(cts$V2)
j <- 10
gene_cts <- c(
dim(cts[cts$`10_c`>j,])[1],
dim(cts[cts$`50_c`>j,])[1],
dim(cts[cts$`100_c`>j,])[1],
dim(cts[cts$`300_c`>j,])[1],
dim(cts[cts$`500_c`>j,])[1],
dim(cts[cts$`10_c_car`>j,])[1],
dim(cts[cts$`50_c_car`>j,])[1],
dim(cts[cts$`100_c_car`>j,])[1],
dim(cts[cts$`300_c_car`>j,])[1],
dim(cts[cts$`500_c_car`>j,])[1])

barplot(gene_cts,las=3, main = "#genes > 10 cts",names.arg = colnames(cts))

plot(density(log(cts$`300_c`,10)))
lines(density(log(cts$`300_c_car`,10)),col="red")
hist(log(cts$`500_c`,10))
hist(log(cts$`300_c`,10))

hist(log(cts$`500_c_car`,10),col="red",add=TRUE)

sampleFiles <- cts_files

sampleNames <- c("10_cells","50_cells","100_cells","300_cells","500_cells",
                 "10_cells_car","50_cells_car","100_cells_car","300_cells_car","500_cells_car")

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleNames)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = "./",
                                  design= ~condition)
keep <- rowSums(counts(dds)) > 0  #here we can go stricter but deseq2 will also filter for lowly expressed genes
dds <- dds[keep,]

# Run DE analysis
dds <- DESeq(dds)


# write normalized count table
data.norm<-counts(dds, normalize=T)


plot(log(data.norm[,5],10), log(data.norm[,10],10), xlab="500_cells", ylab="500_cells_carrier")
plot(log(data.norm[,4],10), log(data.norm[,9],10), xlab="300_cells", ylab="300_cells_carrier")
plot(log(data.norm[,3],10), log(data.norm[,8],10), xlab="100_cells", ylab="100_cells_carrier")
plot(log(data.norm[,2],10), log(data.norm[,7],10), xlab="50_cells", ylab="50_cells_carrier")
plot(log(data.norm[,1],10), log(data.norm[,6],10), xlab="10_cells", ylab="10_cells_carrier")

plot(log(data.norm[,9],10), log(data.norm[,10],10), xlab="300_cells_carrier", ylab="500_cells_carrier", 
     xlim=c(0,5), ylim=c(0,5))
plot(log(data.norm[,4],10), log(data.norm[,5],10), xlab="300_cells", ylab="500_cells", 
     xlim=c(0,5), ylim=c(0,5)) 
plot(log(data.norm[,9],10), log(data.norm[,10],10), xlab="300_cells_carrier", ylab="500_cells_carrier", 
     xlim=c(0,5), ylim=c(0,5),pch=19, cex=0.8)

barplot(colSums(data.norm), las=3)


