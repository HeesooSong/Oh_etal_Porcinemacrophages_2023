
FACS <- read.csv(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/MF/files_3_groups/DeSeq2_Full_Table.csv")
LCM <- read.csv(file = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithPig4_switchAB/files_3_groups/DeSeq2_Full_Table.csv")

FACS$Avg_norm_A <- apply(FACS[,16:19], 1, mean)
FACS$Avg_norm_B <- apply(FACS[,20:23], 1, mean)
FACS$Avg_norm_C <- apply(FACS[,24:27], 1, mean)

LCM$Avg_norm_A <- apply(LCM[,16:19], 1, mean)
LCM$Avg_norm_B <- apply(LCM[,20:23], 1, mean)
LCM$Avg_norm_C <- apply(LCM[,24:27], 1, mean)

A_df <- data.frame(log2(cbind(FACS$Avg_norm_A, LCM$Avg_norm_A)))
B_df <- data.frame(log2(cbind(FACS$Avg_norm_B, LCM$Avg_norm_B)))
C_df <- data.frame(log2(cbind(FACS$Avg_norm_C, LCM$Avg_norm_C)))

A_df <- A_df[!is.infinite(A_df$X1) & !is.infinite(A_df$X2),]
B_df <- B_df[!is.infinite(B_df$X1) & !is.infinite(B_df$X2),]
C_df <- C_df[!is.infinite(C_df$X1) & !is.infinite(C_df$X2),]

Model_A <- lm(X2~X1,data=A_df)
Model_B <- lm(X2~X1,data=B_df)
Model_C <- lm(X2~X1,data=C_df)

format(summary(Model_A)$r.squared,digits=3)
format(summary(Model_B)$r.squared,digits=3)
format(summary(Model_C)$r.squared,digits=3)

#########################################################################
# Plot
#########################################################################

base = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Correlation_Plot/"

pdf(file = paste0(base, "CorrelationPlot_A.pdf"), width = 8, height = 8)
par(mar = c(7, 7, 7, 7), mgp = c(4, 1, 0))
plot(A_df$X1, A_df$X2, main = expression("NaSn"^"-"), xlab = "FACS", ylab = "LCM", cex.lab = 3, cex.axis = 2, cex.main = 3.5, col.main = "#F8766D")
abline(Model_A, col = "red", lwd = 3)
text(x = 3, y = 16,labels=expression(paste(R^2, "=", 0.586)), cex = 3, col = "black")
dev.off()

pdf(file = paste0(base, "CorrelationPlot_B.pdf"), width = 8, height = 8)
par(mar = c(7, 7, 7, 7), mgp = c(4, 1, 0))
plot(B_df$X1, B_df$X2, main = expression("NaSn"^"+"), xlab = "FACS", ylab = "LCM", cex.lab = 3, cex.axis = 2, cex.main = 3.5, col.main = "#00BA38")
abline(Model_B, col = "red", lwd = 3)
text(x = 3, y = 16,labels=expression(paste(R^2, "=", 0.584)), cex = 3, col = "black")
dev.off()

pdf(file = paste0(base, "CorrelationPlot_C.pdf"), width = 8, height = 8)
par(mar = c(7, 7, 7, 7), mgp = c(4, 1, 0))
plot(C_df$X1, C_df$X2, main = expression("LuSn"^"+"), xlab = "FACS", ylab = "LCM", cex.lab = 3, cex.axis = 2, cex.main = 3.5, col.main = "#619CFF")
abline(Model_C, col = "red", lwd = 3)
text(x = 3, y = 15,labels=expression(paste(R^2, "=", 0.56)), cex = 3, col = "black")
dev.off()


