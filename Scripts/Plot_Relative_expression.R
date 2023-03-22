library(reshape2)
library(ggplot2)
library(dplyr)

# 0. Read in data

basedir <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Relative_expression/"

data <- read.csv(paste0(basedir, "Immunological_profile_A_BorB_A.csv"))


# 1. LCM
logFC_df <- data[which(data$RNAseq == "LCM"),c("Ensembl_GeneID_pig", "geneName_human", "Avg_norm_A", "Avg_norm_B", "Avg_norm_C")]
colnames(logFC_df)[3:5] <- c("NaSn-", "NaSn+", "LuSn+")
logFC_df <- melt(logFC_df)

log2AB <- log2(logFC_df$value[1]) - log2(1)
log2CB <- log2(logFC_df$value[3]) - log2(1)
logFC_df$log2FC <- c(log2AB, log2(1), log2CB)
logFC_df_LCM <- logFC_df


# 2. FACS
logFC_df <- data[which(data$RNAseq == "MF"),c("Ensembl_GeneID_pig", "geneName_human", "Avg_norm_A", "Avg_norm_B", "Avg_norm_C")]
logFC_df$minimum <- apply(logFC_df[,3:5], 1, min) # calculate minimum
logFC_df$minimum <- 0.5 * logFC_df$minimum
logFC_df <- logFC_df %>% mutate_at(vars(3:5), list(relative=~./minimum))  #calculate relative expression
logFC_df <- logFC_df[,c(1,2,7,8,9)]
colnames(logFC_df)[3:5] <- c("NaSn-", "NaSn+", "LuSn+")

logFC_df <- melt(logFC_df)
logFC_df$geneName_human <- factor(logFC_df$geneName_human, levels = c("CCR2", "IL1B", "CXCL6", "CD1A", "CD4"))
logFC_df$log2FC <- log2(logFC_df$value)


# 3. Combine LCM & FACS and Plot

logFC_df <- rbind(logFC_df, logFC_df_LCM)

pdf(file = paste(basedir,"Relative_expression_LCM_FACS.pdf", sep=""), width = 12, height = 7)
relative_expression_bulk <- ggplot(logFC_df, aes(x=variable, y=log2FC, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("NaSn-"="#F8766D", "NaSn+"="#00BA38", "LuSn+"="#619CFF"),
                    labels = c(bquote("NaSn"^"low"), bquote("NaSn"^"mid"), bquote("LuSn"^"high")))+
  #geom_text(data = logFC_df, aes(label = variable, y = (log2 + 0.4)), size = 5) +
  theme_classic() +
  labs(y = "Relative Expression") +
  facet_wrap(~ geneName_human, ncol = 6) +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.key.size = unit(1, 'cm'))
print(relative_expression_bulk)
dev.off()