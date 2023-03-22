library(ggplot2)
library(patchwork)


# 0. Read in data

basedir <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/GO_BP_CC_MF_KEGG/"

MF_KEGG <- read.csv(paste0(basedir, "GOterms_KEGG_MF.csv"))
LCM_KEGG <- read.csv(paste0(basedir, "GOterms_KEGG_LCM.csv"))

colnames(MF_KEGG)[1] <- "GS"
colnames(LCM_KEGG)[1] <- "GS"

head(MF_KEGG)


# 1. Plot function

combined_barplot <- function(data, datatype, title, analysis){
  
  plot_df = data[which(data$analysis == analysis & !(duplicated(data$GS))),]
  plot_df$normalized_enrichment_score <- abs(plot_df$normalized_enrichment_score)
  plot_df$celltype <- factor(plot_df$celltype, levels = c("NaSn-", "NaSn+", "LuSn+"))
  plot_df <- plot_df[order(plot_df$celltype),]
  plot_df$GS <- factor(plot_df$GS, levels = rev(as.vector(factor(plot_df$GS))))
  
  pdf(file = paste(basedir,datatype, "_", analysis, ".pdf", sep=""), width = 15, height = 15)
  
  barPlot_ES <- ggplot(plot_df) +
    geom_col(aes(x=GS, y = normalized_enrichment_score, fill = celltype)) +
    geom_hline(yintercept = 0,colour = "grey90") +
    labs(y = "Normalized Enrichment Score", title = title) +
    scale_fill_discrete(labels = c(bquote("NaSn"^"low"), bquote("NaSn"^"mid"), bquote("LuSn"^"high"))) +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    theme(plot.title = element_text(size=22, face = "bold"), legend.position="top", legend.justification = "left",
          plot.margin = unit(c(0.5,0,0.5,0.5), "cm"), panel.grid.major = element_line(color = "grey90"),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 15))
  
  barPlot_size <- ggplot(plot_df) +
    geom_col(aes(x=GS, y = size), fill = "grey") +
    geom_text(aes(label = size, x = GS, y = (size + 20)), size =7) +
    labs(y = "Size") +
    coord_flip() +
    theme_classic() +
    theme(plot.margin = unit(c(2.8,0.5,0.5,0), "cm"), panel.grid.major = element_line(color = "grey90"), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          axis.title.x = element_text(size = 20)) +
    guides(y = "none")
    
  print((barPlot_ES + barPlot_size))

  dev.off()
}


# 2. Plot per analysis

## 1) LCM
combined_barplot(LCM_KEGG, "LCM", "KEGG pathway", "KEGG")
combined_barplot(LCM_KEGG, "LCM", "Biological Process", "BP")
combined_barplot(LCM_KEGG, "LCM", "Cellular Component", "CC")
combined_barplot(LCM_KEGG, "LCM", "Molecular Function", "MF")
combined_barplot(LCM_KEGG, "LCM", "Cell Type Signature Analysis", "CTSA")

## 2) MF
combined_barplot(MF_KEGG, "MF", "KEGG pathway", "KEGG")
combined_barplot(MF_KEGG, "MF", "Biological Process", "BP")
combined_barplot(MF_KEGG, "MF", "Cellular Component", "CC")
combined_barplot(MF_KEGG, "MF", "Molecular Function", "MF")
combined_barplot(MF_KEGG, "MF", "Cell Type Signature Analysis", "CTSA")
