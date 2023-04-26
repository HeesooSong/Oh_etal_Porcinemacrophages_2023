library(ggplot2)
library(ggpubr)
library(scales)

#################################
###   Ver.2
#################################

# 0. Read in data

basedir <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/GO_BP_CC_MF_KEGG/"

GSEA_df <- read.csv(paste0(basedir, "GOtermsKEGGCTSA_list.csv"))

colnames(GSEA_df)[1] <- "method"


# 1. Plot function

combined_barplot <- function(method, category, title){
  
  plot_df <- GSEA_df[which(GSEA_df$method == method & GSEA_df$category == category),]
  plot_df$NES <- abs(plot_df$NES)
  plot_df$CellType <- factor(plot_df$CellType, levels = c("NaSn-", "NaSn+", "LuSn+"))
  plot_df <- plot_df[order(plot_df$CellType),]
  plot_df$content <- factor(plot_df$content, levels = rev(as.vector(factor(plot_df$content))))
  
  if(category == "Goterms"){
    
    # Keep legend
    barPlot_ES <- ggplot(plot_df) +
      geom_col(aes(x=content, y = NES, fill = CellType), width = 0.8) +
      labs(title = title) +
      scale_fill_discrete(labels = c(bquote("NaSn"^"-"), bquote("NaSn"^"+"), bquote("LuSn"^"+"))) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                         expand = expansion(mult = c(0, 0.1), add = c(0.05, 0))) +
      coord_flip() +
      theme_classic() +
      theme(plot.title = element_text(size=22, face = "bold"), 
            axis.text.y = element_text(size = 17),
            axis.text.x = element_text(size = 15),
            axis.title = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 17),
            legend.position = "top") +
    
    pdf(file = paste(basedir,method, "_", category, ".pdf", sep=""), width = 8, height = 5)
    
  } else{
    
    # Remove legend in other plots
    barPlot_ES <- ggplot(plot_df) +
      geom_col(aes(x=content, y = NES, fill = CellType), width = 0.8) +
      labs(title = title) +
      scale_fill_discrete(labels = c(bquote("NaSn"^"-"), bquote("NaSn"^"+"), bquote("LuSn"^"+"))) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                         expand = expansion(mult = c(0, 0.1), add = c(0.05, 0))) +
      coord_flip() +
      theme_classic() +
      theme(plot.title = element_text(size=22, face = "bold"),
            axis.text.y = element_text(size = 19),
            axis.text.x = element_text(size = 15),
            axis.title = element_blank(),
            legend.position = "none",
            legend.title = element_blank()) +
      
    pdf(file = paste(basedir,method, "_", category, ".pdf", sep=""), width = 8, height = 4)
  }
  print(barPlot_ES)
  dev.off()
  
  return(barPlot_ES)
}


# 2. Plot per analysis

## 1) LCM
plot_LCM_KEGG <- combined_barplot("LCM", "KEGG", "KEGG pathway")
plot_LCM_GO <- combined_barplot("LCM", "Goterms", "GO terms")
plot_LCM_CTSA <- combined_barplot("LCM", "CTSA", "CTSA")

## 2) FACS
plot_FACS_KEGG <- combined_barplot("FACS", "KEGG", "KEGG pathway")
plot_FACS_GO <- combined_barplot("FACS", "Goterms", "GO terms")
plot_FACS_CTSA <- combined_barplot("FACS", "CTSA", "CTSA")



# 3. Combine plots
plot_LCM_combined <- ggarrange(plot_LCM_GO, plot_LCM_KEGG, plot_LCM_CTSA, nrow = 3, align = "v", heights = c(0.8, 0.6, 0.6), common.legend = TRUE, legend = "top")
plot_FACS_combined <- ggarrange(plot_FACS_GO, plot_FACS_KEGG, plot_FACS_CTSA, nrow = 3, align = "v", heights = c(0.8, 0.6, 0.6), common.legend = TRUE, legend = "top")

pdf(file = paste(basedir, "KEGG_GO_CTSA_Combined.pdf", sep=""), width = 16, height = 16)
ggarrange(plot_FACS_combined, plot_LCM_combined, labels = c("A", "B"), font.label = list(size = 40), ncol = 2)
dev.off()
