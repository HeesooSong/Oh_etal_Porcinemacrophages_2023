#########################################################################
# Dot Plot
#########################################################################

library(tibble)
library(RColorBrewer)
library(ggpubr)
library(reshape2)


func_dotplot <- function(df, breaks, filename, labels){
  ## melt expression data for ggplot
  dotplot_df <- melt(df)
  dotplot_df$rowname <- factor(dotplot_df$rowname, levels = rev(unique(dotplot_df$rowname)))
  dotplot_df$variable <- factor(dotplot_df$variable)
  
  plot <- ggplot(data = dotplot_df, mapping = aes_string(x = "variable", y = 'rowname')) +
    geom_point(mapping = aes_string(size = 'value', fill = "value"), color = "black", pch = 21) +
    scale_radius(range = c(1, 12), breaks = breaks, trans = "sqrt") +
    scale_fill_gradientn(
      colors = brewer.pal(9, "YlGnBu"),
      breaks = breaks,
      space = "Lab",
      guide = guide_colorbar(title = "Average\nNormalized Counts", order = 1),
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw()+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 40, hjust=1, size = 12, color = "black"), 
          axis.text.y = element_text(size = 15, color = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size=12),
          legend.key.size = unit(0.5, 'cm'),
          panel.border = element_blank()) +
    guides(size = guide_legend(title = ' ', shape = 21))
  
  pdf(file = paste(vis_path, filename, sep=""), width = 4, height = 5)
  plot
  dev.off()
  
  return(plot)
}


# Read in marker gene data
vis_path <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/dotplot_selectedgenes/"
marker_genes_df <- read.csv("C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/dotplot_selectedgenes/RNAexpression_of_protein_leve_confirmed_gene.csv", sep = ",", check.names = FALSE)
colnames(marker_genes_df)[1] <- "FACS_NaSn-"

FACS_df <- marker_genes_df[,1:3] %>% rownames_to_column()
LCM_df <- marker_genes_df[,4:6] %>% rownames_to_column()

dotplot_FACS <- func_dotplot(FACS_df, c(0, 5000, 10000, 15000, 20000), "DotPlot_FACS.pdf", c(expression("FACS NaSn"^"-"), expression("FACS NaSn"^"+"), expression("FACS LuSn"^"+")))
dotplot_LCM <- func_dotplot(LCM_df, c(0, 1000, 2000), "DotPlot_LCM.pdf", c(expression("LCM NaSn"^"-"), expression("LCM NaSn"^"+"), expression("LCM LuSn"^"+")))

pdf(file = paste(vis_path, "DotPlot_combined.pdf", sep=""), width = 8, height = 5)
ggarrange(dotplot_FACS, dotplot_LCM, ncol = 2)
dev.off()
