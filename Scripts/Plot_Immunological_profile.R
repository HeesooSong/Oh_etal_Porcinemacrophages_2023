
# 0. Read in data

basedir <- "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/Immunological_Profile/"

MF_ABC <- read.csv(paste0(basedir, "DeSeq2_MF_Full_Table_ABC_Updated.csv"))
MF_LvsN <- read.csv(paste0(basedir, "DeSeq2_MF_Full_Table_lungvsnose_HnameUpdated2.csv"))
LCM_ABC <- read.csv(paste0(basedir, "DeSeq2_LCM_Full_Table_ABC_Updated.csv"))
LCM_LvsN <- read.csv(paste0(basedir, "DeSeq2_LCM_Full_Table_lungvsnose_HnameUpdated2.csv"))

colnames(MF_ABC)[1] <- "NCBI_GeneID_pig"
colnames(MF_LvsN)[1] <- "NCBI_GeneID_pig"
colnames(LCM_ABC)[1] <- "NCBI_GeneID_pig"
colnames(LCM_LvsN)[1] <- "NCBI_GeneID_pig"

head(MF_ABC)


# 1. Plot function

circular_barplot <- function(data, datatype, analysis, celltype, title){
  
  ## 1) Prepare data frame for plot
  
  data[which(data$geneName_pig == ""), "geneName_pig"] <- data[which(data$geneName_pig == ""), "geneName_human"]
  
  plot_df = data[,c("geneName_pig", "group", paste0("Avg_norm_", celltype))]
  plot_df$group <- factor(plot_df$group)
  colnames(plot_df)[3] <- "value"
  plot_df$value <- log2(plot_df$value)
  plot_df[which(plot_df$value < 0),"value"] <- 0
  
  #### Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 3
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(plot_df$group), ncol(plot_df)) )
  colnames(to_add) <- colnames(plot_df)
  to_add$group <- rep(levels(plot_df$group), each=empty_bar)
  plot_df <- rbind(plot_df, to_add)
  plot_df <- plot_df %>% arrange(group)
  plot_df$id <- seq(1, nrow(plot_df))
  
  #### Get the name and the y position of each label
  label_data <- plot_df
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #### prepare a data frame for base lines
  base_data <- plot_df %>% 
    group_by(group) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(group_title=mean(c(start, end)))
  
  #### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  
  ## 2) Plot

  
  pdf(file = paste(basedir,"CircularPlot_", datatype, "_", analysis, "_", celltype, ".pdf", sep=""), width = 15, height = 15)
  
  # Note: the warning message "missing values" are caused by empty bars that we created above.
  CircularPlot <- ggplot(plot_df, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 12, xend = start, yend = 12), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 9, xend = start, yend = 9), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 6, xend = start, yend = 6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 3, xend = start, yend = 3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(plot_df$id),4), y = c(3, 6, 9, 12), label = c("3", "6", "9", "12") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
    
    geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
    ylim(-10,16.3) +
    theme_minimal() +
    theme(
      legend.position = "none", #"right"
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(c(-1,-1,-1,-1), "cm"),
      plot.title = element_blank()
    ) +
    coord_polar() + 
    geom_text(data=label_data, aes(x=id, y=value+0.1, label=geneName_pig, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
    annotate("text", label = title, x = 10, y = -10, size = 10) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = group_title, y = -1.5, label=group), hjust= 0.5, colour = "black", alpha=0.8, size=7, fontface="bold", inherit.aes = FALSE)
  
  print(CircularPlot)
  
  dev.off()
}

circular_barplot_legend <- function(data, datatype, analysis, celltype, title){
  
  ## 1) Prepare data frame for plot
  
  data[which(data$geneName_pig == ""), "geneName_pig"] <- data[which(data$geneName_pig == ""), "geneName_human"]
  
  plot_df = data[,c("geneName_pig", "group", "groupName", paste0("Avg_norm_", celltype))]
  plot_df$group <- factor(plot_df$group)
  plot_df$groupName <- factor(plot_df$groupName)
  colnames(plot_df)[4] <- "value"
  plot_df$value <- log2(plot_df$value)
  plot_df[which(plot_df$value < 0),"value"] <- 0
  
  #### Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 3
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(plot_df$group), ncol(plot_df)) )
  colnames(to_add) <- colnames(plot_df)
  to_add$group <- rep(levels(plot_df$group), each=empty_bar)
  plot_df <- rbind(plot_df, to_add)
  plot_df <- plot_df %>% arrange(group)
  plot_df$id <- seq(1, nrow(plot_df))
  
  #### Get the name and the y position of each label
  label_data <- plot_df
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #### prepare a data frame for base lines
  base_data <- plot_df %>% 
    group_by(group) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(group_title=mean(c(start, end)))
  
  #### prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  
  ## 2) Plot
  
  
  pdf(file = paste(basedir,"CircularPlot_", datatype, "_", analysis, "_", celltype, ".pdf", sep=""), width = 15, height = 15)
  
  # Note: the warning message "missing values" are caused by empty bars that we created above.
  CircularPlot <- ggplot(plot_df, aes(x=as.factor(id), y=value, fill=groupName)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(id), y=value, fill=groupName), stat="identity", alpha=0.5) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 12, xend = start, yend = 12), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 9, xend = start, yend = 9), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 6, xend = start, yend = 6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 3, xend = start, yend = 3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(plot_df$id),4), y = c(3, 6, 9, 12), label = c("3", "6", "9", "12") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
    
    geom_bar(aes(x=as.factor(id), y=value, fill=groupName), stat="identity", alpha=0.5) +
    ylim(-10,16.3) +
    theme_minimal() +
    theme(
      #legend.position = "right", #"none", #"right"
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(c(-1,-1,-1,-1), "cm"),
      plot.title = element_blank(),
      legend.position=c(0.85,0.8),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    coord_polar() + 
    geom_text(data=label_data, aes(x=id, y=value+0.1, label=geneName_pig, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
    annotate("text", label = title, x = 10, y = -10, size = 10) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = group_title, y = -1.5, label=group), hjust= 0.5, colour = "black", alpha=0.8, size=7, fontface="bold", inherit.aes = FALSE)
  
  print(CircularPlot)
  
  dev.off()
}

# 2. Plot for each analysis
## 1) MF
circular_barplot(MF_ABC, "MF", "AvsBvsC", "A", expression(atop("FACS-RNAseq", paste("NaSn"^"low"))))
circular_barplot(MF_ABC, "MF", "AvsBvsC", "B", expression(atop("FACS-RNAseq", paste("NaSn"^"mid"))))
circular_barplot(MF_ABC, "MF", "AvsBvsC", "C", expression(atop("FACS-RNAseq", paste("LuSn"^"high"))))

circular_barplot(MF_LvsN, "MF", "LungvsNose", "lung", "FACS-RNAseq \n Lung macrophages")
circular_barplot(MF_LvsN, "MF", "LungvsNose", "nose", "FACS-RNAseq \n Nasal macrophages")

## 2) LCM
circular_barplot_legend(LCM_ABC, "LCM", "AvsBvsC", "A", expression(atop("LCM-RNAseq", paste("NaSn"^"low"))))
circular_barplot(LCM_ABC, "LCM", "AvsBvsC", "B", expression(atop("LCM-RNAseq", paste("NaSn"^"mid"))))
circular_barplot(LCM_ABC, "LCM", "AvsBvsC", "C", expression(atop("LCM-RNAseq", paste("LuSn"^"high"))))

circular_barplot(LCM_LvsN, "LCM", "LungvsNose", "lung", "LCM-RNAseq \n Lung macrophages")
circular_barplot(LCM_LvsN, "LCM", "LungvsNose", "nose", "LCM-RNAseq \n Nasal macrophages")

# data = LCM_ABC
# datatype = "LCM"
# celltype = "A"
# title = "NaSn-"
