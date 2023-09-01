
library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
library(tidyverse)
library(factoextra)
library(grid)
library(ggsignif)
library(stringr)
library(reshape2)
library(gtools)
library(tidyr)
library(ggbeeswarm) #https://github.com/eclarke/ggbeeswarm
library(ggrepel)
library(scales)

library(ggplot2)
#library(ggExtra) #https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html

library(cowplot)
library(ggridges)

# library(limma)
# library(writexl)


# install.packages("ggforce")
library(ggforce)
getwd()

# read the variable from the text file
cwd <- readLines("Variable_Storage/folder_path.txt")[1]
cwd
setwd(cwd)



# setwd(cwd)




file_list = list.files(path="results_burda", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)
# 
# 
# 
# # Function to plot the graph for each file
# plot_significant_lipids_1 <- function(df, title_for_plot) {
#   
#   # Filtering based on FDR and categorizing by up or down based on logFC
#   summarized_df <- df %>%
#     filter(FDR < 0.1) %>%
#     mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
#     group_by(type, Direction) %>%
#     summarize(count = n(), .groups = "keep")
#   
#   # Plotting the bar graph
#   ggplot(summarized_df, aes(x = type, y = count, fill = Direction)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     labs(title = title_for_plot, y = "Number of significant lipids") +
#     theme_minimal() +
#     scale_fill_manual(values = c("Up" = "blue", "Down" = "red"))
#   ggsave(filename = paste0("plots/count_", title_for_plot, ".svg"))
#   
# }
# 
# 
# 
# # Function to plot the graph for each file
# plot_significant_lipids <- function(df, title_for_plot) {
#   
#   # Filtering based on FDR and categorizing by up or down based on logFC
#   summarized_df <- df %>%
#     filter(FDR < 0.1) %>%
#     mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
#     group_by(type, Direction) %>%
#     summarize(count = n(), .groups = "keep")
#   
#   # Changing the count for "Down" direction to negative
#   summarized_df$count[summarized_df$Direction == "Down"] <- -summarized_df$count[summarized_df$Direction == "Down"]
#   
#   # Plotting the horizontal bar graph
#   ggplot(summarized_df, aes(y = type, x = count, fill = Direction)) +
#     geom_bar(stat = "identity", position = "stack") +
#     labs(title = title_for_plot, x = "Number of significant lipids") +
#     theme_minimal() +
#     scale_fill_manual(values = c("Up" = "blue", "Down" = "red")) +
#     coord_flip()
#   
#   ggsave(filename = paste0("plots/count_", title_for_plot, ".svg"))
# }
# Function to plot the graph for each file
plot_significant_lipids <- function(df, title_for_plot, bar_width = 0.5) {
  
  # Filtering based on FDR and categorizing by up or down based on logFC
  summarized_df <- df %>%
    filter(FDR < 0.1) %>%
    mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
    group_by(type, Direction) %>%
    summarize(count = n(), .groups = "keep")
  
  # Changing the count for "Down" direction to negative
  summarized_df$count[summarized_df$Direction == "Down"] <- -summarized_df$count[summarized_df$Direction == "Down"]
  
  # Plotting the horizontal bar graph with specified bar width
  plot <- ggplot(summarized_df, aes(y = type, x = count, fill = Direction)) +
    geom_bar(stat = "identity", position = "stack", width = bar_width) +
    labs(title = title_for_plot, x = "Number of significant lipids") +
    theme_minimal() +
    scale_fill_manual(values = c("Up" = "#71889c", "Down" = "#c5d2db"))
  
  print(plot)
  
  ggsave(filename = paste0("plots/count_", title_for_plot, ".svg"), plot = plot)
}






plot_combined_values <- function(df, Title1, Title2) {
  
  # Extract column lengths
  len1 <- as.numeric(df$Length1[1])
  len2 <- as.numeric(df$Length2[1])
  
  # Columns after 'type' are the ones of interest
  all_cols <- colnames(df)
  start_idx <- which(all_cols == "type") + 1
  cols1 <- all_cols[start_idx:(start_idx + len1 - 1)]
  cols2 <- all_cols[(start_idx + len1):(start_idx + len1 + len2 - 1)]
  
  # Compute the means for the groups of columns
  df <- df %>%
    mutate(mean1 = rowMeans(select(., one_of(cols1))),
           mean2 = rowMeans(select(., one_of(cols2))))
  
  # Sum by type for significant values
  df_sig <- df %>%
    filter(FDR < 0.1) %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  df_all <- df %>%
    group_by(type) %>%
    summarise(sum_mean1 = sum(mean1), 
              sum_mean2 = sum(mean2), .groups = "keep")
  
  # Combine the two dataframes to have a 'group' column and a 'value' column for plotting
  df_sig_long <- df_sig %>% 
    select(type, sum_mean1, sum_mean2) %>% 
    gather(key = "group", value = "value", -type)
  
  df_all_long <- df_all %>% 
    select(type, sum_mean1, sum_mean2) %>% 
    gather(key = "group", value = "value", -type)
  
  # Plotting significant values
  plot_sig <- ggplot(df_sig_long, aes(x = type, y = value, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste(Title1, "vs", Title2, " (Significant)"), y = "Sum of Means") +
    theme_minimal() +
    scale_fill_manual(values = c(sum_mean1 = "#71889c", sum_mean2 = "#c5d2db"), name = "Group", labels = c(Title1, Title2))
  
  # Plotting all values
  plot_all <- ggplot(df_all_long, aes(x = type, y = value, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste(Title1, "vs", Title2, " (All values)"), y = "Sum of Means") +
    theme_minimal() +
    scale_fill_manual(values = c(sum_mean1 = "#71889c", sum_mean2 = "#c5d2db"), name = "Group", labels = c(Title1, Title2))
  
  # Saving the plots
  ggsave(filename = paste0("plots/", Title1, "_vs_", Title2, "SUM_sig.svg"), plot = plot_sig)
  ggsave(filename = paste0("plots/", Title1, "_vs_", Title2, "SUM_all.svg"), plot = plot_all)
}



jj <-file_list[1]
file_list[2]
for (jj in file_list){
  excel_file <- read_csv(paste0("results_burda/",jj,sep=""))
  dir.create("plots", F)
  # title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub("[:| ]", "_", excel_file$Title_1[1])
  Title2 <- gsub("[:| ]", "_", excel_file$Title_2[1])
  Title1 <- gsub("__", "_", Title1)
  Title2 <- gsub("__", "_", Title2) 
  # Title1 <-
  
    
  title_for_plot <- paste0(Title1,Title2,sep="_")
  plot_combined_values(excel_file, Title1, Title2)
  # Plotting and saving
  plot_object <- plot_significant_lipids(excel_file, title_for_plot)
  
  # Ensuring the plots directory exists
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  

}
  
  
  # 
  # 
  # make_heatmap <- function(tp, design_mat, gr, contrasts, title, file) {
  #   
  #   DElist <-
  #     tp %>%
  #     get_DE_lipids(design_mat, gr, contrasts)
  #   
  #   if(length(DElist) == 0) {
  #     cat("No significant lipids for ", title)
  #     return()
  #   }
  #   
  #   if(length(DElist) == 1) {
  #     cat("Single significant lipids for ", title, " is ", DElist[1])
  #     return()
  #   }
  #   
  #   Blank <- log2(tp[[blank_name]])
  #   
  #   tp %>%
  #     mutate(lipid = make.unique(lipid)) %>%
  #     filter(lipid %in% DElist) %>%
  #     select( -type) %>%
  #     mutate_if(is.numeric, log2) %>%
  #     mutate_if(is.numeric, list(~ . - Blank)) %>%
  #     select(- blank_name) %>% 
  #     column_to_rownames("lipid") %>%
  #     as.matrix() %>%
  #     pheatmap::pheatmap(main = title,cluster_cols = none,
  #                        cluster_rows = none, filename = file,
  #                        cellheight = 10)
  # }
  # 
  # cells_lipid_expr %>%
  #   make_heatmap(design_expr, gr_expr, contrasts_expr, title_for_plot, 
  #                file = paste("plots/Heatmap_",title_for_plot,".svg",sep=''))
# }




