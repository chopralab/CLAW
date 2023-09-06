library(tidyverse)
library(pheatmap)
library(RColorBrewer)


library(dplyr)
library(RColorBrewer)
# File path
create_heatmap <- function(file_path) {
  
  # Read and preprocess the data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Arrange dataframe
  df <- df %>%
    arrange(lipid)
  
  # Identify intensity columns
  intensity_columns <- setdiff(names(df), c("lipid", "logFC", "logCPM", "LR", "PValue", "FDR", "Length1", "Length2", "Title_1", "Title_2", "type", tail(names(df), 1)))
  
  # Apply transformations
  df <- df %>%
    mutate(across(all_of(intensity_columns), log2)) %>%
    rowwise() %>%
    mutate(mean_intensity = mean(c_across(all_of(intensity_columns)), na.rm = TRUE)) %>%
    mutate(across(all_of(intensity_columns), ~ . - mean_intensity)) %>%
    select(-mean_intensity, -tail(names(df), 1))
  
  # Filter by FDR < 0.1
  df <- df %>% filter(FDR < 0.1)
  
  # Create annotation labels
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))
  
  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=15, height=15) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Extract the unique types
  unique_types <- unique(df$type)
  
  for (current_type in unique_types) {
    df_filtered <- df %>% filter(type == current_type) %>% arrange(lipid)
    
    # Check if only 1 or 0 lipid is present
    if(nrow(df_filtered) <= 1) next
    
    title_prefix <- current_type
    title <- sub("_full.csv", "", basename(file_path))
    heatmap_breaks <- seq(-1, 1, length.out = 100)
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
    heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                      main = paste0(title_prefix, " - ", title),
                                      cluster_rows = FALSE, 
                                      cluster_cols = FALSE, 
                                      show_colnames = TRUE, 
                                      show_rownames = TRUE,
                                      breaks = heatmap_breaks,
                                      border_color = "black", 
                                      fontsize = 10,
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                      labels_row = df_filtered$lipid,
                                      annotation_col = annotation_col_df)
    
    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps/", title_prefix, "_", title, ".pdf"))
  }
}



# You can now use the function as:
# create_heatmap("path_to_your_file.csv")

# Check if heatmaps folder exists, if not, create it
if (!dir.exists("heatmaps")) {
  dir.create("heatmaps")
}

# Get the list of files with "_full.csv"
files <- list.files(path = "results", pattern = "_full.csv$", full.names = TRUE)

# Create heatmap for each file
lapply(files, create_heatmap)


