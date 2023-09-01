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
  
  title <- sub("_full.csv", "", basename(file_path))
  labels_title_1 <- unique(df$Title_1)
  labels_title_2 <- unique(df$Title_2)
  annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))

  # Function to save pheatmap as PDF
  save_pheatmap_pdf <- function(x, filename, width=4, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }

  # Setup conditions for FDR filtering
  conditions <- list(
    list(filter = TRUE, suffix = "_FDR"),
    list(filter = FALSE, suffix = "")
  )

  # Loop through conditions to generate and save heatmaps
  for (condition in conditions) {
    if (condition$filter) {
      df_filtered <- df %>% filter(FDR < 0.1)
    } else {
      df_filtered <- df
    }
    df_filtered <- df_filtered %>% arrange(type)
    
    title_suffix <- condition$suffix
    heatmap_breaks <- seq(-1, 1, length.out = 100)
    annotation_data <- data.frame(type = df_filtered$type)
    first_occurrence <- !duplicated(df_filtered$type)
    df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
    
    annotation_col_df <- data.frame(Labels = annotation_labels)
    rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
    
    heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                      main = paste0(title, title_suffix),
                                      cluster_rows = FALSE, 
                                      cluster_cols = FALSE, 
                                      show_colnames = TRUE, 
                                      show_rownames = TRUE,
                                      breaks = heatmap_breaks,
                                      border_color = "black", 
                                      fontsize = 10,
                                      fontsize_legend = 25,
                                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                      labels_row = df_filtered$Row_Label,
                                      annotation_col = annotation_col_df)

    # Save the heatmap
    save_pheatmap_pdf(heatmap_obj, paste0("heatmaps_large/", title, title_suffix, ".pdf"))
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










##Working function

file_path <- "Genotype_ cKO__time point_ 28d ISCI vs Genotype_ cKO__time point_ healthy8_9_23_full.csv"


df <- read.csv(file_path, stringsAsFactors = FALSE)

# Arrange dataframe (assuming you meant to arrange by a column, e.g., "lipid")
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

title <- sub("_full.csv", "", basename(file_path))
labels_title_1 <- unique(df$Title_1)
labels_title_2 <- unique(df$Title_2)
annotation_labels <- c(rep(labels_title_1, each=df$Length1[1]), rep(labels_title_2, each=df$Length2[1]))


# Setup the two conditions
conditions <- list(
  list(filter = TRUE, suffix = "_FDR"),
  list(filter = FALSE, suffix = "")
)

for (condition in conditions) {
  
  # Check condition for filtering
  if (condition$filter) {
    df_filtered <- df %>% filter(FDR < 0.1)
  } else {
    df_filtered <- df
  }
  
  title_suffix <- condition$suffix
  
  # Calculating range for heatmap color breaks
  heatmap_breaks <- seq(-1, 1, length.out = 100)
  
  # Extract the 'type' column from df_filtered and make it a dataframe
  annotation_data <- data.frame(type = df_filtered$type)
  
  # Find the last occurrence of each unique type
  first_occurrence <- !duplicated(df_filtered$type)
  
  # Change row names based on the first occurrence of each type
  df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
  
  # Assuming annotation_labels has the same length as columns of df_filtered[intensity_columns]
  annotation_col_df <- data.frame(Labels = annotation_labels)
  
  # Make sure the row names of the annotation dataframe match the column names of the data
  rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
  
  # Generate the heatmap
  heatmap_obj <- pheatmap::pheatmap(df_filtered[intensity_columns], 
                                    main = paste0(title, title_suffix),
                                    cluster_rows = FALSE, 
                                    cluster_cols = FALSE, 
                                    show_colnames = TRUE, 
                                    show_rownames = TRUE,
                                    breaks = heatmap_breaks,
                                    border_color = "black", 
                                    fontsize = 10,
                                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                    labels_row = df_filtered$Row_Label,  # Use the custom labels
                                    annotation_col = annotation_col_df  # Use the annotation dataframe
  )
  
  # Save the heatmap as PDF
  save_pheatmap_pdf(heatmap_obj,  paste0("heatmaps/", title, title_suffix, ".pdf"))
}

# The save_pheatmap_pdf function remains unchanged







# 
# # Filter df based on FDR
# df_filtered <- df %>% filter(FDR < 0.1)
# title_suffix <- "_FDR"
# 
# # Calculating range for heatmap color breaks
# heatmap_breaks <- seq(-1, 1, length.out = 100)
# 
# # Extract the 'type' column from df_filtered and make it a dataframe
# annotation_data <- data.frame(type = df_filtered$type)
# 
# 
# # Find the last occurrence of each unique type
# 
# first_occurrence <- !duplicated(df_filtered$type)
# 
# # Change row names based on the first occurrence of each type
# df_filtered$Row_Label <- ifelse(first_occurrence, df_filtered$type, "")
# ##No labels
# # pheatmap::pheatmap(df_filtered[intensity_columns],
# #                    main = paste0(title, title_suffix),
# #                    cluster_rows = FALSE,
# #                    cluster_cols = FALSE,
# #                    show_colnames = TRUE,
# #                    show_rownames = TRUE,
# #                    breaks = heatmap_breaks,
# #                    border_color = "black",
# #                    fontsize = 10,
# #                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
# #                    labels_row = df_filtered$Row_Label  # Use the custom labels
# # )
# 
# # Assuming annotation_labels has the same length as columns of df_filtered[intensity_columns]
# annotation_col_df <- data.frame(Labels = annotation_labels)
# 
# # Make sure the row names of the annotation dataframe match the column names of the data
# rownames(annotation_col_df) <- colnames(df_filtered[intensity_columns])
# 
# 
# xx <-pheatmap::pheatmap(df_filtered[intensity_columns], 
#                    main = paste0(title, title_suffix),
#                    cluster_rows = FALSE, 
#                    cluster_cols = FALSE, 
#                    show_colnames = TRUE, 
#                    show_rownames = TRUE,
#                    breaks = heatmap_breaks,
#                    border_color = "black", 
#                    fontsize = 10,
#                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
#                    labels_row = df_filtered$Row_Label,  # Use the custom labels
#                    annotation_col = annotation_col_df  # Use the annotation dataframe
# )
# 
# 
# 
# save_pheatmap_pdf <- function(x, filename, width=15, height=15) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# save_pheatmap_pdf(xx,  paste0("heatmaps/", title, title_suffix, ".pdf"))
# 



  # Save the plot
# ggsave(filename = paste0("heatmaps/", title, title_suffix, ".svg"), width = 10, height = 10)


# Check if heatmaps folder exists, if not, create it
if (!dir.exists("heatmaps")) {
  dir.create("heatmaps")
}

# Get the list of files with "_full.csv"
files <- list.files(path = "results", pattern = "_full.csv$", full.names = TRUE)

# Create heatmap for each file
lapply(files, create_heatmap)