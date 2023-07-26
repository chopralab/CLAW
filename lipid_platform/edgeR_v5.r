# Loading required libraries for data analysis and visualization
#library(writexl)
library(ggfortify)
library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
library(ggsignif)
library(stringr)
library(gtools)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(cowplot)
library(ggExtra) #https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html
library(limma)
library(reshape2)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(dplyr)
library(ggplot2movies) 


# Getting the current working directory
current_working_dir <- getwd()
current_working_dir

CLAW_working_dir <- readLines("CLAW_config.txt")[1]
setwd(CLAW_working_dir)
CLAW_working_dir

# Reading the first line from the text file "varname3.txt" located inside the "Variable_Storage" folder
project_working_dir <- readLines("project_path/brain_5xFAD_old.txt")[1]
setwd(project_working_dir)
project_working_dir
list.files()

# Listing all files in the "Pre_EdgeR" directory
files_in_pre_edgeR <- list.files(path="Pre_EdgeR")

# Looping through each file in the "Pre_EdgeR" directory
for (current_file in files_in_pre_edgeR){

  # Reading the CSV file and storing it in a variable
  lipid_data_file <- read_csv(paste0("Pre_EdgeR/",current_file))

  # Creating a directory named "plots" if it doesn't exist
  dir.create("plots", showWarnings = FALSE)

  # Formatting the title for the plot by replacing ":" with "_"
  plot_title <- gsub(":", "_",lipid_data_file$Title[1])

  # Formatting other titles and lengths
  title_for_plot <- gsub(":", "_",lipid_data_file$Title[1])
  title1 <- gsub(":", "_",lipid_data_file$Title1[1])
  title2 <- gsub(":", "_",lipid_data_file$Title2[1])
  blank_name <- lipid_data_file$Blank_name[1]
  length1 <- lipid_data_file$length1[1]
  length2 <- lipid_data_file$length2[1]

  # Removing unnecessary columns from the data file
  lipid_data_file <- select(lipid_data_file, -c("Title1", "Title2", "Title", "length1", "length2","Blank_name"))
  lipid_data_file <- dplyr::rename(lipid_data_file, type = Class)


  # Checking if length1 or length2 is equal to 1, if yes then skip the current iteration
  if (length1 == 1 | length2 == 1) {
    next
  }

  # Creating two character vectors of lengths 'length1' and 'length2' containing "GR1" and "GR2" respectively
  group1 <- c(rep("GR1",length1))
  group2 <- c(rep("GR2",length2))

  # Combining the two character vectors to create a single vector for PCA groups
  PCA_groups <- c((rep(title1,length1)),(rep(title2,length2)))

  # Storing the lipid expression data
  lipid_expression_data <- lipid_data_file
  lipid_expression_data <- dplyr::rename(lipid_expression_data, lipid = Lipid)

}

# Create factor levels for edgeR analysis
experiment_groups <- c(group1, group2, blank_name) %>%
  factor(levels = c(blank_name, "GR1", "GR2"))

# Create a design matrix for the model using experiment groups
design_matrix <- model.matrix(~experiment_groups)

# Define contrasts for comparison in the edgeR analysis
contrast_matrix <- makeContrasts(
  group_difference = experiment_groupsGR1 - experiment_groupsGR2,
  levels = design_matrix
)

# Function to perform initial edgeR analysis 
perform_edgeR_analysis <- function(counts, design_matrix, experiment_groups) {
  # Preprocess counts data
  edgeR_data <- DGEList(
    counts = counts %>%
      na.omit %>%
      mutate(lipid = make.unique(lipid)) %>%
      select(-type) %>%
      column_to_rownames("lipid"),
    group = experiment_groups
  )

  # Normalize factors and estimate common dispersion
  edgeR_data <- calcNormFactors(edgeR_data, method="TMM")
  edgeR_data <- estimateCommonDisp(edgeR_data, design=design_matrix)
  
  edgeR_data
}

# Function to fit generalized linear model and compute likelihood ratio tests
calculate_significance <- function(edgeR_data, contrast_matrix) {
  edgeR_data %>%
    glmFit() %>%
    glmLRT(contrast = contrast_matrix)
}

# Function to streamline edgeR analysis
edgeR_pipeline <- function(df) {
  df %>%
    perform_edgeR_analysis(design_matrix, experiment_groups) %>%
    calculate_significance(contrast_matrix) %>%
    topTags(1500000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    as_tibble()
}

# Perform the edgeR analysis and store the results
edgeR_results <- lipid_expression_data %>% edgeR_pipeline
edgeR_results

# Display the names of the columns in the lipid expression data
names(lipid_expression_data)

merged_results <- merge(lipid_expression_data, edgeR_results)

# Define lipid classes and their corresponding colors for visualization
lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')

# Create a named vector to map lipid classes to their colors
lipid_class_colors <- setNames(lipid_colors, lipid_classes)

# Assuming 'files_in_pre_edgeR' is your list of csv files
for (current_file in files_in_pre_edgeR) {
  
  # Load and preprocess the data file
  lipid_data_file <- read_csv(paste0("Pre_EdgeR/", current_file))
  
  # Generate a title for the plot
  title_for_plot <- gsub(":", "_", lipid_data_file$Title[1])  # Or some other appropriate line

####
# RIDGE PLOT
####
if (!dir.exists("plots/ridge_plot")) {
  dir.create("plots/ridge_plot", recursive = TRUE, showWarnings = FALSE)
}
  # Create the plot
  ridge_plot <- merged_results %>%
    ggplot(aes(x=logFC, y = type, fill = type)) +
    geom_density_ridges2(alpha = 0.5, size = .5) +  # Create ridges
    geom_vline(xintercept = 0, linetype = "dashed") +  # Add vertical line at x=0
    theme_classic() +  # Use classic theme
    ggtitle(title_for_plot) +  # Add title
    xlab("Fold change, lipids all") +  # Label x-axis
    ylab("") +  # Remove y-axis label
    scale_alpha(guide = 'none') +  # Remove alpha legend
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +  # Fill colors according to lipid classes
    scale_y_discrete(limits = rev)  # Reverse y-axis limits

  # Save the plot to a pdf file
  ggsave(paste("plots/ridge_plot/Ridge_Plot_All_Lipids_", title_for_plot, ".pdf", sep=''), plot = ridge_plot)
}

####
# PCA PLOT
####

# Select only numeric columns for PCA
numeric_data <- lipid_data_file[sapply(lipid_data_file, is.numeric)]

# Get the standard deviation of each column
column_sd <- apply(numeric_data, 2, sd)

# Identify columns where SD is not 0
variable_columns <- column_sd != 0
# Subset data to keep only variable columns
numeric_data <- numeric_data[, variable_columns]

# Use prcomp for PCA
pca_result <- prcomp(numeric_data, scale. = TRUE)

# Use autoplot to automatically generate a ggplot2 plot
pca_plot <- autoplot(pca_result, 
                         data = lipid_data_file, 
                         #colour = 'type',  # change 'type' to your group column name
                         #aes(color = type)
                         label = TRUE, 
                         label.size = 3, 
                         #size = 'type',  # change 'type' to a column name to size points by
                         ellipse = TRUE, 
                         ellipse.type = 't',
                         alpha = 0.6)

# print the plot
print(pca_plot)
#####

# Create the "plots/pca" directory if it doesn't exist
dir.create("plots/pca", showWarnings = FALSE)

# Loop through each file in the "Pre_EdgeR" directory
for (current_file in files_in_pre_edgeR) {
  
  # Load and preprocess the data file
  lipid_data_file <- read_csv(paste0("Pre_EdgeR/", current_file))
  
  # [Perform necessary preprocessing steps on 'lipid_data_file']

  # Select only numeric columns for PCA
  numeric_data <- lipid_data_file[sapply(lipid_data_file, is.numeric)]

  # Get the standard deviation of each column
  column_sd <- apply(numeric_data, 2, sd)

  # Identify columns where SD is not 0
  variable_columns <- column_sd != 0
  # Subset data to keep only variable columns
  numeric_data <- numeric_data[, variable_columns]

  # Use prcomp for PCA
  pca_result <- prcomp(numeric_data, scale. = TRUE)

  # Generate a title for the plot
  title_for_plot <- gsub(":", "_", lipid_data_file$Title[1])  # Or some other appropriate line

  # Use autoplot to automatically generate a ggplot2 plot
  pca_plot <- autoplot(pca_result, 
                       data = lipid_data_file, 
                       #colour = 'type',  # change 'type' to your group column name
                       label = TRUE, 
                       label.size = 3, 
                       #size = 'type',  # change 'type' to a column name to size points by
                       ellipse = TRUE, 
                       ellipse.type = 't',
                       alpha = 0.6)

  # Generate the file path for saving the PCA plot
  file_path <- paste("plots/pca/PCA_Plot_", title_for_plot, ".pdf", sep = '')
  
  # Print the file path for debugging
  print(paste("Saving plot to:", file_path))
  
  # Save the PCA plot to a PDF file
  ggsave(filename = file_path, plot = pca_plot)
}

####
# HEATMAP Plot
####

# Create the "plots/heatmap" directory if it doesn't exist
dir.create("plots/heatmap", showWarnings = FALSE)
# Loop through each file in the "Pre_EdgeR" directory
for (current_file in files_in_pre_edgeR) {
  
  # Load and preprocess the data file
  lipid_data_file <- read_csv(paste0("Pre_EdgeR/", current_file))
  
  # [Perform necessary preprocessing steps on 'lipid_data_file']

  # Compute the correlation matrix
  cor_matrix <- cor(numeric_data)

  # Create a heatmap of the correlation matrix
  heatmap_plot <- pheatmap(cor_matrix)

  # Generate a title for the plot
  title_for_plot <- gsub(":", "_", lipid_data_file$Title[1])  # Or some other appropriate line

  # Create a folder named "plots/heatmap" if it doesn't exist
  dir.create("plots/heatmap", showWarnings = FALSE)

  # Save the heatmap plot to a PDF file in the "plots/heatmap" folder
  pdf(file = paste("plots/heatmap/Heatmap_", title_for_plot, ".pdf", sep = ''))
  print(heatmap_plot)  # Plot the heatmap to the PDF file
  dev.off()  # Close the current graphics device
}


####
# Summary and Full Results
####

  # Define a function that generates summary and full result files
write_summary_and_results <- function(edgeR_table, lipid_data, experiment_name) {
    
    # Merge edgeR results with the original lipid data, sort by FDR
    edgeR_table %>%
      merge(lipid_data) %>%
      as_tibble() %>%
      arrange(FDR) -> results
    
    # Write the full results to a CSV file in the 'results' directory
    write_csv(results, paste0("results/", experiment_name, "_full.csv"))
    
    # Summarise the number of up- and down-regulated lipids per type
    results %>%
      group_by(type) %>%
      summarise(Down = sum(logFC < 0 & FDR < 0.1),
                Up = sum(logFC > 0 & FDR < 0.1)) %>%
      
    # Write the summary results to a CSV file in the 'results' directory
      write_csv(paste0("results/", experiment_name, "_summary.csv"))
}

# Create 'results' directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Apply the function to the edgeR results
edgeR_results %>% write_summary_and_results(lipid_expression_data, title_for_plot)

# The tables can be printed for verification purposes
edgeR_results



lipid_expression_data




