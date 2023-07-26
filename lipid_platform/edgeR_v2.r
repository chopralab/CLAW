# Loading required libraries for data analysis and visualization
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
library(plyr)
#library(writexl)

# Getting the current working directory
current_working_dir <- getwd()

# Listing all files in the current working directory
list.files()

# Changing the current working directory to the "lipid_platform" inside the "github/lipids/Lipidomics" folder
setwd("github/lipids/Lipidomics/lipid_platform")

# Reading the first line from the text file "varname3.txt" located inside the "Variable_Storage" folder
custom_working_dir <- readLines("Variable_Storage/varname3.txt")[1]

# Setting the current working directory to the value read from the text file
setwd(custom_working_dir)

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

  # # Removing unnecessary columns from the data file
  # lipid_data_file <- select(lipid_data_file, -c(title1, title2, plot_title, length1, length2,blank_name))
  # Removing unnecessary columns from the data file
  lipid_data_file <- select(lipid_data_file, -c("Title1", "Title2", "Title", "length1", "length2","Blank_name"))

  # Renaming the "Lipid" column to "lipid" and "Class" column to "type"
  # lipid_data_file <- lipid_data_file %>%
  #   rename(lipid = `Lipid`) %>%
  #   rename(type= `Class`)
  # lipid_data_file <- rename(lipid_data_file, type = Class)
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

  
  # cells_lipid_expr<-cells_lipid_expr %>%
  #   mutate(type = ifelse((type == "PCandSM" & grepl('PC', lipid)), "PC", 
  #                        ifelse((type == "PCandSM" & grepl('SM', lipid)), "SM", 
  #                               ifelse(type=="PCandSM", "SM", paste0(type)))))
  # 
  # cells_lipid_expr <- cells_lipid_expr %>%
  #   mutate(type = ifelse(type %in% c("TAG1", "TAG2"), "TAG", type))
  # 
  # cl_e1_tbl <-
  #   cells_lipid_expr %>%
  #   experiment_helper
  # cl_e1_tbl
  # Combine edgeR results with the original lipid expression data

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
  
  # [Perform necessary preprocessing steps on 'lipid_data_file']

  # Merge results and calculate 'merged_results'
  # [Perform necessary steps to compute 'merged_results']

  # Generate a title for the plot
  title_for_plot <- gsub(":", "_", lipid_data_file$Title[1])  # Or some other appropriate line

  # Create the plot
  p <- merged_results %>%
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
  ggsave(paste("plots/test3/Ridge_Plot_All_Lipids_", title_for_plot, ".pdf", sep=''), plot = p)
}


#####
# Create a ridge plot of fold changes
merged_results %>%
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
ggsave(paste("plots/test3/Ridge_Plot_All_Lipids_", title_for_plot, ".pdf", sep=''))
#####
# PCA Plot
library(ggbiplot)

merged_results_numeric <- merged_results[, sapply(merged_results, is.numeric)]  # Select only numeric columns for PCA
merged_results_numeric <- na.omit(merged_results_numeric)  # Remove NA values

str(merged_results)

pca_result <- prcomp(merged_results_numeric, center = TRUE, scale = TRUE)  # Perform PCA

# Make PCA biplot
pca_biplot <- ggbiplot(pca_result, obs.scale = 1, var.scale = 1,
                       groups = merged_results$type, ellipse = TRUE, circle = TRUE)
pca_biplot <- pca_biplot + scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

print(pca_biplot)

# Save the plot
ggsave(paste("plots/test/PCA_Plot_", title_for_plot, ".pdf", sep=''))

# Heatmap Plot
library(pheatmap)

# Normalize the data for better visualization in heatmap
merged_results_normalized <- t(scale(t(merged_results_numeric)))

# Make a Heatmap, assuming 'type' as factor of interest.
pheatmap(merged_results_normalized, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "none", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         show_rownames = FALSE, 
         annotation_col = as.data.frame(merged_results$type),
         annotation_colors = list(type = lipid_class_colors),
         fontsize = 8)

# Save the plot
ggsave(paste("plots/test/Heatmap_", title_for_plot, ".pdf", sep=''))

# PCA Plot
library(ggplot2)

# Perform PCA
pca_result <- prcomp(merged_results_numeric, center = TRUE, scale = TRUE)

# Create a data frame for plotting
plot_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  group = merged_results$type
)

# Create PCA plot
pca_plot <- ggplot(data = plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.8) +
  labs(color = "Lipid class") +
  ggtitle("PCA Plot") +
  theme_minimal() +
  scale_color_manual(values = lipid_class_colors)

print(pca_plot)

# Save the plot
ggsave(paste("plots/test/PCA_Plot_", title_for_plot, ".pdf", sep=''))


# Select only numeric columns for PCA
merged_results_numeric <- merged_results[sapply(merged_results, is.numeric)]

# Now try to run PCA
pca_result <- prcomp(merged_results_numeric, center = TRUE, scale = TRUE)

# Print the PCA result
print(pca_result)

# Perform PCA
pca_result <- prcomp(merged_results_numeric, center = TRUE, scale = TRUE)
print(pca_result)

  
  # xx %>%
  #   
  #   filter(FDR < 0.10) %>% #filter by sig
  #   ggplot(aes(x=logFC, y = type, fill = type)) +
  #   geom_density_ridges2(alpha = 0.5, size = .5) + 
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   theme_classic() +
  #   ggtitle(title_for_plot)+
  #   xlab("Fold change, lipids FDR<0.10") +
  #   ylab("") +#xlim(-2, 4)+
  #   scale_alpha(guide = 'none') +
  #   
  #   scale_fill_discrete(name = "Lipid class", guide = 'none') +
  #   scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
  #   scale_y_discrete(limits = rev) #+
  # # facet_wrap(~comparison_LFC)
  # 
  # ggsave(paste("plots/Ridge_Plot_Filtered_FDR_",title_for_plot,".pdf",sep=''))
  
  # make_volcano_plot <- function(df, title) {
  #   df %>%
  #     mutate(sig = factor(FDR < 0.10)) %>%
  #     ggplot(aes(logFC, -log10(FDR), color = sig)) +
  #     geom_point() +
  #     scale_color_manual(values = c("none" = "black", "TRUE" = "red")) +
  #     guides(color = F) +
  #     ggtitle(title)
  # }
  # 
  # cl_e1_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot,sep=''))
  # ggsave(paste("plots/Volcano_Plot_",title_for_plot,".png",sep=''))
  # 
  # 
  
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

  # source("ggbiplot.R")
  
  # 
  # get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.05) {
  #   dls <-
  #     counts %>%
  #     perform_analysis_raw(design_mat, gr) %>%
  #     calculate_significance(contrasts) %>%
  #     decideTestsDGE(p.value = p.value)
  #   
  #   rownames(dls)[dls %>% as.logical()]
  # }
  # # 
  # make_pca_plot <- function(tp, design_mat, gr, contrasts,
  #                           title = "PCA plot",
  #                           ellipse = T, var.axes = F,
  #                           labels = T) {
  # 
  #   print(labels)
  #   print(labels)
  #   tp_edger <-
  #     tp %>%
  #     get_DE_lipids(design_mat, gr, contrasts)
  # 
  #   if(length(tp_edger) == 0) {
  #     cat("No significant lipids for ", title)
  #     return()
  #   }
  # 
  #   if(length(tp_edger) == 1) {
  #     cat("Single significant lipids for ", title, " is ", tp_edger[1])
  #     return()
  #   }
  # 
  #   tp %>%
  #     na.omit %>%
  #     mutate(lipid = make.unique(lipid)) %>%
  #     filter(lipid %in% tp_edger) %>%
  #     select(-Transition, -type, -blank_name) %>%
  #     column_to_rownames("lipid") %>%
  #     as.matrix() %>%
  #     t %>%
  #     prcomp(center = T, scale = T) ->
  #     prcomp_data
  # 
  #   groups = NULL
  # 
  #   tp %>%
  #     select(-Transition, -type, -blank_name, -lipid) %>%
  #     colnames() ->
  #     labels.tmp
  # 
  #   # groups = substr(labels.tmp, 1, 4)
  #   groups =groups_PCA
  #   # labels = groups_PCA
  #   print(groups)
  #   # print(substr(labels.tmp, 1, 5))
  #   if (!is.null(labels)) {
  #     labels = labels.tmp
  #   }
  # 
  #   prcomp_data %>%
  #     ggbiplot(ellipse = ellipse,
  #              labels = labels,
  #              groups = groups,
  #              var.axes = var.axes
  #     ) +
  #     ggtitle(title) +
  #     cowplot::theme_cowplot()
  # }
  
  
  # source("ggbiplot.R")
  
 
  
  # cells_lipid_expr %>%
  #   make_pca_plot(design_expr, gr_expr, contrasts_expr, paste("PCA_Plot_",title_for_plot,sep=''))
  # ggsave(paste("plots/PCA_",title_for_plot,".svg",sep=''))
  
  # contrasts_expr
  # gr_expr
  # design_expr
  
  # 
  # 
  # # Assuming you have a dataframe called cells_lipid_expr and two vectors called group1 and group2 with the column names
  # library(factoextra)
  # # Combine group1 and group2 into a single vector
  # filtered_xx <- xx %>%
  #   filter(FDR < 0.1)
  # 
  # 
  # 
  # 
  # all_groups <- select(excel_file, -c(blank_name,lipid,type))
  # all_groups <- names(all_groups)
  # 
  # 
  # # Extract only the columns that belong to group1 and group2 from the dataframe
  # cells_lipid_expr_subset <- filtered_xx[, all_groups]
  # 
  # # Transpose the data frame
  # cells_lipid_expr_transposed <- t(cells_lipid_expr_subset)
  # 
  # # Perform PCA on the transposed dataframe
  # pca_result <- prcomp(cells_lipid_expr_transposed, center = TRUE, scale = TRUE)
  # 
  # # Create a dataframe for the groups and their corresponding titles
  # group_data <- data.frame(
  #   sample = all_groups,
  #   group = c(rep(Title1, length1), rep(Title2,length2))
  # )
  # 
  # # Calculate the PCA scores
  # pca_scores <- as.data.frame(pca_result$x[, 1:2])
  # 
  # # Create a new data frame with the PCA scores, sample names, and group information
  # plot_data <- data.frame(
  #   PC1 = pca_scores$PC1,
  #   PC2 = pca_scores$PC2,
  #   sample = rownames(pca_scores),
  #   group = group_data$group
  # )
  # 
  # # Create the PCA plot using ggplot2
  # pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
  #   geom_point(size = 3) +
  #   geom_text_repel(size = 3) +
  #   theme_minimal() +
  #   labs(color = "Groups") +
  #   ggtitle(title_for_plot) +
  #   xlab(paste0("PC1: ", round(pca_result$sdev[1] * 100 / sum(pca_result$sdev), 2), "% variance")) +
  #   ylab(paste0("PC2: ", round(pca_result$sdev[2] * 100 / sum(pca_result$sdev), 2), "% variance"))
  # 
  # # Print the PCA plot
  # print(pca_plot)
  # 
  # ggsave(paste("plots/PCA_",title_for_plot,".svg",sep=''))
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




