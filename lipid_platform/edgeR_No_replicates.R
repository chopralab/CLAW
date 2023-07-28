


library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
library(tidyverse)


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

# read the variable from the text file
cwd <- readLines("Variable_Storage/folder_path.txt")[1]
cwd
setwd(cwd)


file_list = list.files(path="Pre_EdgeR", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)


for (jj in file_list){
  excel_file <- read_csv(paste0("Pre_EdgeR/",jj,sep=""))
  dir.create("plots", F)
  title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub(":", "_",excel_file$Title1[1])
  Title2 <- gsub(":", "_",excel_file$Title2[1])
  
  blank_name <- excel_file$Blank_name[1]
  length1 <- excel_file$length1[1]
  length2 <- excel_file$length2[1]
  
  lengther1 = length1
  lengther2 = length2
  
  
  
  
  
  excel_file <- select(excel_file, -c(Title1, Title2, Title, length1, length2,Blank_name))
  
  ##Change Class to Type and change Lipid to lipid
  excel_file <- excel_file %>%
    rename(lipid = `Lipid`) %>%
    rename(type= `Class`)
  
  
  if (isTRUE(length1 > 1 & length2 > 1)) {
    next
  }

  # Rest of the loop code

  
  gr1 = c(rep("GR1",length1))
  gr2 = c(rep("GR2",length2))
  
  
  groups_PCA = c((rep(Title1,length1)),(rep(Title2,length2)))
  
  groups_PCA
  
  
  cells_lipid_expr <- excel_file
  
  
  # Create a dataframe with Length1 and Length2 columns
  df2_other <- data.frame(Length1 = rep(length1, dim(cells_lipid_expr)[1]), 
                          Length2 = rep(length2, dim(cells_lipid_expr)[1]),
                          Title_1 = rep(Title1, dim(cells_lipid_expr)[1]),
                          Title_2 = rep(Title2, dim(cells_lipid_expr)[1]))
  
  ###Need to auto column names - Figure out blank its only for 1 of them not all 3
  ##Including blank
  counts <- cells_lipid_expr
  counts$lipid <- NULL
  # counts$Transition <- NULL
  counts$type <- NULL
  
  counts_matrix <- as.matrix(counts)
  counts <- as.matrix(counts)
  
  rownames(counts) <- cells_lipid_expr$lipid

  ####Creating DGEList Object
  ###New Code
  group <- c(rep("group1", length1), rep("group2", length2), "Blank")
  y <- DGEList(counts=counts, group=factor(group))
  y$common.dispersion <- 0.1
  # Creating the design matrix
  design <- model.matrix(~ 0 + y$samples$group + 1) # Add the intercept term back
  colnames(design) <- c("Intercept", "group1", "group2") # Rename the columns accordingly
  
  
  ##Old code with only 1 1 for intercept
  # y <- DGEList(counts=counts, group=factor(c(rep("group1", length(group1)), rep("group2", length(group2)), "Blank")))
  # ###Manually setting Common Dispersion
  # y$common.dispersion <- 0.1
  # ###Creating Experiment Design Matrix and using Blank as Intercept.
  # design <- model.matrix(~ 0 + y$samples$group) # Use Blank as an intercept
  # colnames(design) <- c("Blank", "group1", "group2") # Use Blank as an intercept
  
  ###Fitting the GLM
  fit <- glmFit(y, design)
  
  ##Creating Contrasts
  lrt <- glmLRT(fit, contrast=makeContrasts(group1 - group2, levels=design))
  lrt_results <- topTags(lrt, n=Inf)
  
  # Convert the lrt_results data frame to a standard data frame
  lrt_results_df <- as.data.frame(lrt_results)
  
  # Add a "lipid" column to the lrt_results_df data frame by extracting the row names
  lrt_results_df$lipid <- rownames(lrt_results_df)
  
  # Merge the lrt_results_df data frame with the cells_lipid_expr data frame, matching the "lipid" column
  merged_data2 <- merge(cells_lipid_expr, lrt_results_df, by="lipid")
  merged_data2
  
  cl_e1_tbl<-merged_data2
  xx <-merged_data2
  
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  
  dir.create("plots", F)
  
  xx %>%
    
    ggplot(aes(x=logFC, y = type, fill = type)) +
    geom_density_ridges2(alpha = 0.5, size = .5) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    xlab("Fold change, lipids all") +
    ylab("") +#xlim(-2, 4)+
    scale_alpha(guide = 'none') +
    scale_fill_discrete(name = "Lipid class", guide = 'none') +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
    scale_y_discrete(limits = rev) #+
  ggsave(paste("plots/Ridge_Plot_All Lipids_",title_for_plot,".pdf",sep=''))
  
  # xx %>%
  #   
  #   filter(FDR < 0.10) %>% #filter by sig
  #   ggplot(aes(x=logFC, y = type, fill = type)) +
  #   geom_density_ridges2(alpha = 0.5, size = .5) + 
  #   geom_vline(xintercept = 0, linetype = "dashed") +
  #   theme_classic() +
  #   xlab("Fold change, lipids FDR<0.10") +
  #   ylab("") +#xlim(-2, 4)+
  #   scale_alpha(guide = 'none') +
  #   
  #   scale_fill_discrete(name = "Lipid class", guide = 'none') +
  #   scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
  #   scale_y_discrete(limits = rev) #+
  # # facet_wrap(~comparison_LFC)
  # 
  # ggsave(paste("plots/Ridge_Plot_Filtered_LogFC_",title_for_plot,".pdf",sep=''))
  
  
  write_summary_and_results <- function(tbl, df,df2, name) {
    df<-cbind(df2, df)
    tbl %>%
      # merge(df2) %>%
      merge(df) %>%
      as_tibble() %>%
      arrange(FDR) -> results
    
    write_csv(results, paste0("results/", name, "_full.csv"))
    
    results %>%
      group_by(type) %>%
      summarise(Down = sum(logFC < 0 & FDR < 0.1),
                Up = sum(logFC > 0 & FDR < 0.1)) %>%
      write_csv(paste0("results/", name, "_summary.csv"))
  }
  
  dir.create("results", F)
  
  cl_e1_tbl %>% write_summary_and_results(cells_lipid_expr,df2_other, title_for_plot)
  
  cl_e1_tbl
  source("ggbiplot.R")
  cl_e1_tbl
  
  make_volcano_plot <- function(df, title) {
    df %>%
      mutate(sig = factor(logFC < -1 | logFC > 1)) %>%
      ggplot(aes(logFC, -log10(FDR), color = sig)) +
      geom_point() +
      scale_color_manual(values = c("none" = "black", "TRUE" = "red")) +
      guides(color = F) +
      ggtitle(title)
  }
  
  cl_e1_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot,sep=''))
  ggsave(paste("plots/Volcano_Plot_",title_for_plot,".png",sep=''))
  
  # 
  # 
  # # Assuming you have a dataframe called cells_lipid_expr and two vectors called group1 and group2 with the column names
  # library(factoextra)
  # # Combine group1 and group2 into a single vector
  # # Combine group1 and group2 into a single vector
  # filtered_xx <- xx %>%
  #   filter(logFC < -1 | logFC > 1)
  # # Combine group1 and group2 into a single vector
  # 
  # 
  # all_groups <- select(excel_file, -c(blank_name,lipid,type))
  # all_groups <- names(all_groups)
  # # Combine group1 and group2 into a single vector
  # # all_groups <- c(group1, group2)
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

  

  

  
  
}
