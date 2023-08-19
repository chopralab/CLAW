
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

# library(limma)
# library(writexl)



getwd()

# read the variable from the text file
cwd <- readLines("Variable_Storage/folder_path.txt")[1]
cwd
setwd(cwd)



# setwd(cwd)




file_list = list.files(path="Pre_EdgeR", pattern=NULL, all.files=FALSE,
                       full.names=FALSE)

jj <-file_list[1]
file_list[2]
for (jj in file_list){
  excel_file <- read_csv(paste0("Pre_EdgeR/",jj,sep=""))
  dir.create("plots", F)
  title_for_plot <- gsub(":", "_",excel_file$Title[1])
  Title1 <- gsub(":", "_",excel_file$Title1[1])
  Title2 <- gsub(":", "_",excel_file$Title2[1])
  
  blank_name <- excel_file$Blank_name[1]
  length1 <- excel_file$length1[1]
  length2 <- excel_file$length2[1]
  
 
  
  
  excel_file <- select(excel_file, -c(Title1, Title2, Title, length1, length2,Blank_name))
  
  ##Change Class to Type and change Lipid to lipid
  excel_file <- excel_file %>%
    rename(lipid = `Lipid`) %>%
    rename(type= `Class`)
  
  if (length1 == 1 | length2 == 1) {
    next
  }
  
  
  
  gr1 = c(rep("GR1",length1))
  gr2 = c(rep("GR2",length2))
  
  
  groups_PCA = c((rep(Title1,length1)),(rep(Title2,length2)))
  
  groups_PCA
  
  
  
  
  blank_name
  
  
  
  
  
  cells_lipid_expr <- excel_file
  
  # Create a dataframe with Length1 and Length2 columns
  df2_other <- data.frame(Length1 = rep(length1, dim(cells_lipid_expr)[1]), 
                          Length2 = rep(length2, dim(cells_lipid_expr)[1]),
                          Title_1 = rep(Title1, dim(cells_lipid_expr)[1]),
                          Title_2 = rep(Title2, dim(cells_lipid_expr)[1]))
  
  #EdgeR groups
  gr_expr = c(gr1,gr2,
              blank_name) %>%
    factor(levels = c(blank_name, "GR1", "GR2"))
  design_expr = model.matrix(~gr_expr)
  
  contrasts_expr = makeContrasts(
    H = gr_exprGR1 - gr_exprGR2,
    levels = design_expr
  )
  
  
  
  ###Functions to do edgeR analysis
  
  perform_analysis_raw <- function(counts, design_mat, gr) {
    
    data.edgeR <- DGEList(counts = counts %>%
                            na.omit %>%
                            mutate(lipid = make.unique(lipid)) %>%
                            select( -type) %>%
                            # select(-Transition) %>%
                            column_to_rownames("lipid"),
                          group = gr
    )
    
    data.edgeR <- calcNormFactors(data.edgeR, method="TMM")
    data.edgeR <- estimateCommonDisp(data.edgeR, design=design_mat)
    data.edgeR
  }
  
  calculate_significance <- function(dge, contrast) {
    dge %>%
      glmFit() %>%
      glmLRT(contrast = contrast)
  }
  
  experiment_helper <- function(df) {
    df %>%
      perform_analysis_raw(design_expr, gr_expr) %>%
      calculate_significance(contrasts_expr) %>%
      topTags(1500000) %>%
      as.data.frame() %>%
      rownames_to_column("lipid") %>%
      as_tibble()
    
  }
  
  cl_e1_tbl <-
    cells_lipid_expr %>%
    experiment_helper
  cl_e1_tbl
  
  names(cells_lipid_expr)
  
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
  
  xx <- merge(cells_lipid_expr,cl_e1_tbl)
  names(cl_e1_tbl)
  names(xx)
  cl_e1_tbl
  # lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TG")
  # lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a")
  # 
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TG",'DAG','TG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  
  xx %>%
    
    ggplot(aes(x=logFC, y = type, fill = type)) +
    geom_density_ridges2(alpha = 0.5, size = .5) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    ggtitle(title_for_plot)+
    xlab("LogFC") +
    ylab("") +#xlim(-2, 4)+
    scale_alpha(guide = 'none') +
    # scale_fill_discrete(name = "Lipid class", guide = 'none') +
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
  
  
  cells_lipid_expr
  
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
}




