
library(plyr)
library(ggridges)
library(tidyverse)
library(readxl)
library(edgeR)
library(tidyverse)
# library(factoextra)
library(grid)
library(ggsignif)
library(stringr)
library(reshape2)
library(gtools)
library(tidyr)
library(ggbeeswarm)
library(ggrepel)
library(scales)

library(ggplot2)

library(cowplot)
library(ggridges)



# library(ggforce)
getwd()

# read the variable from the text file
cwd <- readLines("Variable_Storage/folder_path.txt")[1]
cwd
setwd(cwd)



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

  
  if (isTRUE(!is.na(length1)) && isTRUE(!is.na(length2))) {
    if (isTRUE(length1 == 1) || isTRUE(length2 == 1)) {
      next
    }
  }
  
  
  
  gr1 = c(rep("GR1",length1))
  gr2 = c(rep("GR2",length2))
  
  
  groups_PCA = c((rep(Title1,length1)),(rep(Title2,length2)))
  

  
  
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
  

  
  xx <- merge(cells_lipid_expr,cl_e1_tbl)
  names(cl_e1_tbl)
  names(xx)
  cl_e1_tbl
  # lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG")
  # lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a")
  # 
  lipid_classes <- c("CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE')
  lipid_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3')
  
  
  
  # Create a named vector to map lipid classes to their colors
  lipid_class_colors <- setNames(lipid_colors, lipid_classes)
  
  xx %>%
    
    ggplot(aes(x=logFC, y = type, fill = type)) +
    geom_density_ridges2(alpha = 0.5, size = .5) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    ggtitle(title_for_plot)+
    xlab("Fold change, lipids all") +
    ylab("") +#xlim(-2, 4)+
    scale_alpha(guide = 'none') +
    # scale_fill_discrete(name = "Lipid class", guide = 'none') +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
    scale_y_discrete(limits = rev) #+
  ggsave(paste("plots/Ridge_Plot_All Lipids_",title_for_plot,".png",sep=''))
  
  filtered_xxx <- xx %>%
    filter(FDR < 0.1)

  if (isTRUE(nrow(filtered_xxx) < 2)) {
    next
  }

  xx %>%

    filter(FDR < 0.10) %>% #filter by sig
    ggplot(aes(x=logFC, y = type, fill = type)) +
    geom_density_ridges2(alpha = 0.5, size = .5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    ggtitle(title_for_plot)+
    xlab("Fold change, lipids FDR<0.10") +
    ylab("") +#xlim(-2, 4)+
    scale_alpha(guide = 'none') +

    scale_fill_discrete(name = "Lipid class", guide = 'none') +
    scale_fill_manual(values = lipid_class_colors, name = "Lipid class", guide = 'none') +
    scale_y_discrete(limits = rev) #+
  # facet_wrap(~comparison_LFC)

  ggsave(paste("plots/Ridge_Plot_Filtered_FDR_",title_for_plot,".pdf",sep=''))
  
  make_volcano_plot <- function(df, title) {
    df %>%
      mutate(sig = factor(FDR < 0.10)) %>%
      ggplot(aes(logFC, -log10(FDR), color = sig)) +
      geom_point() +
      scale_color_manual(values = c("none" = "black", "TRUE" = "red")) +
      guides(color = F) +
      ggtitle(title)
  }

  cl_e1_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot,sep=''))
  ggsave(paste("plots/Volcano_Plot_",title_for_plot,".png",sep=''))

  
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
  
  



  
  
  
}




