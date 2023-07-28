###Need to add 
##LogFC of Blood and Plasma
####Relative Int

x <- 5

UMAP_generator_LOG_FC_1_FDR <- function(dataframe,lipid_list,color_vector=c(0,0.5,1),min_max=FALSE,limit_minimum=-1,
                                           limit_maximum=1,legend_title="LogFC 1 to -1",plot_title="Title",save_title="save_temp.png")
{
  # perm_plot =ggplot()
  counter = 0
  print(counter)
  # perm_plot =perm_plot+ somesetup
  for (i in lipid_list) {
    
    print(i)
    temp_df <-dataframe %>%
      filter(lipidtypes == i)
    temp_df2 <-dataframe %>%
      filter(lipidtypes == i,FDR<0.1)
    
    temp_plot<-ggplot(data =dataframe ,aes(x=X,y=Y)) + 
      geom_point(color="gray" ,size = 1)+geom_point(data = temp_df, aes(x=X, y=Y),color="black")+
      geom_point(data = temp_df2, aes(x=X, y=Y,color=LogFC), inherit.aes = FALSE)+
      ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))+xlab("")+ylab("")+
      #theme(legend.position="none")+labs(fill='LogFC')+
      scale_color_gradientn(legend_title,colors = c("blue","#DC9313","red"),limits = c(limit_minimum, limit_maximum),values = color_vector,oob = scales::squish)+
      theme(plot.title = element_text(hjust = 0.5))+ 
      theme(axis.text = element_text(size = axis_font_size))+ theme(axis.title = element_text(size = Axis_title_font))+ 
      theme(plot.title = element_text(size = title_font))          + theme(axis.text = element_text(size = axis_font_size))+ 
      theme(axis.title = element_text(size = Axis_title_font))+ theme(plot.title = element_text(size = title_font))+
      theme(panel.background = element_rect(fill = plot_background_color , color = 'black'), 
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = plot_grid_color ),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour =plot_grid_color )
      )
    
    # perm_plot=perm_plot+temp_plot
    
    if (counter==0){perm_plot=temp_plot}
    
    if (counter>0){perm_plot=perm_plot+temp_plot}
    counter= counter +1 
    print(counter)
    
    
    
    
  }
  # 
  # perm_plot+plot_layout(guides = "collect")+
  #   plot_annotation(title = plot_title,caption=UMAP_Hyperparameters, theme = theme(plot.title = element_text(size = major_title_font))) & 
  #   theme(plot.caption = element_text(size = caption_size),plot.title = element_text(hjust = 0.5))  
  # ##Adjust title to center
  
    perm_plot+plot_layout(guides = "collect")+ plot_annotation(title = "plot_title",caption=UMAP_Hyperparameters, theme = theme(plot.title = element_text(size = major_title_font))) & 
    theme(plot.caption = element_text(size = caption_size),plot.title = element_text(hjust = 0.5))  
  
  
  
  
  
  
  print(perm_plot) 
  ggsave(save_title, width=14, height=9)
  
  } 


















# rm(list=ls())
dot_size <-5
axis_font_size <-10
Axis_title_font <-20
title_font <- 30

major_title_font <-20
plot_background_color <- "white"
plot_grid_color <-"#5b5b5b"
letter_size <- 15
caption_size <- 20

####Plots MONA and DESI and COmparisons


library(plyr)
library(tidyverse)

library(readxl)


library(readxl)

library(ggplot2)

library(svglite)
library(stringi)
library(patchwork)
 
getwd()



setwd("C:/Users/Connor Beveridge/Box/connor_beveridge/New/ruilin_TBI_2_12_23/")





file_list = list.files(path="UMAP", pattern=NULL, all.files=FALSE,
           full.names=FALSE)

# file_list = c("1_vs_2_full_Umap_lipids_profiles_transpose_neighbors_30_minimum_distance_0.25_ruilin_TBI_11_21_22.xlsx")

file_list

file_list[1]

neighbors_ <- (30)
distances_ <- (0.25)
LOGFC_limit <-1
neighbors <-30
distance <- 0.25


lipid_list <- c("AC","FFA","CE","PG","CER","PC","PE","PI","PS","PG","SM","TAG")





for (neighbors in neighbors_)
{
  for (distance in distances_)
  {
    for (jj in file_list) {

      UMAP_Hyperparameters = "Neighbors = 30 Distance = 0.25"
      full_name = stri_sub(jj, 1, nchar(jj)-5)      
      
      begin_name = stri_sub(jj, 1, 10)   
      
      save_name =(paste("UMAP_plots/",begin_name,"UMAP LOGFC1 Starry Eyed.png",sep=""))
      
      Current_PCA_file11 <- read_excel(paste("UMAP/",jj,sep=''))
      PCA_title11<-(paste("Spine ",neighbors," Min Distance ",distance,sep=""))
      Umap1 <-"UMAP 1"
      Umap2 <-"UMAP 2"
    # save_title11 <-paste("UMAP_plots/Umap_lipids_profiles_transpose_neighbors_",neighbors,"_minimum_distance_",distance,"_Patient_15",sep='')
    
      UMAP1_11<-(Current_PCA_file11$`UMAP1`)
      UMAP2_11<-(Current_PCA_file11$`UMAP2`)
      
      # Blank_11<-(Current_PCA_file11$`SolventBlank`) 
      
      types_of_lipids11<-(Current_PCA_file11$`lipid Types`) 
      LOGFC_11<-(Current_PCA_file11$`LogFC`) 
      FDR_11<-(Current_PCA_file11$`FDR`) 
      
     

  
    
      UMAP_transpose_11 <- data.frame(lipidtypes=types_of_lipids11,X=UMAP1_11,Y=UMAP2_11,LogFC=LOGFC_11,FDR=FDR_11)
      
      UMAP_generator_LOG_FC_1_FDR(UMAP_transpose_11,lipid_list,save_title = save_name)
    
    }                           
  }
}
    
    # neighbors_ <- (30)
    # distances_ <- (0.25)
    # for (neighbors in neighbors_)
    # {
    #   for (distance in distances_)
    #   {
#######MAKING VARIABLE PLOT FUNCTION
    ###VARIABLES TO INSERT - dataframe with all dataneeded
    ###List of Lipids that one would like to use for the stary eyed name
    ###Must name lipidtypes for the lipidtypes
    #####Variables that define what X and Y and 
    ###what UMAP is being colored by
    ###Scaling of color bar have default values
    #####Color of color bar
    ###Filter based on lipid name
    ##SAVE Title
    ###Preset all font information
    # limit_minimum <- -1
    # limit_maximum <-1
    # color_vector<- c(0,0.5,1)
    # legend_title<- "LogFC 1 to -1"
    ##COLOR by is chagned in the function
    ##Currently you need to change clot_vs_brain_logfc
    ###To fix read files with just LogFC and FDR and have a list of file names that all get rid in and the app uses them
    

    
    
    
    
    
    
    
    
    
    # 
    # 
    # 
    # 
    # UMAP_generator_LOG_FC_1_FDR <- function(dataframe,lipid_list,color_vector=c(0,0.5,1),min_max=FALSE,limit_minimum=-1,
    #                                            limit_maximum=1,legend_title="LogFC 1 to -1",plot_title="title",save_title="save_temp.png")
    # {
    #   # perm_plot =ggplot()
    #   counter = 0
    #   print(counter)
    #   # perm_plot =perm_plot+ somesetup
    #   for (i in lipid_list) {
    #     
    #     print(i)
    #     temp_df <-dataframe %>%
    #       filter(lipidtypes == i)
    #     temp_df2 <-dataframe %>%
    #       filter(lipidtypes == i,FDR<0.1)   
    #     
    #     temp_plot<-ggplot(data =dataframe ,aes(x=X,y=Y)) + 
    #       geom_point(color="gray" ,size = 1)+geom_point(data = temp_df, aes(x=X, y=Y),color="black")+
    #       geom_point(data = temp_df2, aes(x=X, y=Y,color=clot_vs_brain_logfc), inherit.aes = FALSE)+
    #       ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))+xlab("")+ylab("")+
    #       #theme(legend.position="none")+labs(fill='LogFC')+
    #       scale_color_gradientn(legend_title,colors = c("blue","#DC9313","red"),limits = c(limit_minimum, limit_maximum),values = color_vector,oob = scales::squish)+
    #       theme(plot.title = element_text(hjust = 0.5))+ 
    #       theme(axis.text = element_text(size = axis_font_size))+ theme(axis.title = element_text(size = Axis_title_font))+ 
    #       theme(plot.title = element_text(size = title_font))          + theme(axis.text = element_text(size = axis_font_size))+ 
    #       theme(axis.title = element_text(size = Axis_title_font))+ theme(plot.title = element_text(size = title_font))+
    #       theme(panel.background = element_rect(fill = plot_background_color , color = 'black'), 
    #             panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = plot_grid_color ),
    #             panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
    #                                             colour =plot_grid_color )
    #       )
    #     
    #     # perm_plot=perm_plot+temp_plot
    #     
    #     if (counter==0){perm_plot=temp_plot}
    #     
    #     if (counter>0){perm_plot=perm_plot+temp_plot}
    #     counter= counter +1 
    #     print(counter)
    #     
    #     
    #     
    #     
    #   }
    #   
    #   perm_plot+plot_layout(guides = "collect")+
    #     plot_annotation(title = plot_title,caption=UMAP_Hyperparameters, theme = theme(plot.title = element_text(size = major_title_font))) & 
    #     theme(plot.caption = element_text(size = caption_size),plot.title = element_text(hjust = 0.5))  
    #   ##Adjust title to center
    #   
    #   print(perm_plot) 
    #   ggsave(save_title, width=14, height=9)
    #   
    # }
    # 
    # 
    # 
    # UMAP_generator_LOG_FC_1_FDR(UMAP_transpose_11,lipid_list)
    # 
    
    
    
    
    
    
    
    
    