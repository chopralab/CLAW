# from msilib.schema import Component
import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA

from matplotlib.colors import ListedColormap, LinearSegmentedColormap

threshold = 'FDR'
# threshold = 'PValue'

lst=os.listdir('heatmap_data_files')

print('fdaf')
# lst = []
import re


# Original names
names222 = [
    "DOD100-M2-5xFAD-cerebellum10x_01202023", "DOD100-M2-5xFAD-cortex10x_01202023",
    "DOD100-M2-5xFAD-diencephalon10x_01202023", "DOD100-M2-5xFAD-hippo10x_01202023",
    "DOD100-M3-5xFAD-cerebellum10x_01202023", "DOD100-M3-5xFAD-diencephalon10x_01202023",
    "DOD100-M3-5xFAD-hippo10x_01202023", "DOD99-M1-5xFAD-cerebellum10x_01202023",
    "DOD99-M1-5xFAD-cortex10x_01202023", "DOD99-M1-5xFAD-diencephalon10x_01202023",
    "DOD99-M1-5xFAD-hippo10x_01202023", "M1FAD173-5xFAD-Cerebellum10x_01182023",
    "M1FAD173-5xFAD-Cortex10x_01182023", "M1FAD173-5xFAD-Diencephalon10x_01182023",
    "M1FAD173-5xFAD-Hippo10x_01182023", "M2FAD173-5xFAD-Cerebellum10x_01192023",
    "M2FAD173-5xFAD-Cortex10x_01192023", "M2FAD173-5xFAD-Diencephalon10x_01192023",
    "M2FAD173-5xFAD-Hippo10x_01192023", "DOD100-M1-WT-cerebellum10x_01202023",
    "DOD100-M1-WT-cortex10x_01202023", "DOD100-M1-WT-diencephalon10x_01202023",
    "DOD100-M1-WT-hippo10x_01202023", "DOD99-M2-WT-cerebellum10x_01202023",
    "DOD99-M2-WT-cortex10x_01202023", "DOD99-M2-WT-diencephalon10x_01202023",
    "DOD99-M2-WT-hippo10x_01202023", "M5FAD173-WT-Cerebellum10x_01192023",
    "M5FAD173-WT-Cortex10x_01192023", "M5FAD173-WT-Diencephalon10x_01192023",
    "M5FAD173-WT-Hippo10x_01192023", "10xBlank_01202023"
]

# New hardcoded names
new_names = [
    "5xFAD_Male_2_Cerebellum", "5xFAD_Male_2_Cortex",
    "5xFAD_Male_2_Diencephalon", "5xFAD_Male_2_Hippocampus",
    "5xFAD_Male_3_Cerebellum", "5xFAD_Male_3_Diencephalon",
    "5xFAD_Male_3_Hippocampus", "5xFAD_Male_1_Cerebellum",
    "5xFAD_Male_1_Cortex", "5xFAD_Male_1_Diencephalon",
    "5xFAD_Male_1_Hippocampus", "5xFAD_Male_1_Cerebellum",
    "5xFAD_Male_1_Cortex", "5xFAD_Male_1_Diencephalon",
    "5xFAD_Male_1_Hippocampus", "5xFAD_Male_2_Cerebellum",
    "5xFAD_Male_2_Cortex", "5xFAD_Male_2_Diencephalon",
    "5xFAD_Male_2_Hippocampus", "WT_Male_1_Cerebellum",
    "WT_Male_1_Cortex", "WT_Male_1_Diencephalon",
    "WT_Male_1_Hippocampus", "WT_Male_2_Cerebellum",
    "WT_Male_2_Cortex", "WT_Male_2_Diencephalon",
    "WT_Male_2_Hippocampus", "WT_Male_5_Cerebellum",
    "WT_Male_5_Cortex", "WT_Male_5_Diencephalon",
    "WT_Male_5_Hippocampus", "10xBlank_01202023"
]

# Verify the names:
ij = lst[0]

full_file_path = 'heatmap_data_files/'+ij


try:
    df_big = pd.read_csv(full_file_path)
except:
    df_big = pd.read_excel(full_file_path)

# for ij in lst:


print(ij[:-5])
# exit()

print(df_big['type'].unique())
# Loop through each unique value in the 'type' column
for jjjj in df_big['type'].unique():
    # Filter the dataframe for the current unique type
    df = df_big[df_big['type'] == jjjj]
    
    # Here, you can use 'filtered_df' as per your requirements
    # print(filtered_df)



    # print(df["type"][0])
    # exit()

        
    name_type = " "+jjjj
    my_list = list(df)
    print(my_list)
    

    end_of_list = len(my_list)-1
    print(end_of_list)
    print(my_list[end_of_list])
    # for j,i in enumerate(my_list):
    #     if 'Blank' == i:
    #         end_of_list=j
    for i in range(3,end_of_list):
            print(my_list[i])

    # exit()


    # exit()
    indexes = []

    for i in range(3,end_of_list):
        indexes.append(i)


    print(my_list[end_of_list])
    print(my_list[3])
    # exit()
    ncomponets = 4

    # exit()/
    PCA_title = "PCA_"+str(ncomponets)+"_componets_DIvided_by_blank"+ij[:-4]


    heat_title = "heat_maps_by_class/"+ij[:-4]+name_type+".svg"
    heat_title2 = "heat_maps_by_class/"+ij[:-4]+"_Log10_heatmap_division.png"


    blank_int = np.array(df[my_list[end_of_list]])

    print(blank_int)
    # exit()
    # files_list = [my_list[3],my_list[9],my_list[10],my_list[11]]
    print(full_file_path)
    itense_values = []
    itense_values_divide = []
    names = []
    non_normalized_values = []

    lines_plot = []





    for i in indexes:
        print(i)
        value = np.array(df[my_list[i]]) -blank_int
        value_divide = np.array(df[my_list[i]]) /blank_int
        value2 = df[my_list[i]].tolist()
        names.append(my_list[i])
        value3 = df[my_list[i]].tolist()

        value2.sort(reverse=True)
        lines_plot.append(value3)
        itense_values.append(value)
        itense_values_divide.append(value_divide)
        non_normalized_values.append(value2)
    print()
    print(len(itense_values))
    print(len(names))
    itense_values = np.array(itense_values)
    itense_values[itense_values<0]=0

    # Creating plot

    print(names)
    blank_sort = blank_int.tolist()
    blank_sort.sort(reverse=True)
    print(non_normalized_values[0])
    ###Taking Last numbers for the days and patient# and making my colors and dot types





    x = np.arange(len(non_normalized_values[0])).tolist() 



    colors = [[0,"blue"],[0.5,"yellow"], [1,"red"]]
    colors2 = [(0,"black"),(.149999,"black"),(.15,"purple"),(.2,"purple"),(.2001,"blue"),(.3,"blue"),(0.3000000001,"red"),(.5,"yellow"),(1,"green")]



    cmap1 = LinearSegmentedColormap.from_list("mycmap", colors)
    cmap2 = LinearSegmentedColormap.from_list("mycmap", colors2)
    list_of_lipid_classes = []

    for i in df["type"]:
        if "TAG" in i:
            list_of_lipid_classes.append("TAGs")
        # elif "PCandSM" in i:
        #     list_of_lipid_classes.append("PC/SM")
        # elif "Ceramides" in i:
        #     list_of_lipid_classes.append("Cer")
        else:
            list_of_lipid_classes.append(i)



    # list_of_lipid_classes = df["type"]

    list_of_lipid_classes_fixed = []
    # plt.figure(figsize=(500, 500))
    # plt.figure()
    list_of_lipid_classes_fixed.append(list_of_lipid_classes[0])
    for i in range(1,len(list_of_lipid_classes)):
        print(i)


        if list_of_lipid_classes[i] == list_of_lipid_classes[i-1]:
            list_of_lipid_classes_fixed.append(" ")
        elif list_of_lipid_classes[i] != list_of_lipid_classes[i-1]:
            list_of_lipid_classes_fixed.append(list_of_lipid_classes[i] )

    print('')
    print(len(list_of_lipid_classes_fixed))
    print(len(list_of_lipid_classes))

    # exit()
    plt.figure()










###this part matters
    heatmap = plt.pcolor(np.log10(itense_values),cmap=cmap1,vmin=0,vmax=6)#,vmax=15,vmin=0)
    # plt.xlabel(name_type)
    # plt.ylabel("Sample")
    plt.yticks(ticks=np.arange(len(names)),labels=names)
    plt.colorbar(heatmap)
    plt.xticks([])  # Remove x-axis numbering
    plt.yticks([]) 
    plt.title(jjjj)
    plt.tight_layout()
    plt.show()

    plt.savefig(heat_title,dpi=500)
    plt.savefig(heat_title[:-4]+".svg",dpi=500)

    plt.clf()

    plt.close()




    heatmap = plt.pcolor((itense_values_divide),cmap=cmap2,vmin=0,vmax=10)#,vmax=15,vmin=0)
    plt.xlabel(name_type)
    plt.ylabel("Sample")
    plt.yticks(ticks=np.arange(len(names)),labels=names)
    plt.colorbar(heatmap)
    plt.title("Intesnity/Blank "+name_type)
    plt.tight_layout()
    plt.show()

    plt.savefig(heat_title[:-4]+"_division.png",dpi=500)
    # plt.savefig(heat_title[:-4]+"_division.png",dpi=500)

    plt.clf()

    plt.close()


# exit()


    # print(itense_values)

    # pca = PCA(n_components = ncomponets)
    # print(ij)
    # components = pca.fit_transform(itense_values)
    # print(components)

    # print(components[0])
    # print(components[0][0])


    # X = []
    # Y = []
    # Z = []
    # print(len(components))

    # for i in range(len(components)):
    #     X.append(components[i][0])
    #     Y.append(components[i][1])
    #     Z.append(components[i][2])



    # print(pca.explained_variance_ratio_)

    # plt.scatter(X,Y,color="blue")
    # # plt.scatter(X,Y,color='r',marker="d")
    # xlabel = "PC1_"+ str(100*pca.explained_variance_ratio_[0])[:2]+"%"
    # ylabel = "PC2_"+ str(100*pca.explained_variance_ratio_[1])[:2]+"%"
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # for i, label in enumerate(names):
    #     plt.annotate(label, (X[i], Y[i]))

    # plt.legend()


    # plt.title(PCA_title)
    # # show plot
    # plt.show()
    # save_title_PCA = "prety_PCA/"+PCA_title +str("_PCA.svg")
    # plt.savefig(save_title_PCA)

    # plt.close()
    # plt.cla()
    # plt.clf()


    for i in range(len(lines_plot)):
        
        plt.figure()
        
        plt.plot(x,lines_plot[i], color='red',label=names[i])
        plt.plot(x,blank_int.tolist(), color='blue',label="blank")
        # plt.xticks(ticks=np.arange(len(list_of_lipid_classes_fixed)),labels=list_of_lipid_classes_fixed,rotation=90)

        plt.title(names[i]+name_type+' UnSorted List of Intensities')
        plt.xlabel("Unsorted of Lipid Intensities")
        plt.ylabel("Intensity")
        plt.legend()
        plt.tight_layout()

        plt.savefig("Ion_Distribution_by_class/"+names[i]+name_type+"_Non_sorted_.png")

        plt.clf()

        plt.close()

    print(np.amin(itense_values))
