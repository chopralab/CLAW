
import os, umap
import numpy as np
import pandas as pd
import datatable as dt
import seaborn as sns
from matplotlib import pyplot as plt



neighbors = 500

distance = 0.5



neighbors_ = [30]

distance_ = [0.25]

df = pd.read_excel('/storage/connor_beveridge/Projects/5xFADvsWT_palak_caitlin/WT_vs_5xFAD_LD.xlsx')

lst=os.listdir('lipid_make_up')

for ij in lst:

    pi_title = ij[:-4]
        ##Gets path to open dataframe
    full_file_path = 'lipid_make_up/'+ij

    ###Opens the file and makes it a dataframe and adds the 
    df = pd.read_csv(full_file_path)
    transitions = list(df['Transition'])

    lipid_types_old = list(df['type'])
    lipid_names = list(df['lipid'])


    new_types = []
    for i in range(len(lipid_types_old)):
        if "TAG" in lipid_types_old[i]:
            new_types.append("TAG")
        elif "PCandSM" in lipid_types_old[i]:
            if "PC" in lipid_names[i]:
                new_types.append("PC")
            else:
                new_types.append("SM")
        else:
            new_types.append(lipid_types_old[i])

    df["New Types"] = new_types

    print(lipid_types_old)

    print(df.groupby(['New Types']).sum().reset_index())

    summed_df = df.groupby(['New Types']).sum().reset_index()

    print(list(summed_df))

    my_list = list(summed_df)

    # exit()

    for j,i in enumerate(my_list):
        if 'Blank' in i:
            end_of_list=j

    files_list = []
    lipid_types = summed_df["New Types"]

    ##Intensities start at position 8
    for i in range(1,end_of_list):
        print(my_list[i])
        files_list.append(my_list[i])

    blank_int = np.array(summed_df[my_list[end_of_list]])


    intense_values = []


    for i in files_list[:-1]:
        value = np.array(summed_df[i]) - blank_int
        print(np.array(value))
        value = np.array(value)
        value[value<0]=0 # /max(value)
        print(type(value))
        # value = np.log(value)
        intense_values.append(value)

    intense_values = np.array(intense_values)

    print(np.array(intense_values).shape)
    print(summed_df)
    summed_intense_values = np.sum(intense_values,0)

    print(summed_intense_values.shape)
    print(lipid_types)
            
    lipid_list_for_pie = list(lipid_types)
    intensity_list_for_pie = list(summed_intense_values)


    print(summed_intense_values)





    ###new colors
    colors = []
    # colors.append(["ACs","#E49B0F"])
    colors.append(["AC","#E49B0F"])
    # colors.append(["CEs","#000000"])
    colors.append(["CE","#666699"])
    colors.append(["FFA","#4596a7"])
    # colors.append(["TAG2","#de2e95"])
    colors.append(["TAG","#cc40cc"])
    # colors.append(["TAGs","#de2e95"])
    # colors.append(["TAG1","#de2e95"])
    # colors.append(["TAG","#de2e95"])
    colors.append(["SM","red"])
    # colors.append(["Sm","red"])
    # colors.append(["sm","red"])
    colors.append(["PC","#6e34a4"])
    # colors.append(["PEs","#c3a2e2"])
    colors.append(["PE","#c3a2e2"])
    colors.append(["PG","#9fc5e8"])
    colors.append(["PGs","#9fc5e8"])
    colors.append(["PSs","#e1a7c3"])
    colors.append(["PS","#e1a7c3"])

    colors.append(["CER","#b6d7a8"])

    colors.append(["PI","#cabcab"])



    colors_2_use = []

    print("done")

    # lipid_classes.sort()

    for i in range(len(lipid_list_for_pie)):
        for j in range(len(colors)):
            if lipid_list_for_pie[i] == colors[j][0]:
                colors_2_use.append(colors[j][1])
                print(colors[j][0],"  ",lipid_list_for_pie[i],"    ",colors[j][1])
    print(lipid_list_for_pie)

    # exit()
    print(lipid_list_for_pie)
    print(colors_2_use)
    fig = plt.figure(figsize =(10, 7))
    plt.pie(intensity_list_for_pie, labels = lipid_list_for_pie,autopct='%.1f%%')#,colors=colors_2_use )
    plt.title(pi_title.replace("_"," "))
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    save_title_pi = "lipid_make_up/"+pi_title +str("_PIE_.png")


    # show plot
    plt.show()
    plt.tight_layout()

    plt.savefig(save_title_pi, dpi=600)

    plt.close()
    plt.cla()
    plt.clf()








