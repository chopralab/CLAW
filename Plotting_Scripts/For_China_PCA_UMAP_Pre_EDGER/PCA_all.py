import os, umap

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Load the data

# all = True

# if all == False:
# # Load the data from an Excel file
#     df = pd.read_csv('No_Neurons.csv', index_col=0)
#     folder_name = "No_Neurons/"
#     ij = "No Neurons"
# if all == True:
#     ij = "All"
#     folder_name = "All/"
#     df = pd.read_csv('all.csv', index_col=0)



lst=os.listdir('data_for_PCA')


for ij in lst:

    folder_name = "new_plots_PCA_FUll/"
    full_file_name = 'data_for_PCA/'+ij
    df = pd.read_csv(full_file_name, index_col=0)




    my_list = list(df)

    print(my_list)

    my_list = my_list[:-1]
    print(my_list)

    age = []
    sample_type = []
    cell_type = []


    print(my_list)
    # exit()
    for i in my_list:
        print(i, "Good")
        if i[0] =="5" or i[0]=="7":
            age.append("Young")
            print(i)
        elif i[0] =="4":
            print(i)
            age.append("Old")
        if i[1] =="E" or i[2] =="E":
            print(i)
            sample_type.append("AP1 blocked")
        if i[1] =="C"or i[2] =="C":
            print(i)
            sample_type.append("Control")



        if "GFP" in i:
            print(i)
            sample_type.append("Non-senescent glia")
        if "DP" in i:
            print(i)
            sample_type.append("Senescent glia")
        if "DN" in i:
            print(i)
            sample_type.append("Neuron")
    print(len(sample_type))
    print(len(age))

    # exit()


    # Get the column names for the data
    cols = df.columns.tolist()
    print(len(cols))
    # print(cols[32])
    # exit()


    # Get the blank column name
    blank_col = cols[-1]

    # Subtract the blank column from each data column
    df_subtracted = df.sub(df[blank_col], axis=0)

    del df_subtracted[blank_col]

    print("")
    print("")
    print("")
    print("")
    print(list(df_subtracted))
    # exit()
    # for type_name, type_df in df_subtracted.groupby('type'):
        # Get the data values for this type
    type_data = df_subtracted.values
    
    # Normalize the data by subtracting the mean and dividing by the standard deviation
    # type_data_norm = (type_data - np.mean(type_data, axis=0)) / np.std(type_data, axis=0)
    type_data_norm = type_data #/max(np.array(type_data).any()) 
    
    # Perform PCA
    pca = PCA()
    pca.fit(type_data_norm)
    
    # Get the explained variance of each principal component
    explained_var = pca.explained_variance_ratio_
    
    # Plot the results

    brain_umap = umap.UMAP(random_state=999, n_neighbors=30, min_dist=0.25)

    
    embedding2 = pd.DataFrame(brain_umap.fit_transform(np.array(type_data_norm).transpose()), columns = ['UMAP1','UMAP2'])



    umap1 = []
    umap2 = []

    for i in range(len(embedding2['UMAP1'])):
        umap1.append(embedding2['UMAP1'][i])
        umap2.append(embedding2['UMAP2'][i])

    fig, ax = plt.subplots()
    ax.scatter(pca.components_[0], pca.components_[1])
    ax.set_xlabel('PC1 Variance: {:.2f}%'.format(explained_var[0]*100))
    ax.set_ylabel('PC2 Variance: {:.2f}%'.format(explained_var[1]*100))
    ax.set_title(ij[:-4].replace("_"," "))
    print(len(type_data))
    print(len(cols))
    print(len(df))
    print(cols)
    for i in range(len(cols)-1):
        
        plt.annotate(cols[i], (pca.components_[0][i], pca.components_[1][i]))

    fig.savefig(folder_name+ij[:-4] + '_PCA.png')


    plt.clf()
    plt.close()

    fig, ax = plt.subplots()
    ax.scatter(umap1, umap2)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title(ij[:-4].replace("_"," "))
    print(len(type_data))
    print(len(cols))
    print(len(df))
    print(cols)
    print(len(umap1))
    print(len(umap2))

    temp1 = []
    temp2 = []
    temp3 = []
    for i in range(len(cols)-1):
        print(i,umap1[i],umap2[i],len(cols),cols[i])
        temp1.append('PC1 Variance: {:.2f}%'.format(explained_var[0]*100))
        temp2.append('PC2 Variance: {:.2f}%'.format(explained_var[1]*100))
        temp3.append(ij[:-4].replace("_"," "))
        
        plt.annotate(cols[i], (umap1[i], umap2[i]))

    fig.savefig(folder_name+ij[:-4].replace("_"," ") + '_UMAP.png')

    print(len(cols[:-1]))
    print(len(age))
    print(len(sample_type))
    print(len(temp1))
    print(len(umap1))

    df1 = {"Type":temp3,"Name":cols[:-1],"Age":age,"Sample_Type":sample_type,"PC1 Variance":temp1,"PC2 Variance":temp2,"UMAP1":umap1,"UMAP2":umap2,"PC1":pca.components_[0],"PC2": pca.components_[1]}
    df1 = pd.DataFrame(df1)

    df1.to_csv(folder_name+ij[:-4]+"_"+ij[:-4].replace("_"," ")+"_Dataframe_To_Plot.csv",index=None)

    # Show the plots
