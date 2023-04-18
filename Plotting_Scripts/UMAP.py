


import os, umap
import numpy as np
import pandas as pd
import datatable as dt
import seaborn as sns




neighbors = 500

distance = 0.5



neighbors_ = [30]

distance_ = [0.25]

df = pd.read_excel('/storage/connor_beveridge/Projects/fly_age_part_2/Fly_Brain_Intensities__part4_other_lipids.xlsx')

print(list(df))

print('yo')


print(list(df))
##Actual Transition
lipid_names = list(df['Transition'])


lst=os.listdir('data_files_name_changed')

for ij in lst:


    if ".csv" in ij:

        lipid_names = list(df['Transition'])

        full_file_path = 'data_files_name_changed/'+ij

        name_of_file = (ij[:-4])
        print(name_of_file)
        ###Opens the file and makes it a dataframe and adds the 
        temp_df = pd.read_csv(full_file_path)

        temp_df_lipids_unsorted = list(temp_df["Transition"])
        temp_df_lipids_LOGFC_unsorted = list(temp_df["logFC"])
        temp_df_lipids_FDR_unsorted = list(temp_df["FDR"])
        temp_df_LogFC = []
        temp_df_FDR = []
        temp_df_lipids = []




        for i in range(len(lipid_names)):
            for j in range(len(temp_df_lipids_unsorted)):
                if lipid_names[i] ==temp_df_lipids_unsorted[j]:
                    print(lipid_names[i] ,"     ",temp_df_lipids_unsorted[j])
                    temp_df_LogFC.append(temp_df_lipids_LOGFC_unsorted[j])
                    temp_df_lipids.append(temp_df_lipids_unsorted[j])
                    temp_df_FDR.append(temp_df_lipids_FDR_unsorted[j])
                    break
                continue


        my_list = list(df)
        print(my_list)

        # exit()
        print(df['type'])
        for j,i in enumerate(my_list):
            if 'Blank' in i:
                end_of_list=j



        # lipid_names = df['lipid']
        print(lipid_names)


        indexes = []

        for i in range(3,end_of_list):
            indexes.append(i)

        blank_int = np.array(df[my_list[end_of_list]])

        # print(blank_int)

        itense_values = []
        names = []

        for i in indexes:
            # print(i)
            value = np.array(df[my_list[i]])-blank_int
            names.append(my_list[i])

            value = np.array(value) #/max(value)

            # exit()
            itense_values.append(value)

        itense_values = np.array(itense_values)
        # print(min(itense_values.any()))

        itense_values[itense_values<0 ] = 0
        # print(min(itense_values.any()))
        # exit()
        # # Remove expression features with > 50% zero-valued expression levels     
        # is_expressed = np.apply_along_axis(lambda x: np.mean(x == 0) < .5, arr=matrix, axis=0)
        # matrix = matrix[:,is_expressed.tolist()]

        # # Log2-transform
        # matrix = np.log2(matrix.to_numpy() + 1)

        print(itense_values.shape)
        print(itense_values.transpose().shape)


        intesity_matrix = itense_values.transpose()

        
        lipids_types = list(df["type"])
        lipids_types2 = []

        for j,i in enumerate(lipids_types):
            if "TAG" in i or "Tag" in i:
                lipids_types2.append("TAG")
            elif "PCandSM" in i or "PCSM" in i:
                if "PC" in df['lipid'][j] or "pc" in df['lipid'][j]:
                    lipids_types2.append("PC")
                else: # "PC" in df['lipid'][j] or "pc" in df['lipid'][j]:
                    lipids_types2.append("SM")           
            else:
                lipids_types2.append(i)


                # Define UMAP





        for neighbors in neighbors_:

            for distance in distance_:


                brain_umap = umap.UMAP(random_state=999, n_neighbors=neighbors, min_dist=distance)

                # Fit UMAP and extract latent vars 1-2
                # embedding = pd.DataFrame(brain_umap.fit_transform(itense_values), columns = ['UMAP1','UMAP2'])

                embedding2 = pd.DataFrame(brain_umap.fit_transform(itense_values.transpose()), columns = ['UMAP1','UMAP2'])



                umap1 = []
                umap2 = []

                for i in range(len(embedding2['UMAP1'])):
                    umap1.append(embedding2['UMAP1'][i])
                    umap2.append(embedding2['UMAP2'][i])

                print(len(umap1))
                print(len(umap2))
                print(len(umap2))

                print(len(list(df['lipid'])),"list(df['lipid'])")
                print(len(lipid_names),"lipid_names")
                print(len(temp_df_lipids),"confirm_lipids")

                print(len(lipids_types2),"lipid Types")
                print(len(umap1),"umap1")
                print(len(umap2),"umap2")

                print(len(temp_df_LogFC),"temp_df_LogFC")
                print(len(temp_df_FDR),"temp_df_FDR")
                print(ij)

                df1 = {"Lipid Names":list(df['lipid']),"Transitions":lipid_names,"confirm_lipids":temp_df_lipids,"lipid Types":lipids_types2,"UMAP1":umap1,"UMAP2":umap2,"LogFC":temp_df_LogFC,"FDR":temp_df_FDR}



                df1 = pd.DataFrame(df1)
                print(df1)
                # exit()
                df1.to_excel("UMAP/"+ij[:-4]+"_Umap_lipids_profiles_transpose_neighbors_"+str(neighbors)+"_minimum_distance_"+str(distance)+"_fly_samples_December_22.xlsx",index=None)











