import pandas as pd
import matplotlib.pyplot as plt
import os

import numpy as np

lst=os.listdir('../data_files_name_changed')
print("df")

for ij in lst:
    if ".csv" in ij:
        # try:
        print("Hi")
        pi_title = ij[:-4]
            ##Gets path to open dataframe
        full_file_path = '../data_files_name_changed/'+ij


        df = pd.read_csv(full_file_path)

        old_list = list(df)
        print(old_list)
        # exit()
        limiting_length = len(df[df['FDR'] <0.1])

        lipid_names = list(df['lipid'])[:limiting_length]
        PVALUE = list(df['PValue'])[:limiting_length]
        logFC = list(df['logFC'])[:limiting_length]
        FDR = list(df['FDR'])[:limiting_length]

        lipid_types_old = list(df['type'])[:limiting_length]
        transition = list(df['Transition'])[:limiting_length]

        df = df.drop(['lipid', 'Transition','lipid','PValue','logFC','type',"logCPM","LR"], axis=1)

        df = df[df['FDR'] <0.1]

        df = df.drop(['FDR'], axis=1)


        ###Switches for proper amount of class
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
      # # Store the 'type' column separately

        type_col = new_types

        sub_classes = []
        chain_length = []
        DB_number = []
        db_present = []

        for i in range(len(lipid_types_old)):

            
            xx =(lipid_names[i])

            if ":" in xx:
                position = xx.find(":")

                chainlength = xx[position-2:position]
                dbnumber = xx[position+1:position+2]
                
                chain_length.append(chainlength)
                DB_number.append(dbnumber)

                if dbnumber == 0:
                    db_present.append("No")
                else:
                    db_present.append("Yes")

            else:

                chain_length.append("NA")
                DB_number.append("NA")
                db_present.append("NA")
        # print(abs(x))
        # if abs(x) > abs(threshold_value): ##For LogFC
            if (lipid_types_old[i]) == 'PCandSM':
                # xx =(lipid_names[i])
                print(xx)
                if 'PC' in xx:
                    if "Lyso" in xx:
                        
                        sub_classes.append("Lyso_PC")
    
                    elif "PCp" in xx:
                        
                        sub_classes.append("PCp")


                    elif "PCo" in xx:
                        
                        sub_classes.append("PCo")

                    elif "PC (" in xx:
                        
                        sub_classes.append("PC")

                else:
                    sub_classes.append("SM")



            elif (lipid_types_old[i]) == 'PI':



                if "Lyso" in xx:
                    
                    sub_classes.append("Lyso_PI")

                elif "PIp" in xx:
                    
                    sub_classes.append("PIp")

                elif "PIo" in xx:
                    
                    sub_classes.append("PIo")

                elif "PI (" in xx:
                    
                    sub_classes.append("PI")




            elif (lipid_types_old[i]) == 'PE':



                if "Lyso" in xx:
                    
                    sub_classes.append("Lyso_PE")

                elif "PEp" in xx:
                    
                    sub_classes.append("PEp")

                elif "PEo" in xx:
                    
                    sub_classes.append("PEo")

                elif "PE (" in xx:
                    
                    sub_classes.append("PE")



            elif (lipid_types_old[i]) == 'PG':



                if "Lyso" in xx:
                    
                    sub_classes.append("Lyso_PG")

                elif "PGp" in xx:
                    
                    sub_classes.append("PGp")

                elif "PGo" in xx:
                    
                    sub_classes.append("PGo")


                elif "PG (" in xx:
                    
                    sub_classes.append("PG")


            elif "TAG" in (lipid_types_old[i]):
                sub_classes.append("TAG")



            else:
                sub_classes.append(lipid_types_old[i])






        ##Drop Type
        # df = df.drop(['type'], axis=1)



        # Subtract all value columns (except 'type') by the 'SolventBlank1' column

        # Subtract all value columns (except 'type') by the 'SolventBlank1' column
        df.iloc[:, 0:-1] = df.iloc[:, 0:-1].subtract(df['SolventBlank1'], axis=0)
        # Make all negative numbers 0
        df[df < 0] = 0


        df = df.drop(['SolventBlank1'], axis=1)

        my_list = list(df)

        mid_length = int(len(my_list)/2)
        print(mid_length)
        print(my_list[2:])
        print(my_list[5:])
        first_half = my_list[mid_length:]
        second_half = my_list[:mid_length]
        print(my_list)
        print(first_half)
        print(second_half)

        first_name = first_half[0][:-1]+"_sum"
        second_name = second_half[0][:-1]+"_sum"

        print(first_name)
        print(second_name)


        

        # Add the 'type' column back to the dataframe
        df['type'] = sub_classes

        # Sum GFP columns together
        df[first_name] = df[first_half].sum(axis=1)
        df[second_name] = df[second_half].sum(axis=1)
        
        # Group by 'type' and sum up the GFP values
        grouped_df_1 = df.groupby('type')[first_name].sum()
        grouped_df_2 = df.groupby('type')[second_name].sum()


        print()

        print("")
        print("")
        print(list(grouped_df_2.index))
        print("")
        print("")
        print("")
        print("")
        print(grouped_df_1)
        print(grouped_df_2)
        print(grouped_df_1)
        print(list(grouped_df_1))
        print(type(grouped_df_1))
        # exit()

        list_values_1 = list(grouped_df_1)
        list_values_2 =list(grouped_df_2)
        lipid_names = list(grouped_df_1.index)
        barWidth = 0.25
        fig = plt.subplots(figsize =(12, 8))

        
        # Set position of bar on X axis
        br1 = np.arange(len(lipid_names))
        br2 = [x + barWidth for x in br1]
        # br3 = [x + barWidth for x in br2]

        # Make the plot





        # intensities_1
        # intensities_2




        plt.bar(br1, list_values_1, color ='#4C4E52', width = barWidth,
                edgecolor ='grey', label =first_name)
        plt.bar(br2, list_values_2, color ='gray', width = barWidth,
                edgecolor ='grey', label =second_name)
        # plt.bar(br3, counts_down, color ='b', width = barWidth,
        #         edgecolor ='grey', label ='Down Regulated')

        # Adding Xticks

        bar_plot = ij[:-4]

        plt.xlabel('Lipid Class', fontweight ='bold', fontsize = 15)
        plt.ylabel('Summed Intensity', fontweight ='bold', fontsize = 15)
        plt.xticks([r + barWidth for r in range(len(lipid_names))],lipid_names,rotation=90)
        plt.title(bar_plot.replace("_"," "))
        plt.tight_layout()


        plt.legend()
        plt.show()
        save_title_bar_up = "bar_plots_sub_plots/"+bar_plot +str("_FDR_bar.png")

        plt.savefig(save_title_bar_up ,dpi=700)
        # exit()

