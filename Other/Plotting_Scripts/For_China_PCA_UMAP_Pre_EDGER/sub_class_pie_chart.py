import pandas as pd
import matplotlib.pyplot as plt
import os

# Import the dataframe
df = pd.read_csv('lipid_make_up/40GFP.csv')

lst=os.listdir('lipid_make_up')


for ij in lst:
    if ".csv" in ij:
        # try:
        print("Hi")
        pi_title = ij[:-4]
            ##Gets path to open dataframe
        full_file_path = 'lipid_make_up/'+ij



        df = pd.read_csv(full_file_path)



        lipid_names = list(df['lipid'])

        lipid_types_old = list(df['type'])
        df = df.drop(['lipid', 'Transition'], axis=1)


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

        type_col = new_types



        sub_classes = []


        for i in range(len(lipid_types_old)):

            
            xx =(lipid_names[i])


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










        # # Store the 'type' column separately
        # type_col = df['type']
        df = df.drop(['type'], axis=1)

        my_list = list(df)
        print(my_list)
        print(my_list[:-1])

        my_list = my_list[:-1]
        # exit()


        # Subtract all value columns (except 'type') by the 'SolventBlank1' column
        df.iloc[:,0:-1] = df.iloc[:, 0:-1].subtract(df['SolventBlank1'], axis=0)
        # Make all negative numbers 0


        df = df.drop(['SolventBlank1'], axis=1)

        df[df < 0] = 0
        # Add the 'type' column back to the dataframe
        df['sub_class'] = sub_classes

        # Sum GFP columns together
        df['sum'] = df[my_list].sum(axis=1)
        
        # Group by 'type' and sum up the GFP values
        grouped_df = df.groupby('sub_class')['sum'].sum()

        grouped_df.to_csv("lipid_make_up_DF_sub_plot/"+ij)



        save_title_pi = "lipid_make_up_sub_plots/"+pi_title +str("subplot_PIE_.png")

        # Create a pie chart
        grouped_df.plot(kind='pie', y=grouped_df.index, legend=False, autopct='%1.1f%%', startangle=90)

        plt.title(pi_title.replace("_"," "))
        plt.axis('equal')
        plt.ylabel('')
        plt.savefig(save_title_pi, dpi=600)

        plt.show()
        plt.clf()
        plt.close()