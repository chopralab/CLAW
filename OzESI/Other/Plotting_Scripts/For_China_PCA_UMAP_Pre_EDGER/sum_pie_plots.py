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




        # # Store the 'type' column separately
        # type_col = df['type']
        df = df.drop(['type'], axis=1)

        my_list = list(df)
        print(my_list)
        print(my_list[:-1])

        my_list = my_list[:-1]
        # exit()


        # Subtract all value columns (except 'type') by the 'SolventBlank1' column
        df.iloc[:, 0:-1] = df.iloc[:, 0:-1].subtract(df['SolventBlank1'], axis=0)
        # Make all negative numbers 0

        df = df.drop(['SolventBlank1'], axis=1)


        df[df < 0] = 0
        # Add the 'type' column back to the dataframe
        df['type'] = type_col

        # Sum GFP columns together
        df['sum'] = df[my_list].sum(axis=1)
        
        # Group by 'type' and sum up the GFP values
        grouped_df = df.groupby('type')['sum'].sum()


        grouped_df.to_csv("lipid_make_up_DF/"+ij)

        save_title_pi = "lipid_make_up_plots/"+pi_title +str("_PIE_.png")

        # Create a pie chart
        grouped_df.plot(kind='pie', y=grouped_df.index, legend=False, autopct='%1.1f%%', startangle=90)

        plt.title(pi_title.replace("_"," "))
        plt.axis('equal')
        plt.ylabel('')
        plt.savefig(save_title_pi, dpi=600)

        plt.show()
        plt.clf()
        plt.close()