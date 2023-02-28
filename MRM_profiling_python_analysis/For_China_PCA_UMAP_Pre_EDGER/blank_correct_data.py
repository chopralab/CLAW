import pandas as pd
import matplotlib.pyplot as plt
import os



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
        lipid_names = list(df['lipid'])
        PVALUE = list(df['PValue'])
        logFC = list(df['logFC'])
        FDR = list(df['FDR'])

        lipid_types_old = list(df['type'])
        transition = list(df['Transition'])

        df = df.drop(['lipid', 'Transition','lipid','PValue','logFC','FDR','type',"logCPM","LR"], axis=1)

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

                if dbnumber == str(0):
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
        df.iloc[:, 1:-1] = df.iloc[:, 1:-1].subtract(df['SolventBlank1'], axis=0)
        # Make all negative numbers 0
        df[df < 0] = 0


        df = df.drop(['SolventBlank1'], axis=1)


        print(len(DB_number))
        print(len(db_present))
        print(len(sub_classes))


        # Add the 'type' column back to the dataframe and other columns
        df['type'] = type_col
        df['DB_number'] = DB_number
        df['DB_Present'] = db_present
        df['sub_class'] = sub_classes
        df['lipid'] = lipid_names
        df['PValue'] = PVALUE
        df['logFC'] = logFC
        df['FDR'] = FDR
        df['Transition'] = transition

        df.to_csv("Corrected_for_china/"+ij,index=False)