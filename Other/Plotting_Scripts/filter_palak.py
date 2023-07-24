
import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt

# logFC
###Chooses FDR or P value for making bar plot
threshold = 'FDR'
# threshold = 'logFC'
# threshold_value = 1

# threshold = 'PValue'
threshold_value = 0.1
ncomponets = 4

##Create Directory for 
lst=os.listdir('data_files_name_changed')

lipid_types = ["CE","TAG","CER","FFA","PC","PE","PG","PI","SM","AC"]

lipid_types.sort()

for ij in lst:
    if ".csv" in ij:
        ##Directory and file name for the count of up and down regulated lipids
        new_xcel_file_name_for_palak = "FDR_filtered_xcel_sheets/"+ij[:-4]+'_'+threshold+'_'+str(threshold_value)+'filtered_for_palak_.csv'

        print(ij[:-4])
        # exit()
        ##Gets path to open dataframe
        full_file_path = 'data_files_name_changed/'+ij

        ###Opens the file and makes it a dataframe and adds the 
        df_main = pd.read_csv(full_file_path)

        df_main = df_main.sort_values('FDR')
        df_main = df_main[df_main['FDR'] <threshold_value]

        df_main.to_csv(new_xcel_file_name_for_palak)

        for iiii in range(len(lipid_types)):

            print(lipid_types[iiii])
            if lipid_types[iiii] =="TAG":

                df = df_main[(df_main['type'] == 'TAG2') | (df_main['type'] =="TAG1")]

            else:
                df = df_main[(df_main['type'] == lipid_types[iiii])]


            new_xcel_file_name_for_palak2 = "FDR_filtered_xcel_sheets_by_lipid_class/"+ij[:-4]+'_'+threshold+'_'+str(threshold_value)+'_LIPID_'+lipid_types[iiii]+'_filtered_for_palak_.csv'
            df.to_csv(new_xcel_file_name_for_palak2)
            