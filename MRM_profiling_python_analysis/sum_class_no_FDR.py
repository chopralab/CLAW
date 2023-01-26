

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
        new_xcel_file_name_for_palak = "sig_lipids_count/"+ij[:-4]+'_'+threshold+'_'+str(threshold_value)+'_.xlsx'

        print(ij[:-4])
        # exit()
        ##Gets path to open dataframe
        full_file_path = 'data_files_name_changed/'+ij

        ###Opens the file and makes it a dataframe and adds the 
        df_main = pd.read_csv(full_file_path)

        df_main = df_main.sort_values('FDR')
        # df_main = df_main[df_main['FDR'] <threshold_value]



        intensities1 = []
        intensities2 = []
        lipid_names = []

        stds1 = []
        stds2=[]

        for iiii in range(len(lipid_types)):

            print(lipid_types[iiii])
            if lipid_types[iiii] =="TAG":

                df = df_main[(df_main['type'] == 'TAG2') | (df_main['type'] =="TAG1")]

            else:
                df = df_main[(df_main['type'] == lipid_types[iiii])]



            if len(df)<1:
                continue

            lipid_names.append(lipid_types[iiii])


            

            Log_FC = list(df["logFC"])[:6]



            print(lipid_names)
            print(Log_FC)





        # exit()

        ##Get list of names for finding blank
            my_list = list(df)
            print(my_list)



            ###Finds the blank in order to use it for normilization also is the la
            for j,i in enumerate(my_list):
                if 'Blank' in i:
                    end_of_list=j
            
            print(end_of_list)



            # exit()
            files_list = []



            ##Intensities start at position 8
            for i in range(8,end_of_list+1):
                files_list.append(my_list[i])






            # exit()

            # print(my_list[8])
            # print(my_list[9])
            # print(my_list[10])
            # print(my_list[11])
            # print(my_list[12])
            # # exit()
            bar_plot = ij[:-4]+' _Bar_Plot_summed_relative_intensity_No_Filter'


            blank_int = np.array(df[my_list[end_of_list]])

            print(blank_int)
            # exit()
            # files_list = [my_list[8],my_list[9],my_list[10],my_list[11]]

            intense_values = []
            print(files_list[:-1])
            # exit()
            for i in files_list[:-1]:
                value = np.array(df[i]) - blank_int
                print(np.array(value))
                value = np.array(value) # /max(value)
                print(type(value))
                # value = np.log(value)
                intense_values.append(sum(value))


            intense_values = np.array(intense_values)

            print(intense_values.shape)

            dividing_line = int(intense_values.shape[0]/2)
            print(dividing_line)
            print(intense_values[:dividing_line].shape)
            print(intense_values[dividing_line:].shape)

            first_half = intense_values[:dividing_line]
            second_half = intense_values[dividing_line:]

            mean1 = np.mean(first_half)
            std1 = np.std(first_half)
            mean2 = np.mean(second_half)
            std2 = np.std(second_half)

            print(mean1.shape)
            print(std1.shape)
            print(mean2.shape)
            print(std2.shape)
            print(mean1)
            print(std1)
            print(mean2)
            print(std2)

            max_list = [mean1,mean2]

            mean1=mean1/max(max_list)
            std1=std1/max(max_list)
            mean2=mean2/max(max_list)
            std2=std2/max(max_list)


            intensities1.append(mean1)
            stds1.append(std1)
            intensities2.append(mean2)
            stds2.append(std2)



            # exit()


        barWidth = 0.25
        fig = plt.subplots(figsize =(12, 8))







        # Set position of bar on X axis
        br1 = np.arange(len(lipid_names))
        br2 = [x + barWidth for x in br1]
        # br3 = [x + barWidth for x in br2]

        # Make the plot

        print(len(intensities1))
        print(len(intensities2))
        print(len(lipid_names))
        print((lipid_names))
        print(len(stds2))
        print(len(intensities1))
        plt.bar(br1, intensities1, color ='#4C4E52', width = barWidth,
                edgecolor ='grey', label =files_list[0][:-1],yerr=stds1)
        plt.bar(br2, intensities2, color ='gray', width = barWidth,
                edgecolor ='grey', label =files_list[-2][:-1],yerr=stds2)
        # plt.bar(br3, counts_down, color ='b', width = barWidth,
        #         edgecolor ='grey', label ='Down Regulated')

        # Adding Xticks
        plt.xlabel('Lipid Class', fontweight ='bold', fontsize = 15)
        plt.ylabel('Average Summed Intensity Relative', fontweight ='bold', fontsize = 15)
        plt.xticks([r + barWidth for r in range(len(lipid_names))],lipid_names,rotation=90)
        plt.title(bar_plot.replace("_"," "))
        plt.tight_layout()


        plt.legend()
        plt.show()
        save_title_bar_up = "class_sum/"+bar_plot +str("bar.png")

        plt.savefig(save_title_bar_up ,dpi=700)


        plt.close()
        plt.cla()
        plt.clf()