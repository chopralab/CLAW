import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt


full_file_path = "Fly_Brain_Intensities.xlsx"
full_file_path2 = "Fly_Brain_Intensities_fake_ozone2.xlsx"

df = pd.read_excel(full_file_path)

df_ozone_on= pd.read_excel(full_file_path2)


##Get list of names for finding blank
my_list = list(df)
my_list_ozone = list(df_ozone_on)
print(my_list)



###Finds the blank in order to use it for normilization also is the la
for j,i in enumerate(my_list):
    if 'Blank' == i:
        end_of_list=j

print(end_of_list)

print(df["Transition"])
print(df["lipid"])
print(df["type"])


print()



classes = list(df["type"])
transitions = list(df["Transition"])
transitions_ozone = list(df_ozone_on["Transition"])
# transitions_ozone = list(df_ozone_on["Transition No Ozone"])
transitions_ozone_without_ozone = list(df_ozone_on["Transition No Ozone"])
transitions_ozone_DB_position = list(df_ozone_on["Transition"])
lipids_no_ozone = list(df["lipid"])
lipids_ozone = list(df_ozone_on["lipid"])


# transitions = list(df["lipid"])
# transitions_ozone = list(df_ozone_on["lipid"])


print(my_list)


# exit()
files_list = []



##Intensities start at position 8
for i in range(3,end_of_list+1):
    files_list.append(my_list[i])

print(files_list)

blank_int = np.array(df[my_list[end_of_list]])

print(blank_int)
# exit()
# files_list = [my_list[8],my_list[9],my_list[10],my_list[11]]

itense_values = []
print(files_list[:-1])
# exit()
test = []
file_names = []
# for i in files_list[:-1]:
# for i in range(len(files_list[:-1])):
#     if classes[i] != "Ozid":
#         file_names.append(files_list[i])
#         value = np.array(df[files_list[i]]) - blank_int
#         print(np.array(value))
#         value = np.array(value) # /max(value)
#         print(type(value))
#         # value = np.log(value)
#         itense_values.append(value)

master_list_no_ozone = {}
master_list_no_ozone_lipid_names = {}


lipid_threshold = 100000
print('')
print('')
print('')
for i in range(len(files_list[:-1])):

    # if classes[i] != "Ozid":
        file_names.append(files_list[i])
        value = np.array(df[files_list[i]]) - blank_int

        # print(np.array(value))
        
        value = np.array(value) # /max(value)
        value_list = value.tolist()
        temp = []
        temp2 = []
        for j in range(len(value_list)):
            if classes[j] == "TAG1" or classes[j] ==  "TAG2":
                print(classes[j], "class", type(classes[j]))
                # xxxxx = if :
                print(classes[j])
                print(classes[j] == "TAG1" or classes[j] == "TAG2")
                if value_list[j] > lipid_threshold:
                    temp.append(transitions[j])
                    temp2.append(lipids_no_ozone[j])
            master_list_no_ozone[files_list[i]] = temp
            master_list_no_ozone_lipid_names[files_list[i]] = temp2


        
        # value = np.log(value)
        
        itense_values.append(value)

print('')
print('')
print('')



print('done')
print(file_names)
print(files_list)
print(np.array(itense_values).shape)




for i in range(len(itense_values)):
    

    print(itense_values[i])

print(master_list_no_ozone[file_names[-1]])
print(master_list_no_ozone[file_names[0]])
print(len(master_list_no_ozone[file_names[0]]))
print(len(master_list_no_ozone[file_names[-1]]))

master_list_ozone = {}
master_list_ozone_lipid_names = {}

# ozone_list_df = pd.read_excel("MUFA_info.xlsx")
trans_no_ozone = list(df_ozone_on["Transition No Ozone"])

print(len(master_list_no_ozone))
print(len(master_list_no_ozone))

for i in range(len(file_names)):
    # print(master_list_no_ozone[file_names[i]][5])
    temp = []
    temp2 = []
    print(master_list_no_ozone[file_names[i]])
    # print(master_list_no_ozone[file_names[i]][0])
    for j,jj in enumerate(master_list_no_ozone[file_names[i]]):

    # for j in range(len(master_list_no_ozone[file_names[i]][j])):
        print(j,jj, file_names[i])
        for ij in range(len(trans_no_ozone)):
            if trans_no_ozone[ij] == jj:
                temp.append(jj)
                break

    
    master_list_ozone[file_names[i]] = temp



print(master_list_ozone)
                # print("yahoo it matched")
                # exit()








for j,i in enumerate(my_list_ozone):
    if 'Blank' == i:
        end_of_list_ozone=j

blank_int_ozone = np.array(df_ozone_on[my_list_ozone[end_of_list_ozone]])

files_list_ozone = []

for i in range(4,end_of_list_ozone+1):
    files_list_ozone.append(my_list_ozone[i])


file_names_ozone = []


master_list_ozone_parsed = {}

for i in range(len(files_list_ozone[:-1])):

    # if classes[i] != "Ozid":
        file_names_ozone.append(files_list_ozone[i])
        print(type(blank_int_ozone))
        print((blank_int_ozone))
        print(type(np.array(df_ozone_on[files_list_ozone[i]])))
        print((np.array(df_ozone_on[files_list_ozone[i]])))
        # print()
        value = np.array(df_ozone_on[files_list_ozone[i]]) - blank_int_ozone

        # print(np.array(value))
        
        value = np.array(value) # /max(value)
        value_list = value.tolist()
        temp = []
        for j in range(len(value_list)):

            if value_list[j] > lipid_threshold:


                temp.append(transitions_ozone[j])
                master_list_ozone_parsed[files_list_ozone[i]] = temp


##master_list_ozone is dic of sample and transitions that are present in the no ozone run that match the ozone run
##master_list_ozone_parsed is dic of sample and transitions that are present with ozone on these will be matched to determine what is needed: we need to have sample names match





##master_list_ozone is dic of sample and transitions that are present in the no ozone run that match the ozone run
print(master_list_ozone)
print('')
##master_list_ozone_parsed is dic of sample and transitions that are present with ozone on these will be matched to determine what is needed: we need to have sample names match
print(master_list_ozone_parsed)






# exit()

ozone_present = {}





newest_names = []

for i in range(len(file_names)):

    temp_present_no_ozone = master_list_ozone[file_names[i]]
    temp_present_with_ozone = master_list_ozone_parsed[file_names[i]]

    print(temp_present_no_ozone)
    print(temp_present_with_ozone)
    print('')
    print('')
    temp = []
    for j in temp_present_no_ozone:

        for jj in range(len(transitions_ozone_without_ozone)):
            if j == transitions_ozone_without_ozone[jj]:
                
                
                current_tran = transitions_ozone[jj]
                ### Will hit multiple times
                # temp_present_with_ozone
                
                # temp = []
                # for ij in range(len(transitions_ozone)):

                for iijj in temp_present_with_ozone:
                    if current_tran == iijj:

                        temp.append(current_tran)
                        break
        ozone_present[file_names[i]] = temp
        newest_names.append(file_names[i])


print(ozone_present)
print('')
print(transitions_ozone)
print('')
print('')
print('')
print('')
print('')
print('')
print('')
print('')

print(ozone_present)
print('')

print(master_list_ozone)
print('')
##master_list_ozone_parsed is dic of sample and transitions that are present with ozone on these will be matched to determine what is needed: we need to have sample names match
print(master_list_ozone_parsed)

    # print(master_list_no_ozone[file_names[i]])
    # # print(master_list_no_ozone[file_names[i]][0])
    # for j,jj in enumerate(master_list_no_ozone[file_names[i]]):

    # # for j in range(len(master_list_no_ozone[file_names[i]][j])):
    #     print(j,jj, file_names[i])
    #     for ij in range(len(trans_no_ozone)):
    #         if trans_no_ozone[ij] == jj:
    #             temp.append(jj)
    #             break





# classes = list(df["type"])
# transitions = list(df["Transition"])
# transitions_ozone = list(df_ozone_on["Transition"])
# # transitions_ozone = list(df_ozone_on["Transition No Ozone"])
# transitions_ozone_without_ozone = list(df_ozone_on["Transition No Ozone"])
# transitions_ozone_DB_position = list(df_ozone_on["Transition"])
# lipids_no_ozone = list(df["lipid"])
# lipids_ozone = list(df_ozone_on["lipid"])

###Replace File_names with newer_names[i]
for i in range(len(file_names)):
    try:
        temp_trans = list(ozone_present[file_names[i]])
    except:
        continue
    print(file_names[i])
    temp = []
    for j in range(len(transitions_ozone)):
        for ij in temp_trans:
            if transitions_ozone[j]==ij:
                temp.append(lipids_ozone[j])

    print(temp)









