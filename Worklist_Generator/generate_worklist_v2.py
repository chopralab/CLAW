

# Import pandas
import csv
import pandas as pd



data_frame_read = pd.read_csv('./Path_info.csv')

path_to_methods = data_frame_read["path_to_methods"][0]
path_to_save_data = data_frame_read["path_to_save_data"][0]
clean_file_method = data_frame_read["clean_file_method"][0]


clean_file_data_save = path_to_save_data

clean_every_n = int(input("Clean Every N: "))
n_cleans = int(input("How many cleans every N: "))
n_cleans_after_lipid_class =int(input("How many cleans every lipid class: "))
#clean_solution_position = "p1-a9"

clean_solution_position = input("Clean Solution Position EX p1-a9: ")


biological_replicates = int(input("How Many Technical Replicates?: "))


run_information = input("Run Information (Project Name) EX: Burda_Lab: ")
date_taken = input("Run Date in format 3_3_23: ")



method_dict = {}

with open('methods.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        method_dict[row['annotation']] = row['method name']

print(method_dict)





sample_types = []
lipid_types = []
sample_positions = []


data_frame_read = pd.read_csv('./sample_names.csv')

sample_types = data_frame_read["Sample Name"]
sample_positions = data_frame_read["Position"]



file2 = open('./lipid_classes.csv', 'r')
Lines2 = file2.readlines()

for line in Lines2:
    # print(line.strip())
    lipid_types.append(line.strip())



counter = 0

worklist_to_make = []


for q in range(1,(biological_replicates+1)):
    
    replicate_level = "_N"+str(q)+"_"+run_information+"_"+date_taken 
    for i in lipid_types:

        for index,j in enumerate(sample_types):
            sample_no_lipid = j+replicate_level


            sample_with_lipid = path_to_save_data+i+"_"+sample_no_lipid+".d"
            method_to_use =path_to_methods+ method_dict[i]
            counter = counter +1
            try:
                worklist_to_make.append([sample_positions[index].replace(",","_"),method_to_use.replace(",","_"),sample_with_lipid.replace(",","_")])
            except:
                worklist_to_make.append([sample_positions[index].replace(",","_"),method_to_use.replace(",","_"),sample_with_lipid.replace(",","_")])

            if counter%clean_every_n == 0:
                if n_cleans>1:
                    for ii in range(n_cleans):
                        clean_file = clean_file_data_save+str(counter)+str(ii)+".d"
                        worklist_to_make.append([clean_solution_position,clean_file_method ,clean_file ])


                else:

                    clean_file = clean_file_data_save+str(counter)+".d"
                    worklist_to_make.append([clean_solution_position,clean_file_method ,clean_file ])

        if n_cleans_after_lipid_class>1:
            for ii in range(n_cleans):
                clean_file = clean_file_data_save+str(counter)+str(ii)+"post_lipid_class_cleaning"+".d"
                worklist_to_make.append([clean_solution_position,clean_file_method ,clean_file ])
            else:
                clean_file = clean_file_data_save+str(counter)+"post_lipid_class_cleaning"+".d"

                worklist_to_make.append([clean_solution_position,clean_file_method ,clean_file ])




#print(worklist_to_make)
import csv
file_name_2_save = input("Save File Name(No file extension)?")

with open('./worklists/'+file_name_2_save+'.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    for i in worklist_to_make:
        spamwriter.writerow(i)


