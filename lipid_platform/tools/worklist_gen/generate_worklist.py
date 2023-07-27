
###Sample Names
##Dates it will be run
###number of biological replicates
###What lipid classes will be run
####Datafile save path
###Method file
###Need sample position

###Add clean_X samples
###Add clean after every lipid class
###Add wash and standby mode at the end



# Import pandas
import pandas as pd

# path_to_methods = 'D:\MassHunter\methods\Lipid MRMs\Lipids NEW - 07142022\\LDs-01312023\\'
#path_to_methods = 'D:\MassHunter\methods\Lipid MRMs\Lipids NEW - 07142022\\'
path_to_methods = "D:\\MassHunter\\methods\\Lipid MRMs\\Updated MRM Lipids 03-31-2023\\"
path_project_name = input('Folder for Project Name (You must create it on LCMS computer) (No Slashes):  ')
path_project_subfolder = input("Sub Folder which is located in project name, such as the run data (No Slashes)EX: 3_4_23:  ")

#path_to_save_data = "D:\MassHunter\Data\caitlin_randolph\Lipid MRMs\Amita Lab\\03-24-2023\\"
#path_to_save_data = "D:\MassHunter\Data\caitlin_randolph\Lipid MRMs\\"+path_project_name+"\\"+path_project_subfolder+"\\"
path_to_save_data = "D:\MassHunter\Data\caitlin_randolph\OzESI\\"+path_project_name+"\\"+path_project_subfolder+"\\"




# print("D:\MassHunter\Data\caitlin_randolph\OzESI\\01-31-2023\\")
# exit()
#clean_file_method = "D:\MassHunter\methods\Lipid MRMs\Lipids NEW - 07142022\\Tailored MRM-profile_TAG2_04-18-2019_IS.m"
clean_file_method ="D:\\MassHunter\\methods\\Lipid MRMs\\Updated MRM Lipids 03-31-2023\\MRM profiling_TAG14_0_05-10-2022.m"

clean_file_data_save = path_to_save_data

clean_every_n = int(input("Clean Every N: "))
n_cleans = int(input("How many cleans every N: "))
n_cleans_after_lipid_class =int(input("How many cleans every lipid class: "))
#clean_solution_position = "p1-a9"

clean_solution_position = input("Clean Solution Position EX p1-a9: ")


biological_replicates = int(input("How Many Technical Replicates?: "))


# run_information = "_Burda_"
# date_taken = "3_17_23"


run_information = input("Run Information (Project Name) EX: Burda_Lab: ")
date_taken = input("Run Date in format 3_3_23: ")

method_dict = {}



# method_dict["FFA"] = "Tailored MRM-profile_FFA_02-25-2018_ver6.m"
# method_dict["TAG1"] = "Tailored MRM-profile_TAG1_04-18-2019_IS.m"
# method_dict["TAG2"] = "Tailored MRM-profile_TAG2_04-18-2019_IS.m"
# method_dict["PG"] = "Tailored MRM-profile_Lyso_PG_GP04_12-05-2019_IS_ver6.m"
# method_dict["PCandSM"] = "Tailored MRM-profile_Lyso_PC and SM_12-05-2019_IS.m"
# method_dict["PE"] = "Tailored MRM-profile_Lyso_PE_GP02_12-05-2019_ver6.m"
# method_dict["PI"] = "Tailored MRM-profile_Lyso_PI_GP06_12-05-2019_IS.m"
# method_dict["PS"] = "Tailored MRM-profile_Lyso_PS_GP03_12-05-2019_IS.m"
# method_dict["cardio"] = "Tailored MRM-profiling_Cardiolipin__NEG_01-06-2021.m"

# method_dict["DAG"] = "Tailored MRM-profile_DAGs.m"
# method_dict["CE"] = "Tailored MRM-profile_cholesteryl esters_04-18-2019_IS.m"
# method_dict["CER"] = "Tailored MRM-profile_Ceramides_IS_05-03-2019_ver6.m"
# method_dict["AC"] = "Tailored MRM-profile_Acyl-carnitines_optimized_LIPID PHASE_04-18-2019.m"





method_dict["FFA"] = "MRM profiling_FA_05-10-2022.m"


method_dict["TAG14_0"] = "MRM profiling_TAG14_0_05-10-2022.m"
method_dict["TAG16_0"] = "MRM profiling_TAG16_0_05-10-2022.m"
method_dict["TAG16_1"] = "MRM profiling_TAG16_1_05-10-2022.m"
method_dict["TAG18_0"] = "MRM profiling_TAG18_0_05-10-2022.m"
method_dict["TAG18_1"] = "MRM profiling_TAG18_1_05-10-2022.m"
method_dict["TAG18_2"] = "MRM profiling_TAG18_2_05-10-2022.m"
method_dict["TAG18_3"] = "MRM profiling_TAG18_3_05-10-2022.m"
method_dict["TAG20_4"] = "MRM profiling_TAG20_4_05-10-2022.m"
method_dict["TAG22_5"] = "MRM profiling_TAG22_5_05-10-2022.m"
method_dict["TAG22_6"] = "MRM profiling_TAG22_6_05-10-2022.m"





method_dict["PG"] = "MRM profiling_PG_07-08-2022.m"
method_dict["PCandSM"] = "MRM profiling_LPC_PC-O_PC-P_PC_05-04-2022.m"
method_dict["PE"] = "MRM profiling_LPE_PE-O_PE-P_PE_05-17-2022.m"
method_dict["PI"] = "MRM profiling_PI_07-08-2022.m"
method_dict["PS"] = "MRM profiling_PS_07-08-2022.m"

# method_dict["cardio"] = "Tailored MRM-profiling_Cardiolipin__NEG_01-06-2021.m"

method_dict["DAG16_0"] = "MRM profiling_DAG16_0_06-22-2022.m"
method_dict["DAG16_1"] = "MRM profiling_DAG16_1_06-22-2022.m"
method_dict["DAG18_0"] = "MRM profiling_DAG18_0_06-22-2022.m"
method_dict["DAG18_1"] = "MRM profiling_DAG18_1_06-22-2022.m"
method_dict["DAG18_2"] = "MRM profiling_DAG18_2_06-22-2022.m"




method_dict["CE"] = "MRM profiling_cholesteryl esters_07-08-2022.m"
method_dict["CER"] = "MRM profiling_Ceramides_09-13-2022.m"

method_dict["AC"] = "MRM profiling_Acylcarnitines_05-05-2022.m"




#print(method_dict["FFA"])
#print(method_dict["PCandSM"])
#print(method_dict["PG"])
#print(method_dict["PE"])





sample_types = []
lipid_types = []
sample_positions = []



# file1 = open('/content/drive/MyDrive/automateMRMs/sample_names.csv', 'r')
# Lines = file1.readlines()


data_frame_read = pd.read_csv('/content/drive/MyDrive/automateMRMs/sample_names.csv')

sample_types = data_frame_read["Sample Name"]
sample_positions = data_frame_read["Position"]

# for line in Lines:
#     # print(line.strip())
#     sample_types.append(line.strip())

file2 = open('/content/drive/MyDrive/automateMRMs/lipid_classes.csv', 'r')
Lines2 = file2.readlines()

for line in Lines2:
    # print(line.strip())
    lipid_types.append(line.strip())


#file3 = open('/content/drive/MyDrive/automateMRMs/sample_positions.csv', 'r')
#Lines3 = file3.readlines()

# for line in Lines3:
#     # print(line.strip())
#     sample_positions.append(line.strip())








#print(sample_types)
#print(lipid_types)

# exit()

counter = 0

worklist_to_make = []
#print(lipid_types)

for q in range(1,(biological_replicates+1)):
    
    replicate_level = "_N"+str(q)+"_"+run_information+"_"+date_taken 
    for i in lipid_types:

        for index,j in enumerate(sample_types):
            sample_no_lipid = j+replicate_level

            #print(lipid_types)
            #print(i)
            # print(i)
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
                        # counter = counter +1

                else:

                    clean_file = clean_file_data_save+str(counter)+".d"
                    worklist_to_make.append([clean_solution_position,clean_file_method ,clean_file ])
                    # counter = counter +1
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

with open('/content/drive/MyDrive/automateMRMs/'+file_name_2_save+'.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    for i in worklist_to_make:
        spamwriter.writerow(i)


