import csv
import pandas as pd
import os

# Read the Path_info.csv file to get paths and methods information
data_frame_read = pd.read_csv('./Path_info.csv')

# Extract specific paths and method information from the DataFrame
path_to_methods = data_frame_read["path_to_methods"][0]
path_to_save_data = data_frame_read["path_to_save_data"][0]
clean_file_method = data_frame_read["clean_file_method"][0]

# Set the directory where cleaned file data will be saved
clean_file_data_save = path_to_save_data

# Get user inputs for cleaning parameters and experiment details
clean_every_n = int(input("Clean Every N: "))
n_cleans = int(input("How many cleans every N: "))
n_cleans_after_lipid_class = int(input("How many cleans every lipid class: "))
clean_solution_position = input("Clean Solution Position EX p1-a9: ")
biological_replicates = int(input("How Many Technical Replicates?: "))
run_information = input("Run Information (Project Name) EX: Burda_Lab: ")
date_taken = input("Run Date in format 3_3_23: ")

# Initialize a dictionary to store method information
method_dict = {}

# Read methods.csv to populate the method_dict with annotation-method name pairs
with open('methods.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        method_dict[row['annotation']] = row['method name']

# Print the dictionary to verify correct loading of methods
print(method_dict)

# Initialize lists to store sample types, lipid types, and sample positions
sample_types = []
lipid_types = []
sample_positions = []

# Read sample_names.csv to get sample types and positions
data_frame_read = pd.read_csv('./sample_names.csv')
sample_types = data_frame_read["Sample Name"]
sample_positions = data_frame_read["Position"]

# Read lipid_classes.csv to get the list of lipid types
with open('./lipid_classes.csv', 'r') as file2:
    Lines2 = file2.readlines()
    for line in Lines2:
        lipid_types.append(line.strip())

# Initialize a counter for tracking the number of operations
counter = 0

# Initialize a list to store the worklist entries
worklist_to_make = []

# Loop through each biological replicate
for q in range(1, (biological_replicates + 1)):
    replicate_level = "_N" + str(q) + "_" + run_information + "_" + date_taken

    # Loop through each lipid type
    for i in lipid_types:
        # Loop through each sample type
        for index, j in enumerate(sample_types):
            sample_no_lipid = j + replicate_level
            sample_with_lipid = path_to_save_data + i + "_" + sample_no_lipid + ".d"
            method_to_use = path_to_methods + method_dict[i]
            counter += 1

            # Add the sample information to the worklist
            worklist_to_make.append([sample_positions[index].replace(",", "_"), 
                                     method_to_use.replace(",", "_"), 
                                     sample_with_lipid.replace(",", "_")])

            # Add cleaning operations after every 'clean_every_n' samples
            if counter % clean_every_n == 0:
                for ii in range(n_cleans):
                    clean_file = clean_file_data_save + str(counter) + (str(ii) if n_cleans > 1 else "") + ".d"
                    worklist_to_make.append([clean_solution_position, clean_file_method, clean_file])

        # Add cleaning operations after each lipid class
        for ii in range(n_cleans_after_lipid_class):
            clean_file = clean_file_data_save + str(counter) + (str(ii) if n_cleans_after_lipid_class > 1 else "") + "post_lipid_class_cleaning" + ".d"
            worklist_to_make.append([clean_solution_position, clean_file_method, clean_file])

# Ensure the worklists directory exists
os.makedirs('./worklists', exist_ok=True)

# Get the file name for saving the worklist from the user
file_name_2_save = input("Save File Name (No file extension)?")

# Save the generated worklist to a CSV file
with open('./worklists/' + file_name_2_save + '.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    for i in worklist_to_make:
        spamwriter.writerow(i)
