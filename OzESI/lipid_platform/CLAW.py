#Import all the necessary libraries
import pymzml
import csv
import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
import re
import plotly.express as px
from collections import defaultdict
import plotly.io as pio
import json
import plotly.graph_objs as go
import matplotlib.colors as mcolors
import json
import ipywidgets as widgets
import warnings
import time
import shutil


def pre_parsing_setup():
    # Defaults
    data_base_name_location = 'Databases/lipid_database/Lipid_Database.xlsx'
    Project = './Projects/'
    Project_Name = 'test_py'
    Project_Folder_data = Project + Project_Name + '/o3on/'
    Project_results = Project + Project_Name + '/results/'
    file_name_to_save = 'test_py'
    tolerance = 0.3
    remove_std = True
    save_data = True
    
    # Ask if user wants to use all default settings
    use_defaults = input("Use default settings? (y/n): ").strip().lower()

    # If not using defaults, ask for custom input for each variable
    if use_defaults != 'y':
        custom_db = input(f"Use default database location '{data_base_name_location}'? (y/n): ").strip().lower()
        if custom_db == 'n':
            data_base_name_location = input("Enter custom database location: ").strip()

        custom_project = input(f"Use default project '{Project}'? (y/n): ").strip().lower()
        if custom_project == 'n':
            Project = input("Enter custom project: ").strip()

        custom_project_name = input(f"Use default project name '{Project_Name}'? (y/n): ").strip().lower()
        if custom_project_name == 'n':
            Project_Name = input("Enter custom project name: ").strip()

        Project_Folder_data = Project + Project_Name + '/mzml/o3on/'
        Project_results = Project + Project_Name + '/results/'

        custom_folder_save = input(f"Use default project results '{Project_results}'? (y/n): ").strip().lower()
        if custom_folder_save == 'n':
            Project_results = input("Enter custom folder name for project results: ").strip()

        custom_file_save = input(f"Use default file name to save results '{file_name_to_save}'? (y/n): ").strip().lower()
        if custom_file_save == 'n':
            file_name_to_save = input("Enter custom file name to save results: ").strip()

        custom_tolerance = input(f"Use default tolerance '{tolerance}'? (y/n): ").strip().lower()
        if custom_tolerance == 'n':
            tolerance = float(input("Enter custom tolerance: ").strip())

        custom_remove_std = input(f"Use default remove_std setting '{remove_std}'? (y/n): ").strip().lower()
        if custom_remove_std == 'n':
            remove_std = input("Enter custom remove_std (True/False): ").strip().lower() == 'true'

        custom_save_data = input(f"Use default save_data setting '{save_data}'? (y/n): ").strip().lower()
        if custom_save_data == 'n':
            save_data = input("Enter custom save_data (True/False): ").strip().lower() == 'true'

    # Return all configurations as a dictionary
    # Return all configurations as a dictionary
    configs = {
        "data_base_name_location": data_base_name_location,
        "Project": Project,
        "Project_Name": Project_Name,
        "Project_Folder_data": Project_Folder_data,
        "Project_results": Project_results,
        "file_name_to_save": file_name_to_save,
        "tolerance": tolerance,
        "remove_std": remove_std,
        "save_data": save_data
    }
    for key, value in configs.items():
        print(f"{key}: {value}")
    return data_base_name_location, Project_Folder_data, Project_results, file_name_to_save, tolerance

###All functions


def read_mrm_list(filename,remove_std = True):
    mrm_list_new = pd.read_excel(filename, sheet_name=None)
    mrm_list_new = pd.concat(mrm_list_new, ignore_index=True)
    mrm_list_offical = mrm_list_new[['Compound Name', 'Parent Ion', 'Product Ion', 'Class']]
    # Add underscore to middle of columns names
    mrm_list_offical.columns = mrm_list_offical.columns.str.replace(' ', '_')
    # Round Parent Ion and Product Ion to 1 decimal place
    mrm_list_offical['Parent_Ion'] = np.round(mrm_list_offical['Parent_Ion'],1)
    mrm_list_offical['Product_Ion'] = np.round(mrm_list_offical['Product_Ion'],1)
    # Create transition column by combining Parent Ion and Product Ion with arrow between numbers
    mrm_list_offical['Transition'] = mrm_list_offical['Parent_Ion'].astype(str) + ' -> ' + mrm_list_offical['Product_Ion'].astype(str)
    # Change column compound name to lipid
    mrm_list_offical = mrm_list_offical.rename(columns={'Compound_Name': 'Lipid'})
    # Make a column called Class match lipid column to lipid types
    if remove_std == True:
        lipid_class = mrm_list_offical['Class'].unique()
        lipid_class_to_keep = ['PS','PG','CE','PC', 'DAG', 'PE', 'TAG', 'FA', 'Cer', 'CAR', 'PI','SM']
        mrm_list_offical = mrm_list_offical[mrm_list_offical['Class'].isin(lipid_class_to_keep)]
    return mrm_list_offical


# Function to create an ion dictionary from an MRM database DataFrame
def create_ion_dict(mrm_database):
    ion_dict = defaultdict(list)
    # Iterate through the rows of the MRM database DataFrame
    for index, row in mrm_database.iterrows():
        # Add a tuple with Lipid and Class to the ion dictionary using Parent_Ion and Product_Ion as the key
        ion_dict[(row['Parent_Ion'], row['Product_Ion'])].append((row['Lipid'], row['Class']))
    return ion_dict

### New way to parse OzESI data
OzESI_time_df = pd.DataFrame(columns=['Lipid', 'Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Class', 'Sample_ID', 'Retention_Time', 'OzESI_Intensity'])


def mzml_parser(file_name):
    global OzESI_time_df  # Declare OzESI_time_df as a global variable
    
    rows = []
    ozesi_rows = []
    
    data_folder = os.listdir(file_name)
    data_folder.sort()
    path_to_mzml_files = file_name

    for file in data_folder:
        if file.endswith('.mzML'):
            run = pymzml.run.Reader(path_to_mzml_files + file, skip_chromatogram=False)
            q1_mz = 0
            q3_mz = 0

            for spectrum in run:
                for element in spectrum.ID.split(' '):
                    
                    if 'Q1' in element:
                        q1 = element.split('=')
                        q1_mz = np.round(float(q1[1]), 1)
                    
                    if 'Q3' in element:
                        q3 = element.split('=')
                        q3_mz = np.round(float(q3[1]), 1)
                        
                        intensity_store = np.array([intensity for _, intensity in spectrum.peaks()])
                        intensity_sum = np.sum(intensity_store)
                        
                        transition = f"{q1_mz} -> {q3_mz}"
                        sample_id = file[:-5]
                        
                        rows.append({
                            'Parent_Ion': q1_mz,
                            'Product_Ion': q3_mz,
                            'Intensity': intensity_sum,
                            'Transition': transition,
                            'Sample_ID': sample_id
                        })
                        
                        for time, intensity in spectrum.peaks():
                            ozesi_rows.append({
                                'Parent_Ion': q1_mz,
                                'Product_Ion': q3_mz,
                                'Retention_Time': time,
                                'OzESI_Intensity': intensity,
                                'Sample_ID': sample_id,
                                'Transition': transition
                            })

    df = pd.DataFrame(rows)
    OzESI_time_df = pd.DataFrame(ozesi_rows)
    print('Finished parsing mzML files\n')
    return df, OzESI_time_df


# Function to check if the absolute difference between two values is within a given tolerance
def within_tolerance(a, b, tolerance=0.3):
    return abs(a - b) <= tolerance

# Function to match the ions in a DataFrame row with the ions in an ion dictionary
def match_ions(row, ion_dict, tolerance=0.3):
    ions = (row['Parent_Ion'], row['Product_Ion'])
    matched_lipids = []
    matched_classes = []

    # Iterate through the ion dictionary
    for key, value in ion_dict.items():
        # Check if both the Parent_Ion and Product_Ion values are within the specified tolerance
        if within_tolerance(ions[0], key[0], tolerance) and within_tolerance(ions[1], key[1], tolerance):
            # If within tolerance, extend the matched_lipids and matched_classes lists with the corresponding values
            matched_lipids.extend([match[0] for match in value])
            matched_classes.extend([match[1] for match in value])

    # If any matches were found, update the Lipid and Class columns in the row
    if matched_lipids and matched_classes:
        row['Lipid'] = ' | '.join(matched_lipids)
        row['Class'] = ' | '.join(matched_classes)

    return row

####Combined functions for Matching

def match_lipids_parser(mrm_database,df, tolerance=0.3):
    ion_dict = create_ion_dict(mrm_database)
    # Assuming you have the df DataFrame to apply the match_ions function
    df_matched = df.apply(lambda row: match_ions(row, ion_dict=ion_dict, tolerance=tolerance), axis=1)


    # df_matched = df_matched.dropna()
    
    return df_matched


def save_dataframe(df, Project_results, file_name_to_save, max_attempts=5):
    folder_path = f'data_results/data/data_matching/{Project_results}'
    os.makedirs(folder_path, exist_ok=True)

    for i in range(max_attempts):
        file_path = f'{folder_path}/{file_name_to_save}.csv'
        if not os.path.isfile(file_path):
            df.to_csv(file_path, index=False)
            print(f"Saved DataFrame to {file_path}")
            break
    else:
        print(f"Failed to save DataFrame after {max_attempts} attempts.")
        return None



# def full_parse(data_base_name_location, Project_Folder_data, Project_results, file_name_to_save, tolerance, remove_std=True, save_data=False):
#     mrm_database = read_mrm_list(data_base_name_location, remove_std=remove_std)
#     df = mzml_parser(Project_Folder_data)
#     df_matched = match_lipids_parser(mrm_database, df, tolerance=tolerance)
    
#     if save_data:
#         save_dataframe(df_matched, Project_results, file_name_to_save)

#     return df_matched


def full_parse(data_base_name_location, Project_Folder_data, Project_results, file_name_to_save, tolerance, remove_std=True, save_data=False):
    mrm_database = read_mrm_list(data_base_name_location, remove_std=remove_std)
    
    # Capture both dataframes returned by mzml_parser
    df, OzESI_time_df = mzml_parser(Project_Folder_data)
    
    # Use the df dataframe in match_lipids_parser
    df_matched = match_lipids_parser(mrm_database, df, tolerance=tolerance)
    
    if save_data:
        save_dataframe(df_matched, Project_results, file_name_to_save)

    # If you need to use OzESI_time_df elsewhere in this function, you can do so here

    return df_matched, OzESI_time_df




def filter_rt(df):
    """
    Filters the DataFrame based on retention times and aggregates by max intensity for unique 'Sample_ID' and 'Transition' combinations.
    
    Parameters:
        df (pd.DataFrame): Input DataFrame with columns 'Retention_Time' and 'OzESI_Intensity'.
        
    Returns:
        pd.DataFrame: Filtered and aggregated DataFrame.
    """
    # Filter based on retention time
    filtered_df = df[(df['Retention_Time'] > 10.0) & (df['Retention_Time'] < 16.0)].copy()

    # Round the values
    filtered_df['Retention_Time'] = filtered_df['Retention_Time'].round(2)
    filtered_df['OzESI_Intensity'] = filtered_df['OzESI_Intensity'].round(0)

    # Aggregate by max intensity for unique combinations of 'Sample_ID' and 'Transition'
    filtered_df = filtered_df.groupby(['Sample_ID', 'Transition']).apply(
        lambda x: x.loc[x['OzESI_Intensity'].idxmax()]).reset_index(drop=True)

    return filtered_df

def concat_dataframes(df_matched, filtered_df):
    """
    Concatenates two DataFrames along the columns.
    
    Parameters:
        df_matched (pd.DataFrame): First DataFrame.
        filtered_df (pd.DataFrame): Second DataFrame, only the 'Retention_Time' and 'OzESI_Intensity' columns will be used.
        
    Returns:
        pd.DataFrame: Concatenated DataFrame.
    """
    return pd.concat([df_matched, filtered_df[['Retention_Time', 'OzESI_Intensity']]], axis=1)

def DB_Position_df(df_matched_2, OzESI_list=[7,9,12]):
    """
    Creates a new DataFrame to store the DB_Position and Aldehyde_Ion values,
    and calculate n-i values for the given OzESI_list.
    
    Parameters:
        df_matched_2 (pd.DataFrame): Input DataFrame.
        OzESI_list (list): List of OzESI positions.
        
    Returns:
        pd.DataFrame: Modified DataFrame with new calculated columns.
    """
    # Create a new DataFrame to store the DB_Position and Aldehyde_Ion values
    df_DB_position = pd.DataFrame(columns=['DB_Position','Aldehyde_Ion'])

    # Loop over the range of DB_Position values and calculate the corresponding Aldehyde_Ion values
    for i in range(3, 21):
        df_DB_position.loc[i, 'DB_Position'] = i
        df_DB_position.loc[i, 'Aldehyde_Ion'] = 26 + (14 * (i-3))

    # Loop through OzESI_list
    for i in OzESI_list:
        # Retrieve the Aldehyde_Ion value for the current DB_Position
        aldehyde_ion = df_DB_position.loc[df_DB_position["DB_Position"] == i, "Aldehyde_Ion"].values[0]

        # Calculate n-i values
        df_matched_2["n-{}".format(i)] = df_matched_2["Parent_Ion"] - aldehyde_ion

    return df_matched_2


def within_tolerance(a, b, tolerance=0.3):
    return abs(a - b) <= tolerance
columns = [
    'Lipid', 'Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Class',
    'Sample_ID', 'Retention_Time', 'Intensity_OzESI', 'Mean_Retention_Time',
    'Mean_Intensity_OzESI', 'n-7', 'n-9', 'n-12', 'Labels'
]
df_OzESI_n = pd.DataFrame(columns=columns)

#Function to add lipid name
def add_lipid_info(df_matched_2, OzESI_list, tolerance=0.3):
    df_test = df_matched_2.copy()
    df_test_2 = df_matched_2.copy()
    global df_OzESI_n
    
    for i in OzESI_list:
        df_test['n-' + str(i)] = df_test['n-' + str(i)].astype(float)
    
    for i in range(len(df_test)):
        if pd.isna(df_test.loc[i, 'Lipid']):
            parent_ion = df_test.loc[i, 'Parent_Ion']
            
            for j in range(len(df_test)):
                row_data = df_test.loc[j].copy()
                if within_tolerance(parent_ion, row_data['n-7'], tolerance) and isinstance(row_data['Lipid'], str):
                    df_test.loc[i, 'Lipid'] = row_data['Lipid']
                    df_test.loc[i, 'Labels'] = 'n-7' + row_data['Labels']
                    
                    # Append to df_test_2
                    appended_row = df_test.loc[i].copy()
                    appended_row['Labels'] = 'n-7' + row_data['Labels']
                    df_test_2 = df_test_2.append(appended_row, ignore_index=True)
                    
                elif within_tolerance(parent_ion, row_data['n-9'], tolerance) and isinstance(row_data['Lipid'], str):
                    df_test.loc[i, 'Lipid'] = row_data['Lipid']
                    df_test.loc[i, 'Labels'] = 'n-9' + row_data['Labels']
                    
                    # Append to df_test_2
                    appended_row = df_test.loc[i].copy()
                    appended_row['Labels'] = 'n-9' + row_data['Labels']
                    df_test_2 = df_test_2.append(appended_row, ignore_index=True)
                    
                elif within_tolerance(parent_ion, row_data['n-12'], tolerance) and isinstance(row_data['Lipid'], str):
                    df_test.loc[i, 'Lipid'] = row_data['Lipid']
                    df_test.loc[i, 'Labels'] = 'n-12' + row_data['Labels']
                    
                    # Append to df_test_2
                    appended_row = df_test.loc[i].copy()
                    appended_row['Labels'] = 'n-12' + row_data['Labels']
                    df_test_2 = df_test_2.append(appended_row, ignore_index=True)
    
    df_test_2.dropna(subset=['Lipid'], inplace=True)
    return df_test_2




def calculate_intensity_ratio(df):
    # Create a new column for ratios
    df['Ratios'] = pd.Series(dtype='float64')

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        lipid = row['Lipid']
        label = row['Labels']
        intensity = row['OzESI_Intensity']
        sample_id = row['Sample_ID']

        # Check if the label is n-9
        if label == 'n-9':
            # Find the corresponding row with n-7 label and same lipid name
            n7_row = df[(df['Lipid'] == lipid) & (df['Labels'] == 'n-7')& (df['Sample_ID'] == sample_id)]

            # If a matching row is found, calculate the intensity ratio
            if not n7_row.empty:
                n7_intensity = n7_row['OzESI_Intensity'].values[0]
                ratio = intensity / n7_intensity

                # Assign the ratio to the 'Ratios' column
                df.at[index, 'Ratios'] = ratio

    return df

#filtering df_matched for name but not removing the n-7 ratios
def sort_by_second_tg(lipid):
    tgs = lipid.split(',')
    if len(tgs) > 1:
        return tgs[1]
    else:
        return lipid

def filter_highest_ratios(df):
    # Sort the DataFrame by ratios in descending order
    df_sorted = df.sort_values(by='Ratios', ascending=False)

    # Drop duplicates keeping the first occurrence (highest ratio)
    df_filtered = df_sorted.drop_duplicates(subset=['Sample_ID', 'Lipid','Labels'], keep='first')
    df_filtered = df_filtered.sort_values(by=['Sample_ID', 'Lipid'], ascending=[True, True])

    return df_filtered
