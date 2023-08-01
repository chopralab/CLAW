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
import os



def create_dataframes():
    time_and_intensity_df = pd.DataFrame(columns=['Time', 'Intensity'])
    master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
    OzESI_time_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])
    
    return time_and_intensity_df, master_df, OzESI_time_df



def pre_parsing_setup(data_base_name_location, Project, Project_Name, Project_Folder_data, Project_results, file_name_to_save, tolerance, remove_std, save_data):
    """
    A function to setup and check for the necessary project directories. 
    It also prints and returns the received configurations. 

    :param data_base_name_location: The path of the database name location.
    :param Project: The project path.
    :param Project_Name: The name of the project.
    :param Project_Folder_data: The project data folder path.
    :param Project_results: The path where project results are stored.
    :param file_name_to_save: The name of the file where data is to be saved.
    :param tolerance: The accepted tolerance level.
    :param remove_std: A boolean indicating whether or not to remove standard deviations.
    :param save_data: A boolean indicating whether or not to save the data.

    :return: The given configurations as a dictionary.
    """

    # Check and create folders if they do not exist
    os.makedirs(os.path.dirname(data_base_name_location), exist_ok=True)
    os.makedirs(Project, exist_ok=True)
    os.makedirs(Project_Folder_data, exist_ok=True)
    os.makedirs(Project_results, exist_ok=True)

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
    return data_base_name_location, Project_Folder_data, Project_results, file_name_to_save, tolerance, remove_std, save_data


def read_mrm_list(filename,remove_std = True):
    """
    A function that reads a Multi Reaction Monitoring (MRM) list from an Excel file, formats the data, 
    and filters out certain lipid classes if required.

    :param filename: The path of the Excel file containing the MRM list.
    :param remove_std: A boolean indicating whether or not to filter out certain lipid classes. Defaults to True.

    :return: The formatted and filtered MRM list as a pandas DataFrame.
    """
    
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
    
    # Filter the DataFrame based on lipid class, if required
    if remove_std == True:
        lipid_class_to_keep = ['PS','PG','CE','PC', 'DAG', 'PE', 'TAG', 'FA', 'Cer', 'CAR', 'PI','SM']
        mrm_list_offical = mrm_list_offical[mrm_list_offical['Class'].isin(lipid_class_to_keep)]
    
    return mrm_list_offical


from collections import defaultdict
import os
import numpy as np
import pandas as pd
import pymzml


def create_ion_dict(mrm_database):
    """
    Creates a dictionary of ions from an MRM database DataFrame.
    
    :param mrm_database: DataFrame containing MRM database information.
    
    :return: A dictionary with ion pairs as keys, and a list of tuples containing corresponding lipid and class as values.
    """
    ion_dict = defaultdict(list)
    for index, row in mrm_database.iterrows():
        ion_dict[(row['Parent_Ion'], row['Product_Ion'])].append((row['Lipid'], row['Class']))
    return ion_dict



# Declare the DataFrame globally if it's used across multiple functions
time_and_intensity_df = pd.DataFrame(columns=['Time', 'Intensity'])
master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
OzESI_time_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])

def mzml_parser(file_path):
    global master_df
    global OzESI_time_df
    global time_and_intensity_df
    
    rows = []
    ozesi_rows = []
    
    run = pymzml.run.Reader(file_path, skip_chromatogram=False)
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

                #####checking (TG 52:2) 876.6 -> 577.6 specifically for the OzESI data
                #####################################################################
                # if within_tolerance(q1_mz, 876.6) and within_tolerance(q3_mz, 577.6):
                #     for time, intensity in spectrum.peaks():
                #         time_and_intensity_df = time_and_intensity_df.append({
                #             'Time': time,
                #             'Intensity': intensity
                #         }, ignore_index=True)
                #     times, intensities = zip(*spectrum.peaks())
                #     plt.plot(times, intensities)
                #     plt.xlabel('Time')
                #     plt.ylabel('Intensity')
                #     plt.title('Chromatogram for 876.6 -> 577.6')
                #     plt.show()
                ####################################################################

                intensity_store = np.array([intensity for _, intensity in spectrum.peaks()])
                intensity_sum = np.sum(intensity_store)
                
                transition = f"{q1_mz} -> {q3_mz}"
                sample_id = os.path.basename(file_path)[:-5]
                
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
    # OzESI_time_df = pd.DataFrame(ozesi_rows)
    # Append the new rows to the existing global OzESI_time_df
    OzESI_time_df = OzESI_time_df.append(pd.DataFrame(ozesi_rows), ignore_index=True)
    master_df = master_df.append(df, ignore_index=True)
    print(f'Finished parsing mzML file: {file_path}\n')


def mzml_parser_batch(folder_name):
    global master_df
    global time_and_intensity_df
    
    data_folder = os.listdir(folder_name)
    data_folder.sort()

    for file in data_folder:
            if file.endswith('.mzML'):
                file_path = os.path.join(folder_name, file)
                mzml_parser(file_path)
    
    print('Finished parsing all mzML files\n')




def within_tolerance(a, b, tolerance=0.3):
    """
    Checks if the absolute difference between two values is within a given tolerance.
    
    :param a: First value to compare.
    :param b: Second value to compare.
    :param tolerance: The acceptable difference between the two values. Defaults to 0.3.
    
    :return: Boolean indicating whether the difference is within the given tolerance.
    """
    return abs(a - b) <= tolerance


def match_ions(row, ion_dict, tolerance=0.3):
    """
    Matches the ions in a DataFrame row with the ions in an ion dictionary.
    
    :param row: A DataFrame row containing 'Parent_Ion' and 'Product_Ion' columns.
    :param ion_dict: A dictionary of ion pairs and their corresponding lipid and class information.
    :param tolerance: The acceptable difference between ion values to be considered a match. Defaults to 0.3.
    
    :return: The original row updated with matched lipid and class information if matches were found.
    """
    ions = (row['Parent_Ion'], row['Product_Ion'])
    matched_lipids = []
    matched_classes = []

    for key, value in ion_dict.items():
        if within_tolerance(ions[0], key[0], tolerance) and within_tolerance(ions[1], key[1], tolerance):
            matched_lipids.extend([match[0] for match in value])
            matched_classes.extend([match[1] for match in value])

    if matched_lipids and matched_classes:
        row['Lipid'] = ' | '.join(matched_lipids)
        row['Class'] = ' | '.join(matched_classes)

    return row


def match_lipids_parser(mrm_database, df, tolerance=0.3):
    """
    Performs lipid matching by creating an ion dictionary from the MRM database and applying the match_ions function to each row of a DataFrame.
    
    :param mrm_database: DataFrame containing MRM database information.
    :param df: DataFrame containing ion information to be matched.
    :param tolerance: The acceptable difference between ion values to be considered a match. Defaults to 0.3.
    
    :return: DataFrame with matched lipid and class information if matches were found.
    """
    ion_dict = create_ion_dict(mrm_database)
    df_matched = df.apply(lambda row: match_ions(row, ion_dict=ion_dict, tolerance=tolerance), axis=1)
    return df_matched



def save_dataframe(df, Project_results, file_name_to_save, max_attempts=5):
    """
    Saves a given DataFrame to a CSV file within a specified directory.
    
    :param df: DataFrame to be saved.
    :param Project_results: The project directory to save results in.
    :param file_name_to_save: The desired filename for the saved DataFrame.
    :param max_attempts: The maximum number of attempts to save the DataFrame. Defaults to 5.
    
    :return: None
    """
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



def full_parse(data_base_name_location, Project_Folder_data, Project_results, file_name_to_save, tolerance, remove_std=True, save_data=False, batch_processing=True):
    """
    Performs the complete parsing and data matching process for given inputs.
    
    :param data_base_name_location: Location of the MRM database file to be read.
    :param Project_Folder_data: The project folder containing data to be parsed (folder path or single file path).
    :param Project_results: The project directory to save results in.
    :param file_name_to_save: The desired filename for the saved DataFrame.
    :param tolerance: The acceptable difference between ion values to be considered a match.
    :param remove_std: A flag to indicate whether standard lipid classes should be removed. Defaults to True.
    :param save_data: A flag to indicate whether the matched data should be saved. Defaults to False.
    :param batch_processing: A flag to indicate whether to process a batch of files (directory) or a single file. Defaults to True.
    
    :return: Tuple containing matched DataFrame and OzESI DataFrame.
    """
    global master_df
    mrm_database = read_mrm_list(data_base_name_location, remove_std=remove_std)
    
    if batch_processing:
        mzml_parser_batch(Project_Folder_data)  # Call batch function
    else:
        mzml_parser(Project_Folder_data)  # Call single file function
    
    df_matched = match_lipids_parser(mrm_database, master_df, tolerance=tolerance)
    
    if save_data:
        save_dataframe(df_matched, Project_results, file_name_to_save)

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
    filtered_df = df[(df['Retention_Time'] > 10.0) & (df['Retention_Time'] < 20.0)].copy()

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


###delete this?
# columns = [
#     'Lipid', 'Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Class',
#     'Sample_ID', 'Retention_Time', 'Intensity_OzESI', 'Mean_Retention_Time',
#     'Mean_Intensity_OzESI', 'n-7', 'n-9', 'n-12', 'db_pos'
# ]
# df_OzESI_n = pd.DataFrame(columns=columns)


def add_lipid_info(df_matched_2, OzESI_list, tolerance=0.3):
    """
    Adds lipid information to the data frame based on matched ions within a certain tolerance.

    :param df_matched_2: DataFrame containing matched lipids and ion data.
    :param OzESI_list: List of integer values representing the positions in the OzESI list to be checked.
    :param tolerance: The acceptable difference between ion values to be considered a match.

    :return: Updated DataFrame with added lipid information.
    """
    df_test = df_matched_2.copy()  # create a copy of the input DataFrame to avoid modifying the original data
    df_test_2 = df_matched_2.copy()  # additional copy for final output

    # Iterate over the provided OzESI_list and convert respective column values to float
    for i in OzESI_list:
        df_test['n-' + str(i)] = df_test['n-' + str(i)].astype(float)

    # Iterate over the rows of the DataFrame
    for i in range(len(df_test)):
        # If the 'Lipid' column value is NaN for the current row
        if pd.isna(df_test.loc[i, 'Lipid']):
            parent_ion = df_test.loc[i, 'Parent_Ion']

            # Iterate over the rows again to find a match
            for j in range(len(df_test)):
                row_data = df_test.loc[j].copy()
                
                # If the parent ion is within tolerance and the Lipid column is a string
                for n in [7, 9, 12]:
                    if within_tolerance(parent_ion, row_data[f'n-{n}'], tolerance) and isinstance(row_data['Lipid'], str):
                        df_test.loc[i, 'Lipid'] = row_data['Lipid']
                        df_test.loc[i, 'db_pos'] = f'n-{n}' + row_data['db_pos']

                        # Append to df_test_2
                        appended_row = df_test.loc[i].copy()
                        appended_row['db_pos'] = f'n-{n}' + row_data['db_pos']
                        df_test_2 = df_test_2.append(appended_row, ignore_index=True)

    # Drop rows in df_test_2 where 'Lipid' column value is NaN
    df_test_2.dropna(subset=['Lipid'], inplace=True)

    return df_test_2


def calculate_intensity_ratio(df):
    """
    Calculates the intensity ratio for each lipid in the DataFrame, based on their 'OzESI_Intensity'.

    :param df: DataFrame containing lipid information and intensity values.
    :return: Updated DataFrame with added 'Ratios' column.
    """
    # Create a new column for ratios
    df['Ratios'] = pd.Series(dtype='float64')

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        lipid = row['Lipid']
        label = row['db_pos']
        intensity = row['OzESI_Intensity']
        sample_id = row['Sample_ID']

        # Check if the label is n-9
        if label == 'n-9':
            # Find the corresponding row with n-7 label and same lipid name and Sample_ID
            n7_row = df[(df['Lipid'] == lipid) & (df['db_pos'] == 'n-7')& (df['Sample_ID'] == sample_id)]

            # If a matching row is found, calculate the intensity ratio
            if not n7_row.empty:
                n7_intensity = n7_row['OzESI_Intensity'].values[0]
                ratio = intensity / n7_intensity

                # Assign the ratio to the 'Ratios' column
                df.at[index, 'Ratios'] = ratio

    return df


def sort_by_second_tg(lipid):
    """
    Helper function that sorts lipid names by second triglyceride, if present.

    :param lipid: string of lipid names.
    :return: Second triglyceride if present, else returns original lipid.
    """
    if pd.isna(lipid):
        return lipid
    tgs = lipid.split(',')
    if len(tgs) > 1:
        return tgs[1]
    else:
        return lipid



def filter_highest_ratios(df):
    """
    Filters the DataFrame to keep only rows with the highest ratio value for each unique Sample_ID and lipid.

    :param df: DataFrame containing lipid information and intensity ratios.
    :return: Filtered DataFrame.
    """
    # Sort the DataFrame by ratios in descending order
    df_sorted = df.sort_values(by='Ratios', ascending=False)

    # Drop duplicates keeping the first occurrence (highest ratio)
    df_filtered = df_sorted.drop_duplicates(subset=['Sample_ID', 'Lipid','db_pos'], keep='first')
    df_filtered = df_filtered.sort_values(by=['Sample_ID', 'Lipid'], ascending=[True, True])

    return df_filtered
