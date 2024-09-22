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




def create_analysis_dataframes():
    """
    Creates and returns three DataFrames for storing time-intensity data, master data, and OzESI time data.

    Returns:
        pd.DataFrame: DataFrame for storing time and intensity values.
        pd.DataFrame: Master DataFrame for storing parent ion, product ion, intensity, transition, and sample ID.
        pd.DataFrame: DataFrame for storing OzESI (OzID Electron Spray Ionization) data including parent ion, product ion, retention time, intensity, sample ID, and transition.
    """
    time_intensity_dataframe = pd.DataFrame(columns=['Time', 'Intensity'])
    master_lipid_dataframe = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
    OzESI_time_dataframe = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])
    
    return time_intensity_dataframe, master_lipid_dataframe, OzESI_time_dataframe




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


def read_mrm_list(filename, remove_std=True, deuterated=False):
    """
    Reads a Multiple Reaction Monitoring (MRM) lipid database from an Excel file and processes the data.

    Parameters:
        filename (str): The path to the Excel file containing the MRM list.
        remove_std (bool): Whether to exclude lipid classes not in a predefined list.
        deuterated (bool): Whether to adjust ion values for deuterated lipids.

    Returns:
        pd.DataFrame: A DataFrame containing processed MRM lipid data.
    """
    # Concatenate all sheets from the Excel file into one DataFrame
    raw_mrm_data = pd.read_excel(filename, sheet_name=None)
    concatenated_mrm_data = pd.concat(raw_mrm_data, ignore_index=True)

    # Extract the required columns and rename them
    lipid_MRM_data = concatenated_mrm_data[['Compound Name', 'Parent Ion', 'Product Ion', 'Class']]
    lipid_MRM_data.columns = lipid_MRM_data.columns.str.replace(' ', '_')
    lipid_MRM_data['Parent_Ion'] = np.round(lipid_MRM_data['Parent_Ion'], 1)
    lipid_MRM_data['Product_Ion'] = np.round(lipid_MRM_data['Product_Ion'], 1)
    lipid_MRM_data['Transition'] = lipid_MRM_data['Parent_Ion'].astype(str) + ' -> ' + lipid_MRM_data['Product_Ion'].astype(str)
    lipid_MRM_data = lipid_MRM_data.rename(columns={'Compound_Name': 'Lipid'})

    # Optionally filter the data to keep only specific lipid classes
    if remove_std:
        lipid_classes_to_keep = ['PS', 'PG', 'CE', 'PC', 'DAG', 'PE', 'TAG', 'FA', 'Cer', 'CAR', 'PI', 'SM']
        lipid_MRM_data = lipid_MRM_data[lipid_MRM_data['Class'].isin(lipid_classes_to_keep)]

    # Optionally adjust the ion values for deuterated lipids
    if deuterated:
        lipid_MRM_data['Parent_Ion'] += 1
        lipid_MRM_data['Product_Ion'] += 1
        # Update the Transition column with the updated values
        lipid_MRM_data['Transition'] = lipid_MRM_data['Parent_Ion'].astype(str) + ' -> ' + lipid_MRM_data['Product_Ion'].astype(str)
    
    return lipid_MRM_data



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
OzESI_time_df = pd.DataFrame(columns=['Lipid','Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])

import os
import numpy as np
import pandas as pd
import pymzml
import matplotlib.pyplot as plt

# Initialize global DataFrames
time_and_intensity_df = pd.DataFrame(columns=['Time', 'Intensity'])
master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
OzESI_time_df = pd.DataFrame(columns=['Lipid', 'Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])

def mzml_parser(file_path, plot_chromatogram=False):
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
                q1_mz = float(q1[1])  # Assign correctly
                
            if 'Q3' in element:
                q3 = element.split('=')
                q3_mz = float(q3[1])  # Assign correctly
                intensity_store = np.array([intensity for _, intensity in spectrum.peaks()])
                intensity_sum = np.sum(intensity_store)
                transition = f"{q1_mz} -> {q3_mz}"
                sample_id = os.path.basename(file_path).replace('.mzML', '')  # More robust handling
        
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
    OzESI_time_df = OzESI_time_df.append(pd.DataFrame(ozesi_rows), ignore_index=True)
    master_df = master_df.append(df, ignore_index=True)
    print(f'Finished parsing mzML file: {file_path}\n')



def mzml_parser_batch(folder_name, plot_chromatogram=False):
    global master_df
    global time_and_intensity_df
    
    data_folder = os.listdir(folder_name)
    data_folder.sort()

    for file in data_folder:
        if file.endswith('.mzML'):
            file_path = os.path.join(folder_name, file)
            mzml_parser(file_path, plot_chromatogram=plot_chromatogram)  # Pass the flag here
    
    print('Finished parsing all mzML files\n')



def within_tolerance(a, b, tolerance=0.5):
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



def full_parse(data_base_name_location, 
               Project_Folder_data, 
               Project_results, 
               file_name_to_save, 
               tolerance, 
               remove_std=True, 
               save_data=False, 
               batch_processing=True,
               plot_chromatogram=False):
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
    
    # Batch processing
    if batch_processing:
        mzml_parser_batch(Project_Folder_data, plot_chromatogram=plot_chromatogram)  
    # Single file processing
    else:
        mzml_parser(Project_Folder_data, plot_chromatogram=plot_chromatogram) 
    
    df_matched = match_lipids_parser(mrm_database, master_df, tolerance=tolerance)
    
    if save_data:
        save_dataframe(df_matched, Project_results, file_name_to_save)

    return df_matched, OzESI_time_df



def filter_rt(df, min_rt=10.0, max_rt=22.0, min_intensity=100, 
                                special_min_rt=19.5, special_max_rt=21.5, 
                                special_sample_id='DegummedCanola_O3on_150gN3_02082023', 
                                special_parent_ion=794.6):
    """
    Filters the DataFrame based on retention times and aggregates by max intensity for unique 
    'Sample_ID' and 'Transition' combinations, with a special case that must meet specific criteria 
    before applying the max intensity aggregation.
    
    Parameters:
        df (pd.DataFrame): Input DataFrame with columns 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition', 'Parent_Ion'.
        min_rt (float, optional): Minimum retention time for general case. Defaults to 10.0.
        max_rt (float, optional): Maximum retention time for general case. Defaults to 22.0.
        min_intensity (float, optional): Minimum intensity value for general case. Defaults to 100.
        special_min_rt (float, optional): Minimum retention time for special case. Defaults to 19.5.
        special_max_rt (float, optional): Maximum retention time for special case. Defaults to 21.5.
        special_sample_id (str, optional): Specific sample to filter in special case. Defaults to 'DegummedCanola_O3on_150gN3_02082023'.
        special_parent_ion (float, optional): Specific parent ion to filter in special case. Defaults to 794.6.
        
    Returns:
        pd.DataFrame: Filtered and aggregated DataFrame.
    """
    
    # Apply general filter based on retention time and intensity
    general_filter = (df['Retention_Time'] > min_rt) & \
                     (df['Retention_Time'] < max_rt) & \
                     (df['OzESI_Intensity'] > min_intensity)
    
    filtered_df = df[general_filter].copy()

    # Special case filter: apply additional conditions for the special case
    special_case_filter = (
        (df['Sample_ID'] == special_sample_id) & 
        (df['Parent_Ion'] == special_parent_ion) & 
        (df['Retention_Time'] >= special_min_rt) & 
        (df['Retention_Time'] <= special_max_rt)
    )

    # Mark rows that meet the special case condition
    filtered_df['is_special_case'] = special_case_filter.astype(int)

    # Round the values for 'Retention_Time' and 'OzESI_Intensity'
    filtered_df['Retention_Time'] = filtered_df['Retention_Time'].round(2)
    filtered_df['OzESI_Intensity'] = filtered_df['OzESI_Intensity'].round(0)

    # Aggregate: for non-special cases, take max intensity; for special cases, use special filter first
    def apply_aggregation(group):
        if group['is_special_case'].sum() > 0:
            # For special case, only consider rows that meet special case criteria
            special_case_rows = group[group['is_special_case'] == 1]
            return special_case_rows.loc[special_case_rows['OzESI_Intensity'].idxmax()]
        else:
            # For all other cases, just take the max intensity
            return group.loc[group['OzESI_Intensity'].idxmax()]
    
    # Group by 'Sample_ID' and 'Transition' and apply aggregation
    result_df = filtered_df.groupby(['Sample_ID', 'Transition']).apply(apply_aggregation).reset_index(drop=True)
    
    return result_df


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

def calculate_DB_Position(df_matched_ions, db_pos_list=[7,9,12]):
    """
    Creates a new DataFrame to store the DB_Position and Aldehyde_Ion values,
    and calculate n-i values for the given db_pos_list.
    
    Parameters:
        df_matched_ions (pd.DataFrame): Input DataFrame containing matched ions.
        db_pos_list (list): List of OzESI positions to calculate n-i values.
        
    Returns:
        pd.DataFrame: Modified DataFrame with new calculated columns for n-i values.
    """
    # Create a DataFrame to store the DB_Position and corresponding Aldehyde_Ion values
    df_DB_aldehyde = pd.DataFrame(columns=['DB_Position','Aldehyde_Ion'])

    # Loop through the range of DB_Position values to calculate the corresponding Aldehyde_Ion values
    for position in range(3, 21):
        df_DB_aldehyde.loc[position, 'DB_Position'] = position
        df_DB_aldehyde.loc[position, 'Aldehyde_Ion'] = 26 + (14 * (position-3))

    # Loop through the specified db_pos_list
    for ozesi_position in db_pos_list:
        # Retrieve the corresponding Aldehyde_Ion value for the current DB_Position
        aldehyde_ion = df_DB_aldehyde.loc[df_DB_aldehyde["DB_Position"] == ozesi_position, "Aldehyde_Ion"].values[0]

        # Calculate and store the n-i value for the current OzESI position
        df_matched_ions["n-{}".format(ozesi_position)] = df_matched_ions["Parent_Ion"] - aldehyde_ion

    return df_matched_ions



def add_lipid_info(matched_dataframe, db_pos, tolerance=0.3):
    """
    Adds lipid information to the data frame based on matched ions within a certain tolerance.

    :param matched_dataframe: DataFrame containing matched lipids and ion data.
    :param db_pos: List of integer values representing the positions in the OzESI list to be checked.
    :param tolerance: The acceptable difference between ion values to be considered a match.

    :return: Updated DataFrame with added lipid information.
    """
    working_dataframe = matched_dataframe.copy()  # Create a copy for processing
    final_dataframe = matched_dataframe.copy()    # Create a copy for final output

    # Convert respective column values to float for given db_pos
    for position in db_pos:
        working_dataframe['n-' + str(position)] = working_dataframe['n-' + str(position)].astype(float)

    # Iterate over the rows of the DataFrame to match lipids
    for i in range(len(working_dataframe)):
        if pd.isna(working_dataframe.loc[i, 'Lipid']):
            parent_ion = working_dataframe.loc[i, 'Parent_Ion']

            # Look for matching ions within tolerance
            for j in range(len(working_dataframe)):
                current_row = working_dataframe.loc[j].copy()

                # If the parent ion is within tolerance and the Lipid column is a string
                for n in [9]:
                    if within_tolerance(parent_ion, current_row[f'n-{n}'], tolerance) and isinstance(current_row['Lipid'], str):
                        working_dataframe.loc[i, 'Lipid'] = current_row['Lipid']
                        working_dataframe.loc[i, 'db_pos'] = f'n-{n}' + current_row['db_pos']

                        # Append to the final_dataframe
                        appended_row = working_dataframe.loc[i].copy()
                        appended_row['db_pos'] = f'n-{n}' + current_row['db_pos']
                        final_dataframe = final_dataframe.append(appended_row, ignore_index=True)
                for n in [7]:
                    if within_tolerance(parent_ion, current_row[f'n-{n}'], tolerance) and isinstance(current_row['Lipid'], str):
                        working_dataframe.loc[i, 'Lipid'] = current_row['Lipid']
                        working_dataframe.loc[i, 'db_pos'] = f'n-{n}' + current_row['db_pos']

                        # Append to the final_dataframe
                        appended_row = working_dataframe.loc[i].copy()
                        appended_row['db_pos'] = f'n-{n}' + current_row['db_pos']
                        final_dataframe = final_dataframe.append(appended_row, ignore_index=True)

    # Drop rows in the final_dataframe where 'Lipid' column value is NaN
    final_dataframe.dropna(subset=['Lipid'], inplace=True)

    return final_dataframe



def calculate_intensity_ratio(df):
    """
    Calculates the intensity ratio for each lipid in the DataFrame, based on their 'OzESI_Intensity'.

    :param df: DataFrame containing lipid information and intensity values.
    :return: Updated DataFrame with added 'Ratio' column.
    """
    # Create a new column for ratio
    df['Ratio'] = pd.Series(dtype='float64')

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        lipid = row['Lipid']
        label = row['db_pos']
        intensity = row['OzESI_Intensity']
        sample_id = row['Sample_ID']
        retention_time = row['Retention_Time']

        # Check if the label is n-9
        if label == 'n-9':
            # Check for the special case for TG(54:2)_FA18:1
            if lipid == 'TG(54:2)_FA18:1':
                # Apply retention time filter for the special case (around 20.00 ± 0.5)
                n7_row = df[
                    (df['Lipid'] == lipid) &
                    (df['db_pos'] == 'n-7') &
                    (df['Sample_ID'] == sample_id) &
                    (df['Retention_Time'].between(19.2, 21.5))
                ]
            else:
                # General case: Find the corresponding row with n-7 label and same lipid name and Sample_ID
                n7_row = df[
                    (df['Lipid'] == lipid) &
                    (df['db_pos'] == 'n-7') &
                    (df['Sample_ID'] == sample_id)
                ]

            # If a matching row is found, calculate the intensity ratio
            if not n7_row.empty:
                n7_intensity = n7_row['OzESI_Intensity'].values[0]
                ratio = intensity / n7_intensity

                # Assign the ratio to the 'Ratio' column
                df.at[index, 'Ratio'] = ratio

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



def filter_highest_ratio(df):
    """
    Filters the DataFrame to keep only rows with the highest ratio value for each unique Sample_ID and lipid.

    :param df: DataFrame containing lipid information and intensity ratio.
    :return: Filtered DataFrame.
    """
    # Sort the DataFrame by ratio in descending order
    df_sorted = df.sort_values(by='Ratio', ascending=False)

    # Drop duplicates keeping the first occurrence (highest ratio)
    df_filtered = df_sorted.drop_duplicates(subset=['Sample_ID', 'Lipid','db_pos'], keep='first')
    df_filtered = df_filtered.sort_values(by=['Sample_ID', 'Lipid'], ascending=[True, True])

    return df_filtered

