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

# Add this at the beginning to enable more detailed debugging
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set pandas to show more detailed warnings
pd.set_option('mode.chained_assignment', 'warn')

def create_analysis_dataframes():
    """
    Creates and returns three DataFrames for storing time-intensity data, master data, and OzESI time data.

    Returns:
        pd.DataFrame: DataFrame for storing time and intensity values.
        pd.DataFrame: Master DataFrame for storing parent ion, product ion, intensity, transition, and sample ID.
        pd.DataFrame: DataFrame for storing OzESI (OzID Electron Spray Ionization) data including parent ion, product ion, retention time, intensity, sample ID, and transition.
    """
    logger.debug("Creating analysis dataframes")
    time_intensity_dataframe = pd.DataFrame(columns=['Time', 'Intensity'])
    master_lipid_dataframe = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
    OzESI_time_dataframe = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])
    
    logger.debug(f"Created dataframes - shapes: time_intensity={time_intensity_dataframe.shape}, master_lipid={master_lipid_dataframe.shape}, OzESI_time={OzESI_time_dataframe.shape}")
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
    logger.debug("Starting pre-parsing setup")

    # Check and create folders if they do not exist
    logger.debug(f"Creating directories if they don't exist: {data_base_name_location}, {Project}, {Project_Folder_data}, {Project_results}")
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
    
    logger.debug("Pre-parsing setup completed")
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
    logger.debug(f"Reading MRM list from {filename}")
    
    try:
        # Concatenate all sheets from the Excel file into one DataFrame
        raw_mrm_data = pd.read_excel(filename, sheet_name=None)
        logger.debug(f"Read Excel file with sheets: {list(raw_mrm_data.keys())}")
        
        concatenated_mrm_data = pd.concat(raw_mrm_data, ignore_index=True)
        logger.debug(f"Concatenated data shape: {concatenated_mrm_data.shape}")

        # Extract the required columns and create a proper copy
        lipid_MRM_data = concatenated_mrm_data[['Compound Name', 'Parent Ion', 'Product Ion', 'Class']].copy()
        logger.debug(f"Extracted columns shape: {lipid_MRM_data.shape}")
        
        # Rename columns
        lipid_MRM_data.columns = lipid_MRM_data.columns.str.replace(' ', '_')
        
        # Use .loc to modify columns
        logger.debug("Rounding and formatting data")
        lipid_MRM_data.loc[:, 'Parent_Ion'] = np.round(lipid_MRM_data['Parent_Ion'], 1)
        lipid_MRM_data.loc[:, 'Product_Ion'] = np.round(lipid_MRM_data['Product_Ion'], 1)
        lipid_MRM_data.loc[:, 'Transition'] = lipid_MRM_data['Parent_Ion'].astype(str) + ' -> ' + lipid_MRM_data['Product_Ion'].astype(str)
        lipid_MRM_data = lipid_MRM_data.rename(columns={'Compound_Name': 'Lipid'})

        # Optionally filter the data to keep only specific lipid classes
        if remove_std:
            lipid_classes_to_keep = ['PS', 'PG', 'CE', 'PC', 'DAG', 'PE', 'TAG', 'FA', 'Cer', 'CAR', 'PI', 'SM']
            logger.debug(f"Filtering to keep only these classes: {lipid_classes_to_keep}")
            initial_size = len(lipid_MRM_data)
            lipid_MRM_data = lipid_MRM_data[lipid_MRM_data['Class'].isin(lipid_classes_to_keep)]
            logger.debug(f"Filtered data size: {len(lipid_MRM_data)} (removed {initial_size - len(lipid_MRM_data)} rows)")

        # Optionally adjust the ion values for deuterated lipids
        if deuterated:
            logger.debug("Adjusting for deuterated lipids")
            lipid_MRM_data.loc[:, 'Parent_Ion'] += 1
            lipid_MRM_data.loc[:, 'Product_Ion'] += 1
            # Update the Transition column with the updated values
            lipid_MRM_data.loc[:, 'Transition'] = lipid_MRM_data['Parent_Ion'].astype(str) + ' -> ' + lipid_MRM_data['Product_Ion'].astype(str)
        
        logger.debug(f"MRM list processing complete. Final shape: {lipid_MRM_data.shape}")
        return lipid_MRM_data
    
    except Exception as e:
        logger.error(f"Error in read_mrm_list: {str(e)}")
        raise


def create_ion_dict(mrm_database):
    """
    Creates a dictionary of ions from an MRM database DataFrame.
    
    :param mrm_database: DataFrame containing MRM database information.
    
    :return: A dictionary with ion pairs as keys, and a list of tuples containing corresponding lipid and class as values.
    """
    logger.debug("Creating ion dictionary")
    ion_dict = defaultdict(list)
    try:
        for index, row in mrm_database.iterrows():
            ion_dict[(row['Parent_Ion'], row['Product_Ion'])].append((row['Lipid'], row['Class']))
        
        logger.debug(f"Ion dictionary created with {len(ion_dict)} unique ion pairs")
        return ion_dict
    except Exception as e:
        logger.error(f"Error in create_ion_dict: {str(e)}")
        raise


# Declare the DataFrame globally if it's used across multiple functions
time_and_intensity_df = pd.DataFrame(columns=['Time', 'Intensity'])
master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
OzESI_time_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])
logger.debug("Global DataFrames initialized")

def mzml_parser(file_path, plot_chromatogram=False):
    global master_df
    global OzESI_time_df
    global time_and_intensity_df
    
    logger.debug(f"Starting mzML parsing for file: {file_path}")
    
    rows = []
    ozesi_rows = []
    
    try:
        run = pymzml.run.Reader(file_path, skip_chromatogram=False)
        logger.debug(f"Successfully opened mzML file: {file_path}")
        q1_mz = 0
        q3_mz = 0

        spectrum_count = 0
        for spectrum in run:
            spectrum_count += 1
            if spectrum_count % 100 == 0:
                logger.debug(f"Processed {spectrum_count} spectra")
                
            for element in spectrum.ID.split(' '):
                if 'Q1' in element:
                    q1 = element.split('=')
                    q1_mz = np.round(float(q1[1]), 1)

                if 'Q3' in element:
                    q3 = element.split('=')
                    q3_mz = np.round(float(q3[1]), 1)

                    ##############
                    # Plotting chromatogram if the condition is met and plot_chromatogram is True
                    if plot_chromatogram and within_tolerance(q1_mz, 876.6) and within_tolerance(q3_mz, 577.6):
                        logger.debug(f"Plotting chromatogram for transition 876.6 -> 577.6")
                        times, intensities = zip(*spectrum.peaks())
                        plt.plot(times, intensities)
                        plt.xlabel('Time')
                        plt.ylabel('Intensity')
                        plt.title('Chromatogram for 876.6 -> 577.6')
                        plt.show()
                    ###########

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
        
        logger.debug(f"Collected {len(rows)} rows and {len(ozesi_rows)} OzESI rows")
        df = pd.DataFrame(rows)
        
        logger.debug(f"Before concat - master_df shape: {master_df.shape}, OzESI_time_df shape: {OzESI_time_df.shape}")
        # Use concat instead of append
        if len(ozesi_rows) > 0:
            OzESI_time_df = pd.concat([OzESI_time_df, pd.DataFrame(ozesi_rows)], ignore_index=True)
        if len(rows) > 0:
            master_df = pd.concat([master_df, df], ignore_index=True)
        
        logger.debug(f"After concat - master_df shape: {master_df.shape}, OzESI_time_df shape: {OzESI_time_df.shape}")
        print(f'Finished parsing mzML file: {file_path}\n')
        
    except Exception as e:
        logger.error(f"Error in mzml_parser for file {file_path}: {str(e)}")
        raise


def mzml_parser_batch(folder_name, plot_chromatogram=False):
    global master_df
    global time_and_intensity_df
    
    logger.debug(f"Starting batch parsing from folder: {folder_name}")
    
    try:
        data_folder = os.listdir(folder_name)
        data_folder.sort()
        logger.debug(f"Found {len(data_folder)} files in folder")

        mzml_files = [file for file in data_folder if file.endswith('.mzML')]
        logger.debug(f"Found {len(mzml_files)} mzML files to process")

        for i, file in enumerate(mzml_files):
            logger.debug(f"Processing file {i+1}/{len(mzml_files)}: {file}")
            file_path = os.path.join(folder_name, file)
            mzml_parser(file_path, plot_chromatogram=plot_chromatogram)
        
        logger.debug(f"Batch parsing complete. master_df shape: {master_df.shape}")
        print('Finished parsing all mzML files\n')
    
    except Exception as e:
        logger.error(f"Error in mzml_parser_batch: {str(e)}")
        raise


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

    try:
        for key, value in ion_dict.items():
            if within_tolerance(ions[0], key[0], tolerance) and within_tolerance(ions[1], key[1], tolerance):
                matched_lipids.extend([match[0] for match in value])
                matched_classes.extend([match[1] for match in value])

        if matched_lipids and matched_classes:
            row['Lipid'] = ' | '.join(matched_lipids)
            row['Class'] = ' | '.join(matched_classes)

        return row
    
    except Exception as e:
        logger.error(f"Error in match_ions for ions {ions}: {str(e)}")
        return row


def match_lipids_parser(mrm_database, df, tolerance=0.3):
    """
    Performs lipid matching by creating an ion dictionary from the MRM database and applying the match_ions function to each row of a DataFrame.
    
    :param mrm_database: DataFrame containing MRM database information.
    :param df: DataFrame containing ion information to be matched.
    :param tolerance: The acceptable difference between ion values to be considered a match. Defaults to 0.3.
    
    :return: DataFrame with matched lipid and class information if matches were found.
    """
    logger.debug(f"Starting lipid matching with tolerance={tolerance}")
    logger.debug(f"mrm_database shape: {mrm_database.shape}, df shape: {df.shape}")
    
    try:
        # Make sure the necessary columns exist in df
        if not all(col in df.columns for col in ['Parent_Ion', 'Product_Ion']):
            missing_cols = [col for col in ['Parent_Ion', 'Product_Ion'] if col not in df.columns]
            logger.error(f"Missing required columns in df: {missing_cols}")
            raise ValueError(f"Missing required columns in df: {missing_cols}")
            
        ion_dict = create_ion_dict(mrm_database)
        
        # Add Lipid and Class columns to df if they don't exist
        if 'Lipid' not in df.columns:
            df['Lipid'] = None
        if 'Class' not in df.columns:
            df['Class'] = None
        
        logger.debug("Applying match_ions to each row")
        df_matched = df.apply(lambda row: match_ions(row, ion_dict=ion_dict, tolerance=tolerance), axis=1)
        
        # Count matches
        match_count = df_matched['Lipid'].count()
        logger.debug(f"Matching complete. Found matches for {match_count} out of {len(df_matched)} rows")
        
        return df_matched
    
    except Exception as e:
        logger.error(f"Error in match_lipids_parser: {str(e)}")
        raise


def save_dataframe(df, Project_results, file_name_to_save, max_attempts=5):
    """
    Saves a given DataFrame to a CSV file within a specified directory.
    
    :param df: DataFrame to be saved.
    :param Project_results: The project directory to save results in.
    :param file_name_to_save: The desired filename for the saved DataFrame.
    :param max_attempts: The maximum number of attempts to save the DataFrame. Defaults to 5.
    
    :return: None
    """
    logger.debug(f"Attempting to save DataFrame to {Project_results}/{file_name_to_save}")
    
    folder_path = f'data_results/data/data_matching/{Project_results}'
    try:
        os.makedirs(folder_path, exist_ok=True)
        logger.debug(f"Created folder: {folder_path}")

        for i in range(max_attempts):
            file_path = f'{folder_path}/{file_name_to_save}.csv'
            logger.debug(f"Attempt {i+1}: Trying to save to {file_path}")
            
            if not os.path.isfile(file_path):
                df.to_csv(file_path, index=False)
                logger.debug(f"Successfully saved DataFrame to {file_path}")
                print(f"Saved DataFrame to {file_path}")
                break
            else:
                logger.debug(f"File already exists: {file_path}")
        else:
            logger.error(f"Failed to save DataFrame after {max_attempts} attempts.")
            print(f"Failed to save DataFrame after {max_attempts} attempts.")
            return None
            
    except Exception as e:
        logger.error(f"Error in save_dataframe: {str(e)}")
        raise


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
    global OzESI_time_df
    
    logger.debug(f"Starting full parse with parameters: data_base={data_base_name_location}, folder={Project_Folder_data}, " + 
                 f"results={Project_results}, filename={file_name_to_save}, tolerance={tolerance}, " + 
                 f"remove_std={remove_std}, save_data={save_data}, batch_processing={batch_processing}, " +
                 f"plot_chromatogram={plot_chromatogram}")
    
    try:
        # Reset global DataFrames to ensure clean state
        master_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Intensity', 'Transition', 'Sample_ID'])
        OzESI_time_df = pd.DataFrame(columns=['Parent_Ion', 'Product_Ion', 'Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition'])
        logger.debug("Reset global DataFrames")
        
        # Read MRM database
        mrm_database = read_mrm_list(data_base_name_location, remove_std=remove_std)
        logger.debug(f"MRM database read with shape: {mrm_database.shape}")
        
        # Batch processing
        if batch_processing:
            logger.debug(f"Starting batch processing from folder: {Project_Folder_data}")
            mzml_parser_batch(Project_Folder_data, plot_chromatogram=plot_chromatogram)  
        # Single file processing
        else:
            logger.debug(f"Processing single file: {Project_Folder_data}")
            mzml_parser(Project_Folder_data, plot_chromatogram=plot_chromatogram) 
        
        logger.debug(f"Parsing complete. master_df shape: {master_df.shape}, OzESI_time_df shape: {OzESI_time_df.shape}")
        
        if master_df.empty:
            logger.error("No data was parsed from the mzML files. master_df is empty.")
            raise ValueError("No data was parsed from the mzML files")
            
        # Match lipids
        logger.debug("Starting lipid matching")
        df_matched = match_lipids_parser(mrm_database, master_df, tolerance=tolerance)
        logger.debug(f"Matching complete. df_matched shape: {df_matched.shape}")
        
        # Save data if requested
        if save_data:
            logger.debug(f"Saving matched data to {file_name_to_save}")
            save_dataframe(df_matched, Project_results, file_name_to_save)

        logger.debug("Full parse complete")
        return df_matched, OzESI_time_df
    
    except Exception as e:
        logger.error(f"Error in full_parse: {str(e)}")
        raise


def filter_rt(df, min_rt=10.0, max_rt=20.0, min_intensity=None):
    """
    Filters the DataFrame based on retention times and aggregates by max intensity for unique 'Sample_ID' and 'Transition' combinations.
    
    Parameters:
        df (pd.DataFrame): Input DataFrame with columns 'Retention_Time' and 'OzESI_Intensity'.
        min_rt (float, optional): Minimum retention time for filtering. Defaults to 10.0.
        max_rt (float, optional): Maximum retention time for filtering. Defaults to 20.0.
        min_intensity (float, optional): Minimum intensity for filtering. If None, no filtering by intensity is done.
        
    Returns:
        pd.DataFrame: Filtered and aggregated DataFrame.
    """
    logger.debug(f"Filtering by retention time: min_rt={min_rt}, max_rt={max_rt}, min_intensity={min_intensity}")
    logger.debug(f"Input df shape: {df.shape}")
    
    try:
        # Check if the required columns exist
        required_columns = ['Retention_Time', 'OzESI_Intensity', 'Sample_ID', 'Transition']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"Missing required columns in df: {missing_columns}")
            logger.debug(f"Available columns: {df.columns.tolist()}")
            raise ValueError(f"Missing required columns in df: {missing_columns}")
        
        # Filter based on retention time
        filtered_df = df[(df['Retention_Time'] >= min_rt) & (df['Retention_Time'] <= max_rt)].copy()
        logger.debug(f"After RT filtering: {len(filtered_df)} rows (from {len(df)})")

        # Filter based on intensity if min_intensity is provided
        if min_intensity is not None:
            filtered_df = filtered_df[filtered_df['OzESI_Intensity'] >= min_intensity]
            logger.debug(f"After intensity filtering: {len(filtered_df)} rows")

        # Round the values
        filtered_df.loc[:, 'Retention_Time'] = filtered_df['Retention_Time'].round(2)
        filtered_df.loc[:, 'OzESI_Intensity'] = filtered_df['OzESI_Intensity'].round(0)
        logger.debug("Rounded Retention_Time and OzESI_Intensity values")

        # Check for empty dataframe before aggregation
        if filtered_df.empty:
            logger.warning("No data remained after filtering")
            return filtered_df

        # Check unique combinations
        unique_combinations = filtered_df.groupby(['Sample_ID', 'Transition']).size().reset_index(name='count')
        logger.debug(f"Found {len(unique_combinations)} unique Sample_ID/Transition combinations")

        # Aggregate by max intensity for unique combinations of 'Sample_ID' and 'Transition'
        logger.debug("Aggregating by max intensity")
        filtered_df = filtered_df.groupby(['Sample_ID', 'Transition']).apply(
            lambda x: x.loc[x['OzESI_Intensity'].idxmax()]).reset_index(drop=True)
        logger.debug(f"After aggregation: {len(filtered_df)} rows")

        return filtered_df
    
    except Exception as e:
        logger.error(f"Error in filter_rt: {str(e)}")
        raise


def concat_dataframes(df_matched, filtered_df):
    """
    Concatenates two DataFrames along the columns.
    
    Parameters:
        df_matched (pd.DataFrame): First DataFrame.
        filtered_df (pd.DataFrame): Second DataFrame, only the 'Retention_Time' and 'OzESI_Intensity' columns will be used.
        
    Returns:
        pd.DataFrame: Concatenated DataFrame.
    """
    logger.debug(f"Concatenating DataFrames: df_matched shape={df_matched.shape}, filtered_df shape={filtered_df.shape}")
    
    try:
        # Check if the required columns exist in filtered_df
        required_columns = ['Retention_Time', 'OzESI_Intensity']
        missing_columns = [col for col in required_columns if col not in filtered_df.columns]
        if missing_columns:
            logger.error(f"Missing required columns in filtered_df: {missing_columns}")
            logger.debug(f"Available columns in filtered_df: {filtered_df.columns.tolist()}")
            raise ValueError(f"Missing required columns in filtered_df: {missing_columns}")
            
        # Check for matching indexes or common columns to merge on
        logger.debug(f"df_matched columns: {df_matched.columns.tolist()}")
        logger.debug(f"filtered_df columns: {filtered_df.columns.tolist()}")
        
        # Try to find common columns to merge on
        common_columns = [col for col in df_matched.columns if col in filtered_df.columns]
        logger.debug(f"Common columns for potential merge: {common_columns}")
        
        # Perform the concatenation
        result = pd.concat([df_matched, filtered_df[['Retention_Time', 'OzESI_Intensity']]], axis=1)
        logger.debug(f"Concatenation result shape: {result.shape}")
        
        return result
    
    except Exception as e:
        logger.error(f"Error in concat_dataframes: {str(e)}")
        raise

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
    logger.debug(f"Calculating DB positions for positions: {db_pos_list}")
    logger.debug(f"Input df shape: {df_matched_ions.shape}")
    
    try:
        # Check if Parent_Ion column exists
        if 'Parent_Ion' not in df_matched_ions.columns:
            logger.error("Missing required column 'Parent_Ion' in df_matched_ions")
            logger.debug(f"Available columns: {df_matched_ions.columns.tolist()}")
            raise ValueError("Missing required column 'Parent_Ion' in df_matched_ions")
            
        # Create a DataFrame to store the DB_Position and corresponding Aldehyde_Ion values
        df_DB_aldehyde = pd.DataFrame(columns=['DB_Position','Aldehyde_Ion'])
        logger.debug("Created DB_aldehyde DataFrame")

        # Loop through the range of DB_Position values to calculate the corresponding Aldehyde_Ion values
        for position in range(3, 21):
            df_DB_aldehyde.loc[position, 'DB_Position'] = position
            df_DB_aldehyde.loc[position, 'Aldehyde_Ion'] = 26 + (14 * (position-3))
        
        logger.debug(f"Filled DB_aldehyde DataFrame with positions 3-20: {df_DB_aldehyde.head()}")

        # Loop through the specified db_pos_list
        for ozesi_position in db_pos_list:
            logger.debug(f"Processing DB position: {ozesi_position}")
            
            # Retrieve the corresponding Aldehyde_Ion value for the current DB_Position
            aldehyde_ion_row = df_DB_aldehyde.loc[df_DB_aldehyde["DB_Position"] == ozesi_position]
            
            if aldehyde_ion_row.empty:
                logger.warning(f"No Aldehyde_Ion found for DB_Position {ozesi_position}")
                continue
                
            aldehyde_ion = aldehyde_ion_row["Aldehyde_Ion"].values[0]
            logger.debug(f"Aldehyde_Ion for position {ozesi_position}: {aldehyde_ion}")

            # Calculate and store the n-i value for the current OzESI position
            column_name = f"n-{ozesi_position}"
            df_matched_ions[column_name] = df_matched_ions["Parent_Ion"] - aldehyde_ion
            logger.debug(f"Added column {column_name} to df_matched_ions")

        logger.debug(f"DB position calculation complete. Output df shape: {df_matched_ions.shape}")
        return df_matched_ions
    
    except Exception as e:
        logger.error(f"Error in calculate_DB_Position: {str(e)}")
        raise