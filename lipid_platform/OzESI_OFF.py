#Import all the necessary python libraries
import pymzml
import csv
import os
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import json
from scipy.integrate import trapz
from tqdm import tqdm

#Import all the necessary CLAW libraries
import create_directory
import CLAW
import matplotlib.pyplot as plt
import warnings

import re
from sklearn.mixture import GaussianMixture
from tqdm import tqdm

# Suppress all warnings
warnings.simplefilter(action='ignore', category=FutureWarning)




def RT_filter(df: pd.DataFrame, retention_time_range: tuple) -> pd.DataFrame:
    """
    Filters the DataFrame based on a specified retention time range.

    Args:
        df (pd.DataFrame): The DataFrame to be filtered.
        retention_time_range (tuple): A tuple containing the lower and upper bounds of the retention time range.

    Returns:
        pd.DataFrame: The filtered DataFrame containing only rows where the retention time is within the specified range.
    """
    # Ensure column 'Retention_Time' exists in the DataFrame
    if 'Retention_Time' not in df.columns:
        raise ValueError("DataFrame does not contain the column 'Retention_Time'.")

    # Filter the DataFrame to keep only rows where Retention_Time is within the specified range
    filtered_df = df[(df['Retention_Time'] >= retention_time_range[0]) & (df['Retention_Time'] <= retention_time_range[1])]

    return filtered_df



def create_match_group(df: pd.DataFrame, group_columns: list) -> pd.DataFrame:
    """
    Adds a 'Match_Group' column to the DataFrame based on grouping by specified columns.
    
    Args:
        df (pd.DataFrame): The DataFrame to modify.
        group_columns (list): A list of column names to group by.
        
    Returns:
        pd.DataFrame: The modified DataFrame with an added 'Match_Group' column.
    """
    # Check if the specified columns exist in the DataFrame
    missing_columns = [col for col in group_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame: {missing_columns}")

    # Group by the specified columns and assign group numbers
    df['Match_Group'] = df.groupby(group_columns).ngroup()

    return df


def match_lipid(d1a: pd.DataFrame, df_MRM: pd.DataFrame, tolerance: float = 0.3) -> pd.DataFrame:
    """
    Matches lipids from df_MRM to d1a based on Parent_Ion and Product_Ion within a specified tolerance.
    Adds a 'Lipid' column to d1a indicating the matched lipid. Includes a progress bar for tracking execution time.

    Args:
        d1a (pd.DataFrame): The DataFrame containing ions to match.
        df_MRM (pd.DataFrame): The DataFrame containing the reference lipids and their ions.
        tolerance (float, optional): The tolerance within which two ions are considered a match. Defaults to 0.3.

    Returns:
        pd.DataFrame: A copy of d1a with an additional 'Lipid' column indicating the matched lipid for each ion pair.
    """
    # Copy d1a to avoid modifying the original DataFrame
    d1b = d1a.copy()

    # Create Match_Group in d1a if it doesn't exist
    if 'Match_Group' not in d1a.columns:
        d1a['Match_Group'] = d1a.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()

    # Copy Match_Group to d1b
    d1b['Match_Group'] = d1a['Match_Group']

    # Helper function to check if two ions are within the tolerance
    def is_within_tolerance(ion1, ion2, tolerance=0.3):
        return abs(ion1 - ion2) <= tolerance

    # Add a new column for Lipid in d1b
    d1b['Lipid'] = None

    # Iterate through each unique Match_Group in d1a with a progress bar
    for group in tqdm(d1a['Match_Group'].unique(), desc="Matching Lipids"):
        # Extract a representative row for the current group from d1a
        group_row = d1a[d1a['Match_Group'] == group].iloc[0]

        # Find a matching lipid in df_MRM for the representative row
        for _, mrm_row in df_MRM.iterrows():
            if is_within_tolerance(group_row['Parent_Ion'], mrm_row['Parent_Ion'], tolerance) and is_within_tolerance(group_row['Product_Ion'], mrm_row['Product_Ion'], tolerance):
                # Assign the lipid to all rows in the corresponding group in d1b
                d1b.loc[d1b['Match_Group'] == group, 'Lipid'] = mrm_row['Lipid']
                break  # Stop searching once a match is found

    return d1b


import pandas as pd
from tqdm.auto import tqdm

def assign_correct_rt(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assigns the correct retention time (Correct_RT) for each group in a DataFrame based on the maximum OzESI_Intensity.

    Args:
        df (pd.DataFrame): The DataFrame to process, which must contain 'Match_Group', 'OzESI_Intensity', and 'Retention_Time' columns.

    Returns:
        pd.DataFrame: The DataFrame with an additional 'Correct_RT' column indicating the correct retention time for each group.
    """
    tqdm.pandas(desc="Assigning Correct RT")

    # Ensure required columns are present
    required_columns = {'Match_Group', 'OzESI_Intensity', 'Retention_Time'}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns in DataFrame: {missing_columns}")

    # Group by Match_Group and find the Retention_Time corresponding to the max OzESI_Intensity for each group
    max_rt_per_group = df.groupby('Match_Group').progress_apply(
        lambda x: x.loc[x['OzESI_Intensity'].idxmax(), 'Retention_Time']
    )

    # Map the max retention time to the Correct_RT column for each group
    df['Correct_RT'] = df['Match_Group'].map(max_rt_per_group)

    return df


import pandas as pd
from tqdm.auto import tqdm

def extract_and_process_details(df: pd.DataFrame, intensity_threshold: int = 2000) -> pd.DataFrame:
    """
    Processes the DataFrame by filtering groups, mapping max intensities, rounding values,
    dropping duplicates, extracting details from Sample_ID, and filtering by intensity threshold.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        intensity_threshold (int): The minimum max intensity for a Match_Group to be included.

    Returns:
        pd.DataFrame: The processed DataFrame.
    """
    # Filter out groups where all 'Lipid' values are NaN
    df_filtered = df.groupby('Match_Group').filter(lambda x: not x['Lipid'].isna().all())

    # Find the max OzESI_Intensity for each Match_Group
    max_intensity_per_group = df_filtered.groupby('Match_Group')['OzESI_Intensity'].max()

    # Map the max intensity to a new column Max_Intensity for each group
    df_filtered['Max_Intensity'] = df_filtered['Match_Group'].map(max_intensity_per_group)

    # Round specified columns
    df_filtered['Max_Intensity'] = df_filtered['Max_Intensity'].round(0)
    df_filtered['OzESI_Intensity'] = df_filtered['OzESI_Intensity'].round(0)
    df_filtered['Retention_Time'] = df_filtered['Retention_Time'].round(2)

    # # Keep only the row with the highest OzESI_Intensity in each Match_Group
    # df_filtered = df_filtered.sort_values('OzESI_Intensity', ascending=False).drop_duplicates('Match_Group')

    # Extract details from Sample_ID
    df_filtered = extract_details_from_sample_id(df_filtered)

    # Filter out groups where the max Max_Intensity is under the threshold
    df_final = df_filtered.groupby('Match_Group').filter(lambda x: x['Max_Intensity'].max() >= intensity_threshold)

    return df_final

def extract_details_from_sample_id(df: pd.DataFrame, column_name: str = 'Sample_ID') -> pd.DataFrame:
    """
    Extracts details from the Sample_ID column and adds them as new columns: Cage, Mouse, Genotype, and Biology.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        column_name (str): The column from which to extract details.

    Returns:
        pd.DataFrame: The DataFrame with additional columns extracted from Sample_ID.
    """
    tqdm.pandas(desc="Extracting Sample Details")
    pattern = r'^[^_]*_(?P<Cage>[^_]+)_(?P<Mouse>[^_]+)_(?P<Genotype>[^_]+)_(?P<Biology>[^_]+)'
    df_extracted = df[column_name].progress_apply(lambda x: pd.Series(re.match(pattern, x).groups(), index=['Cage', 'Mouse', 'Genotype', 'Biology']) if re.match(pattern, x) else pd.Series([None, None, None, None], index=['Cage', 'Mouse', 'Genotype', 'Biology']))
    df = pd.concat([df, df_extracted], axis=1)

    return df




def add_group_sample_column(df: pd.DataFrame, group_columns: list) -> pd.DataFrame:
    """
    Adds a new column 'Group_Sample' to the DataFrame, assigning a unique group number 
    for each unique combination of specified columns.

    Args:
        df (pandas.DataFrame): The DataFrame to process.
        group_columns (list): A list of column names to group by for assigning unique group numbers.

    Returns:
        pandas.DataFrame: The DataFrame with the added 'Group_Sample' column.
    """
    # Validate if all specified columns exist in the DataFrame
    missing_columns = [col for col in group_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing columns in DataFrame: {missing_columns}")

    # Create the 'Group_Sample' column by assigning a group number for each unique combination of specified columns
    df['Group_Sample'] = df.groupby(group_columns).ngroup()

    return df


def add_group_sample_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a new column 'Group_Sample' to the DataFrame, assigning a unique group number 
    for each combination of Cage, Mouse, Genotype, Biology, and Lipid.

    Args:
        df (pd.DataFrame): The DataFrame to process.

    Returns:
        pd.DataFrame: The DataFrame with the added 'Group_Sample' column.
    """
    # Create the 'Group_Sample' column by assigning a group number for each combination
    df['Group_Sample'] = df.groupby(['Cage', 'Mouse', 'Genotype', 'Biology', 'Lipid']).ngroup()
    return df

### OzON correct RTs


def count_exact_matches_and_create_id(df1, df2, columns):
    """
    Counts the number of exact matches for specified columns between two DataFrames
    and creates a new column 'New_ID' in each DataFrame for matched rows.

    Args:
    df1 (pandas.DataFrame): The first DataFrame.
    df2 (pandas.DataFrame): The second DataFrame.
    columns (list): A list of column names to compare.

    Returns:
    int: The number of exact matches.
    """
    # Create a temporary 'ID' column by concatenating the specified columns
    df1['ID'] = df1[columns].dropna().apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    df2['ID'] = df2[columns].dropna().apply(lambda row: '-'.join(row.values.astype(str)), axis=1)

    # Create sets of these IDs
    set_df1 = set(df1['ID'])
    set_df2 = set(df2['ID'])

    # Find the intersection of these sets
    matches = set_df1.intersection(set_df2)

    # Create 'New_ID' column based on whether 'ID' is in the matches
    df1['New_ID'] = df1['ID'].apply(lambda x: x if x in matches else None)
    df2['New_ID'] = df2['ID'].apply(lambda x: x if x in matches else None)

    # Drop the temporary 'ID' column
    df1.drop(columns=['ID'], inplace=True)
    df2.drop(columns=['ID'], inplace=True)

    # Return the number of matches
    return len(matches)


def Correct_OzON_RT(df):
    """
    Assigns the Retention_Time of the row with the highest OzESI_Intensity
    for each group of Biology, Cage, Genotype, and Lipid as Correct_RT.

    Args:
        df (pandas.DataFrame): The DataFrame to process.

    Returns:
        pandas.DataFrame: The DataFrame with the new Correct_RT column.
    """
    # Group by the specified columns and find the max intensity for each group
    max_intensity_per_group = df.groupby(['Biology', 'Cage', 'Genotype', 'Lipid'])['OzESI_Intensity'].transform('max')

    # Mark the rows with the highest intensity in each group
    is_max_intensity = df['OzESI_Intensity'] == max_intensity_per_group

    # Copy the Retention_Time of the row with the max intensity to Correct_RT
    df['Correct_RT'] = df[is_max_intensity]['Retention_Time'].round(2)

    # Forward fill and backward fill Correct_RT within each group
    df['Correct_RT'] = df.groupby(['Biology', 'Cage', 'Genotype', 'Lipid'])['Correct_RT'].apply(lambda x: x.ffill().bfill())

    return df