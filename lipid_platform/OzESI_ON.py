import pandas as pd
from tqdm import tqdm
import numpy as np
import re

# Define the function to process lipid data
def process_lipid_data(df):
    """
    Processes lipid data from a given DataFrame.

    Parameters:
    - df: pandas DataFrame containing the lipid data.

    Returns:
    - df_db_pos: pandas DataFrame with processed lipid data.
    """

    # Initialize df_DB_aldehyde
    df_DB_aldehyde = pd.DataFrame(columns=['DB_Position', 'Aldehyde_Ion'])
    for position in range(3, 21):
        df_DB_aldehyde.loc[position, 'DB_Position'] = position
        df_DB_aldehyde.loc[position, 'Aldehyde_Ion'] = 26 + (14 * (position - 3))

    # Copy df to df_db_pos and add n-# columns
    df_db_pos = df.copy()
    db_pos_list = [7, 9, 10, 12]
    for number in db_pos_list:
        df_db_pos[f'n-{number}'] = pd.NA

    # Calculate the n-# values
    for ozesi_position in db_pos_list:
        aldehyde_ion = df_DB_aldehyde.loc[df_DB_aldehyde["DB_Position"] == ozesi_position, "Aldehyde_Ion"].values[0]
        df_db_pos[f"n-{ozesi_position}"] = df_db_pos["Parent_Ion"] - aldehyde_ion

    # Determine 'FAC'
    def determine_fac(lipid):
        if 'FA18:1' in lipid:
            return 'FA18:1'
        elif 'FA16:1' in lipid:
            return 'FA16:1'
        else:
            return None
    df_db_pos['FAC'] = df_db_pos['Lipid'].apply(determine_fac)

    # Extract and sort 'TG' values
    def extract_tg(lipid):
        pattern = r'TG\(\d+:\d+\)'
        matches = re.findall(pattern, lipid)
        return ', '.join(matches) if matches else None

    def sort_tg_values(tg_cell):
        if not pd.isna(tg_cell) and ', ' in tg_cell:
            tg_list = tg_cell.split(', ')
            tg_list_sorted = sorted(tg_list, key=lambda x: [int(i) for i in x[3:-1].split(':')])
            return ', '.join(tg_list_sorted)
        else:
            return tg_cell
    
    def remove_tg_zero(lipid):
        pattern = r'TG\(\d+:0\)'
        new_lipid = re.sub(pattern, '', lipid)
        new_lipid = re.sub(r',\s+', ', ', new_lipid).strip(', ')
        return new_lipid

    df_db_pos['TG'] = df_db_pos['Lipid'].apply(extract_tg)
    df_db_pos['Lipid'] = df_db_pos['Lipid'].apply(remove_tg_zero)
    df_db_pos['TG'] = df_db_pos['TG'].apply(sort_tg_values)

    # Final DataFrame adjustments
    df_db_pos.reset_index(drop=True, inplace=True)
    df_db_pos.drop(columns=['Unnamed: 0', 'New_ID', 'Match_Group', 'Group_Sample'], inplace=True)

    return df_db_pos

# Note: The function calls are commented out to adhere to the instructions. This code is meant to be placed in a .py file for execution.


#### 
#OzESI match dbs
# Redefine the function after code execution state reset
def match_db_pos(df_db_pos, d1d, tolerance=0.3, retention_time_tolerance=0.275, columns_to_match=['Cage', 'Mouse', 'Genotype', 'Biology']):
    """
    Matches data between two DataFrames based on retention time, ion tolerance, and specific columns.

    Parameters:
    - df_db_pos: DataFrame to process.
    - d1d: DataFrame to match against.
    - tolerance: Tolerance for ion matching.
    - retention_time_tolerance: Tolerance for retention time matching.
    - columns_to_match: List of column names to match.

    Returns:
    - d2: DataFrame with matched rows.
    """
 # Assuming tqdm is available, if not, remove or replace with a simple loop

    # Re-initializing an empty DataFrame for filtered_d2
    d2 = pd.DataFrame()

    # Iterating through df_db_pos
    for index, row in tqdm(df_db_pos.iterrows(), total=df_db_pos.shape[0]):
        ground_truth_retention_time = row['Retention_Time']  # Ground truth retention time
        product_ion = row['Product_Ion']
        parent_ion_n7 = row['n-7']
        parent_ion_n9 = row['n-9']
        parent_ion_n10 = row['n-10']
        parent_ion_n12 = row['n-12']
        lipid_name = row['Lipid']  # Extracting Lipid name from df_db_pos

        # Defining the retention time window based on ground truth
        lower_bound_time = ground_truth_retention_time - retention_time_tolerance
        upper_bound_time = ground_truth_retention_time + retention_time_tolerance

        # Defining the tolerance for Parent and Product ions
        lower_bound_ion = product_ion - tolerance
        upper_bound_ion = product_ion + tolerance

        # Filtering d1d within the specified window, matching ions with tolerance, and matching additional columns
        matches = d1d[(d1d['Retention_Time'] >= lower_bound_time) & 
                         (d1d['Retention_Time'] <= upper_bound_time) &
                         (d1d['OzESI_Intensity'] >= 5) &
                         (d1d['Product_Ion'] >= lower_bound_ion) &
                         (d1d['Product_Ion'] <= upper_bound_ion) &
                         (d1d[columns_to_match].eq(row[columns_to_match])).all(axis=1)]

        # Adding 'db' and 'Lipid' columns based on the tolerance matching
        matches['db'] = ''
        matches['Lipid'] = lipid_name  # Adding the Lipid name to all matches
        matches.loc[(matches['Parent_Ion'] >= parent_ion_n7 - tolerance) & 
                    (matches['Parent_Ion'] <= parent_ion_n7 + tolerance), 'db'] = 'n-7'
        matches.loc[(matches['Parent_Ion'] >= parent_ion_n9 - tolerance) & 
                    (matches['Parent_Ion'] <= parent_ion_n9 + tolerance), 'db'] = 'n-9'
        matches.loc[(matches['Parent_Ion'] >= parent_ion_n10 - tolerance) & 
                    (matches['Parent_Ion'] <= parent_ion_n10 + tolerance), 'db'] = 'n-10'
        matches.loc[(matches['Parent_Ion'] >= parent_ion_n12 - tolerance) & 
                    (matches['Parent_Ion'] <= parent_ion_n12 + tolerance), 'db'] = 'n-12'

        # Append the matching rows to d2
        d2 = d2.append(matches)

    return d2

# Note: The function calls are commented out to adhere to the instructions. This code is meant to be placed in a .py file for execution.
