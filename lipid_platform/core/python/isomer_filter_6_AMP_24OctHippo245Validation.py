import pandas as pd
import os

def extract_sample_name(file_name):
    """
    Extracts the sample name based on the first four underscore-separated parts in the file name.
    Ignores any part after the fourth underscore.
    
    Parameters:
    file_name (str): The full file name.
    
    Returns:
    str: The extracted sample name (first four parts joined by underscores).
    """
    return '_'.join(file_name.split('_')[2:6])

def filter_ozON_by_ozOFF(input_dir, off_possible_dir, retention_time_tolerance=0.3, isomer_filter_output='Projects/AMP/isomer_filter_6/'):
    """
    Iterates over all files in the input directory, filtering each df_analysis (OzON data) based on Retention_Time, 
    Species, and Isomer from the corresponding off_possible file. Matches files based on the sample name 
    extracted from the file names. Keeps the whole Sample value for the output file name and adds `_isomer_filtered_6`.
    
    Parameters:
    input_dir (str): Directory containing parquet files to process.
    off_possible_dir (str): Directory containing parquet files for off_possible data.
    retention_time_tolerance (float): The tolerance for filtering retention times (default: 0.3).
    isomer_filter_output (str): Path where the filtered parquet files will be saved.

    Returns:
    None: The function saves the filtered DataFrame as parquet files in the specified output path.
    """
    # Create the output directory if it does not exist
    os.makedirs(isomer_filter_output, exist_ok=True)

    # Get a list of all parquet files in the off_possible directory
    off_possible_files = {extract_sample_name(f): f for f in os.listdir(off_possible_dir) if f.endswith(".parquet")}

    # Iterate over all parquet files in the input directory (OzON data)
    for parquet_file in os.listdir(input_dir):
        if parquet_file.endswith(".parquet"):
            # Read the df_analysis file
            file_path = os.path.join(input_dir, parquet_file)
            df_analysis = pd.read_parquet(file_path)

            # Extract the sample name from the df_analysis file
            sample_name = extract_sample_name(parquet_file)
            print(f"Processing df_analysis file: {parquet_file}, Sample Name: {sample_name}")

            # Look for the corresponding off_possible file
            if sample_name in off_possible_files:
                off_possible_file = off_possible_files[sample_name]
                off_possible_path = os.path.join(off_possible_dir, off_possible_file)
                print(f"Matching off_possible file found: {off_possible_file}")

                # Load the off_possible DataFrame
                off_possible = pd.read_parquet(off_possible_path)

                # Create a new column Adjusted_RT based on Retention_Time and STD_RT_Dif
                if 'STD_RT_Dif' in df_analysis.columns:
                    df_analysis['Adjusted_RT'] = df_analysis['Retention_Time'] + df_analysis['STD_RT_Dif']
                else:
                    df_analysis['Adjusted_RT'] = df_analysis['Retention_Time']

                # Create the Isomer column in df_analysis if it does not exist
                if 'Isomer' not in df_analysis.columns:
                    df_analysis['Isomer'] = df_analysis.groupby('Lipid')['Retention_Time'].rank(method='dense').astype(int)

                # Initialize an empty DataFrame to store the filtered results
                filtered_df = pd.DataFrame()

                # Iterate over each row in the off_possible (OzOFF data)
                for _, row in off_possible.iterrows():
                    species = row['Species']
                    isomer = row['Isomer']
                    retention_time_off = row['Retention_Time']

                    # Define the retention time window based on the tolerance
                    retention_time_start = retention_time_off - retention_time_tolerance
                    retention_time_end = retention_time_off + retention_time_tolerance

                    # Filter df_analysis (OzON data) based on the species, isomer, and the adjusted retention time window
                    filtered_rows = df_analysis[
                        (df_analysis['Species'] == species) &
                        (df_analysis['Isomer'] == isomer) &
                        (df_analysis['Adjusted_RT'] >= retention_time_start) &
                        (df_analysis['Adjusted_RT'] <= retention_time_end)
                    ]

                    # Concatenate the filtered results to the overall filtered_df
                    filtered_df = pd.concat([filtered_df, filtered_rows], ignore_index=True)

                # Construct the output file path for the parquet file using the full Sample value and append `_isomer_filtered_6`
                sample_full_name = df_analysis['Sample'].iloc[0] if 'Sample' in df_analysis.columns else parquet_file.replace('.parquet', '')
                output_file = os.path.join(isomer_filter_output, f'{sample_full_name}_isomer_filtered_6.parquet')

                # Save the filtered DataFrame to a parquet file
                filtered_df.to_parquet(output_file, index=False)
                print(f"Filtered data for {sample_full_name} saved to: {output_file}")
            else:
                print(f"No matching off_possible file found for: {sample_name}")
