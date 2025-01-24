import pandas as pd
import os
import glob
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def list_files_in_dirs(dir1, dir2, extension='*.parquet', specific_files1=None, specific_files2=None):
    """
    Lists file names from two directories and creates keys for matching.

    Parameters:
    dir1 (str): Path to the first directory.
    dir2 (str): Path to the second directory.
    extension (str): File extension to match (default: '*.parquet').
    specific_files1 (list, optional): Specific file names from dir1 to include.
    specific_files2 (list, optional): Specific file names from dir2 to include.

    Returns:
    DataFrame: DataFrame containing file names from both directories and their keys.
    """
    # Get list of files in both directories
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))

    # If specific_files1 is provided, filter files1
    if specific_files1 is not None:
        # Create full paths for specific files
        specific_files1_full = [os.path.join(dir1, f) for f in specific_files1]
        # Filter files1 to include only those specified
        files1 = [f for f in specific_files1_full if f in files1]
        # Identify missing files
        missing_files1 = set(specific_files1) - set([os.path.basename(f) for f in files1])
        if missing_files1:
            logging.warning(f"The following specified files in dir1 were not found: {missing_files1}")

    # If specific_files2 is provided, filter files2
    if specific_files2 is not None:
        # Create full paths for specific files
        specific_files2_full = [os.path.join(dir2, f) for f in specific_files2]
        # Filter files2 to include only those specified
        files2 = [f for f in specific_files2_full if f in files2]
        # Identify missing files
        missing_files2 = set(specific_files2) - set([os.path.basename(f) for f in files2])
        if missing_files2:
            logging.warning(f"The following specified files in dir2 were not found: {missing_files2}")

    # Function to generate keys safely
    def generate_key(filename):
        parts = os.path.basename(filename).split('_')
        # Extract 'FAME' or other lipid names as needed
        for part in parts:
            if part.upper() == 'FAME':
                return 'FAME'
        # Fallback to existing key generation logic
        if len(parts) >= 7:
            key = "_".join(parts[2:7]).replace('.parquet', '')
        else:
            key = "_".join(parts[2:]).replace('.parquet', '')
            logging.warning(f"Filename '{os.path.basename(filename)}' does not have enough parts for key generation.")
        return key

    # Create separate DataFrames for each directory
    df1 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files1],
        'Key': [generate_key(f) for f in files1]
    })
    df2 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files2],
        'Key': [generate_key(f) for f in files2]
    })

    # Convert 'Key' to string type to avoid AttributeError
    df1['Key'] = df1['Key'].astype(str).str.replace('^5_', '', regex=True)
    df2['Key'] = df2['Key'].astype(str).str.replace('^5_', '', regex=True)

    # Print all keys before merging
    logging.info("\nOzOFF Keys:")
    for idx, row in df1.iterrows():
        logging.info(f"{idx + 1}. File: {row['File']} | Key: {row['Key']}")

    logging.info("\nOzON Keys:")
    for idx, row in df2.iterrows():
        logging.info(f"{idx + 1}. File: {row['File']} | Key: {row['Key']}")

    # Merge the DataFrames on the Key column
    merged_df = pd.merge(
        df1,
        df2,
        on='Key',
        how='outer',
        suffixes=('_OzOFF', '_OzON')
    )

    return merged_df

def filter_ozON_by_ozOFF(input_dir, off_possible_dir, retention_time_tolerance=0.3, isomer_filter_output='Projects/AMP/isomer_filter_6/',
                         specific_off_files=None, specific_on_files=None):
    """
    Iterates over all files in the input directory, filtering each df_analysis (OzON data) based on Retention_Time
    and Species from the corresponding off_possible file. Uses key-based matching logic.

    Additionally, allows specifying specific file names from each input directory to match.

    Parameters:
    input_dir (str): Directory containing parquet files to process.
    off_possible_dir (str): Directory containing parquet files for off_possible data.
    retention_time_tolerance (float): The tolerance for filtering retention times (default: 0.3).
    isomer_filter_output (str): Path where the filtered parquet files will be saved.
    specific_off_files (list, optional): List of specific OzOFF file names to include.
    specific_on_files (list, optional): List of specific OzON file names to include.

    Returns:
    None: The function saves the filtered DataFrame as parquet files in the specified output path.
    """
    # Create the output directory if it does not exist
    os.makedirs(isomer_filter_output, exist_ok=True)

    logging.info("\nStarting file matching process...")
    logging.info(f"Looking for parquet files in:")
    logging.info(f"OzOFF directory: {off_possible_dir}")
    logging.info(f"OzON directory: {input_dir}")

    # Get matched files using the updated function with specific file parameters
    file_matches = list_files_in_dirs(off_possible_dir, input_dir, specific_files1=specific_off_files, specific_files2=specific_on_files)
    
    # Process only rows where we have both OzON and OzOFF files
    matched_files = file_matches.dropna()

    logging.info(f"\nFound {len(matched_files)} matched file pairs")
    
    # Print matched pairs
    if not matched_files.empty:
        logging.info("\nMatched Pairs:")
        for idx, row in matched_files.iterrows():
            logging.info(f"\nMatch {idx + 1}:")
            logging.info(f"Key: {row['Key']}")
            logging.info(f"OzOFF: {row['File_OzOFF']}")
            logging.info(f"OzON: {row['File_OzON']}")
    else:
        logging.info("\nNo matched file pairs found.")

    for _, row in matched_files.iterrows():
        ozon_file = row['File_OzON']
        ozoff_file = row['File_OzOFF']
        key = row['Key']
        
        logging.info(f"\nProcessing matched pair with key: {key}")
        logging.info(f"OzON file: {ozon_file}")
        logging.info(f"OzOFF file: {ozoff_file}")
        
        try:
            # Read the df_analysis file (OzON data)
            ozon_path = os.path.join(input_dir, ozon_file)
            df_analysis = pd.read_parquet(ozon_path)
            
            # Read the off_possible file
            off_possible_path = os.path.join(off_possible_dir, ozoff_file)
            off_possible = pd.read_parquet(off_possible_path)
            
            # Create Adjusted_RT column
            if 'STD_RT_Dif' in df_analysis.columns:
                df_analysis['Adjusted_RT'] = df_analysis['Retention_Time'] + df_analysis['STD_RT_Dif']
            else:
                df_analysis['Adjusted_RT'] = df_analysis['Retention_Time']

            # Extract Species from Lipid names if not already present
            if 'Species' not in df_analysis.columns:
                df_analysis['Species'] = df_analysis['Lipid'].str.extract(r'\((\d+:\d+)\)')[0]

            # Extract n-position from Lipid names
            df_analysis['n_position'] = df_analysis['Lipid'].str.extract(r'n-(\d+)')[0].astype(int)

            # Initialize an empty DataFrame to store the filtered results
            filtered_df = pd.DataFrame()

            # Iterate over each row in the off_possible (OzOFF data)
            for _, off_row in off_possible.iterrows():
                species = off_row['Species']
                isomer_off = off_row['Isomer']
                retention_time_off = off_row['Retention_Time']

                # Define the retention time window
                retention_time_start = retention_time_off - retention_time_tolerance
                retention_time_end = retention_time_off + retention_time_tolerance

                # Filter df_analysis based on species and retention time window
                filtered_rows = df_analysis[
                    (df_analysis['Species'] == species) &
                    (df_analysis['Adjusted_RT'] >= retention_time_start) &
                    (df_analysis['Adjusted_RT'] <= retention_time_end)
                ]

                # Add the OzOFF isomer number to the matched rows
                if not filtered_rows.empty:
                    filtered_rows = filtered_rows.assign(OzOFF_Isomer=isomer_off)
                
                # Print matching details for debugging
                for _, analysis_row in df_analysis.iterrows():
                    if analysis_row['Species'] == species:
                        lipid = analysis_row['Lipid']
                        retention_time_on = analysis_row['Adjusted_RT']
                        n_position = analysis_row['n_position']
                        
                        is_match = (retention_time_start <= retention_time_on <= retention_time_end)
                        
                        logging.info(f"Matching OzON: Lipid={lipid} (n-{n_position}), Adjusted_RT={retention_time_on} "
                                     f"with OzOFF: Species={species}, Isomer={isomer_off}, Retention_Time={retention_time_off}")
                        
                        if is_match:
                            logging.info("Result: Match")
                        else:
                            logging.info("Result: No Match")

                # Concatenate the filtered results
                filtered_df = pd.concat([filtered_df, filtered_rows], ignore_index=True)

            # Construct the output file path
            sample_full_name = df_analysis['Sample'].iloc[0] if 'Sample' in df_analysis.columns else ozon_file.replace('.parquet', '')
            output_file = os.path.join(isomer_filter_output, f'{sample_full_name}_isomer_filtered_6.parquet')

            # Save the filtered DataFrame
            filtered_df.to_parquet(output_file, index=False)
            logging.info(f"Filtered data saved to: {output_file}")
            
        except Exception as e:
            logging.error(f"Error processing files with key {key}: {str(e)}")
            continue

    # Report unmatched files
    unmatched = file_matches[file_matches.isnull().any(axis=1)]
    if not unmatched.empty:
        logging.info("\nUnmatched files:")
        for _, row in unmatched.iterrows():
            if pd.isna(row['File_OzON']):
                logging.info(f"OzOFF file without match: {row['File_OzOFF']} (Key: {row['Key']})")
            if pd.isna(row['File_OzOFF']):
                logging.info(f"OzON file without match: {row['File_OzON']} (Key: {row['Key']})")

# Example Usage:


