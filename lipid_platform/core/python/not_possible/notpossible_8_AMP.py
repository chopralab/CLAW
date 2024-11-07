import pandas as pd
import os
import glob

def list_files_in_dirs(dir1, dir2, extension='*.parquet'):
    """
    Lists file names from two directories and creates keys for matching.

    Parameters:
    dir1 (str): Path to the first directory.
    dir2 (str): Path to the second directory.
    extension (str): File extension to match (default: '*.parquet').

    Returns:
    DataFrame: DataFrame containing file names from both directories and their keys.
    """
    # Get list of files in both directories
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))

    # Create separate DataFrames for each directory
    df1 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files1],
        'Key': ["_".join(os.path.basename(f).split('_')[2:7]).replace('.parquet', '') for f in files1]
    })
    df2 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files2],
        'Key': ["_".join(os.path.basename(f).split('_')[0:4]).replace('.parquet', '') for f in files2]
    })

    # Remove '5_' prefix from keys if present
    df1['Key'] = df1['Key'].str.replace('^5_', '', regex=True)
    df2['Key'] = df2['Key'].str.replace('^5_', '', regex=True)

    # Print out the files and keys from both directories
    print("\nFiles in Directory 1 (OzOFF):")
    for _, row in df1.iterrows():
        print(f"File: {row['File']}, Key: {row['Key']}")

    print("\nFiles in Directory 2 (OzON):")
    for _, row in df2.iterrows():
        print(f"File: {row['File']}, Key: {row['Key']}")

    # Merge the DataFrames on the Key column
    merged_df = pd.merge(
        df1, 
        df2, 
        on='Key', 
        how='outer',
        suffixes=('_OzOFF', '_OzON')
    )

    return merged_df

def match_lipids_with_adjusted_rt_to_rt(input_dir, off_possible_dir, rt_window=0.5, output_matched_dir='', output_unmatched_dir='', output_log_dir=''):
    """
    Matches lipids between OzON and OzOFF DataFrames based on adjusted retention time by matching file names.

    Parameters:
    - input_dir: str, directory containing OzON parquet files
    - off_possible_dir: str, directory containing OzOFF parquet files
    - rt_window: float, the retention time window for matching
    - output_matched_dir: str, directory to save matched DataFrames
    - output_unmatched_dir: str, directory to save unmatched DataFrames
    - output_log_dir: str, directory to save log files

    Returns:
    - None
    """
    # Ensure output directories exist
    if output_matched_dir:
        os.makedirs(output_matched_dir, exist_ok=True)
        print(f"Output directory for matched files: {output_matched_dir}")
    if output_unmatched_dir:
        os.makedirs(output_unmatched_dir, exist_ok=True)
        print(f"Output directory for unmatched files: {output_unmatched_dir}")
    if output_log_dir:
        os.makedirs(output_log_dir, exist_ok=True)
        print(f"Output directory for log files: {output_log_dir}")

    # Get matched files using key-based matching
    file_matches = list_files_in_dirs(off_possible_dir, input_dir)

    # Process only rows where we have both OzON and OzOFF files
    matched_files = file_matches.dropna()
    print(f"Total matched files: {len(matched_files)}")

    # Print out matched and unmatched files
    unmatched_files = file_matches[file_matches.isna().any(axis=1)]
    if not unmatched_files.empty:
        print("\nUnmatched Files:")
        for _, row in unmatched_files.iterrows():
            if pd.isna(row['File_OzON']):
                print(f"OzOFF file without match: {row['File_OzOFF']} (Key: {row['Key']})")
            if pd.isna(row['File_OzOFF']):
                print(f"OzON file without match: {row['File_OzON']} (Key: {row['Key']})")

    for _, row in matched_files.iterrows():
        ozon_file = row['File_OzON']
        ozoff_file = row['File_OzOFF']
        key = row['Key']

        log_messages = []

        try:
            # Read the OzON (ON) data
            ozon_path = os.path.join(input_dir, ozon_file)
            ozon_test = pd.read_parquet(ozon_path)
            log_messages.append(f"Successfully read OzON file: {ozon_path}")

            # Read the OzOFF (OFF) data
            off_possible_path = os.path.join(off_possible_dir, ozoff_file)
            sorted_notpossible_lipids = pd.read_parquet(off_possible_path)
            log_messages.append(f"Successfully read OzOFF file: {off_possible_path}")

            # Create new columns in the ON data to store matched information from OFF
            ozon_test['Matched_Lipid_OFF'] = None
            ozon_test['Intensity_OFF'] = None
            ozon_test['Retention_Time_OFF'] = None

            # Iterate through each row in the ON data
            for index_on, row_on in ozon_test.iterrows():
                on_lipid = row_on['Lipid'].strip()  # Trim spaces for accurate matching
                on_adjusted_rt = row_on['Adjusted_RT']  # Adjusted RT of the ON lipid
                match_found = False

                # Iterate through each row in the OFF data
                for index_off, row_off in sorted_notpossible_lipids.iterrows():
                    # Split the OFF Lipid values by '|' and trim each part
                    off_lipids = [lipid.strip() for lipid in row_off['Lipid'].split('|')]
                    off_retention_time = row_off['Retention_Time']  # Retention time of the OFF lipid

                    # Check if the ON lipid matches any of the possible OFF lipids
                    # and if the Adjusted_RT is within the specified window of the Retention_Time
                    if any(on_lipid.lower() == off_lipid.lower() for off_lipid in off_lipids) and abs(on_adjusted_rt - off_retention_time) <= rt_window:
                        # If a match is found, update the 'Matched_Lipid_OFF', 'Intensity_OFF', and 'Retention_Time_OFF' columns in the ON DataFrame
                        ozon_test.at[index_on, 'Matched_Lipid_OFF'] = row_off['Lipid']
                        ozon_test.at[index_on, 'Intensity_OFF'] = row_off['OzESI_Intensity']
                        ozon_test.at[index_on, 'Retention_Time_OFF'] = off_retention_time
                        match_found = True
                        break  # Stop searching for this ON lipid once a match is found

                # If no match is found, log or handle this scenario
                if not match_found:
                    log_messages.append(f"No match found for ON lipid: {on_lipid} within Adjusted_RT window {on_adjusted_rt} ± {rt_window}")

            # Create DataFrames for matched and unmatched lipids
            matched_lipids = ozon_test.dropna(subset=['Matched_Lipid_OFF']).copy()
            unmatched_lipids = ozon_test[ozon_test['Matched_Lipid_OFF'].isna()].copy()

            # Save matched and unmatched DataFrames to the specified output directories
            if output_matched_dir and not matched_lipids.empty:
                matched_filename = ozon_file.replace('_fame_filter_7.parquet', '_notpossible_maybe')
                matched_output_path_parquet = os.path.join(output_matched_dir, matched_filename + '.parquet')
                matched_output_path_csv = os.path.join(output_matched_dir, matched_filename + '.csv')
                try:
                    matched_lipids.to_parquet(matched_output_path_parquet, index=False)
                    matched_lipids.to_csv(matched_output_path_csv, index=False)
                    log_messages.append(f"Matched data saved to: {matched_output_path_parquet} and {matched_output_path_csv}")
                except Exception as e:
                    log_messages.append(f"Failed to save matched data, Error: {str(e)}")

            if output_unmatched_dir and not unmatched_lipids.empty:
                unmatched_filename = ozon_file.replace('_fame_filter_7.parquet', '_notpossible_yes')
                unmatched_output_path_parquet = os.path.join(output_unmatched_dir, unmatched_filename + '.parquet')
                unmatched_output_path_csv = os.path.join(output_unmatched_dir, unmatched_filename + '.csv')
                try:
                    unmatched_lipids.to_parquet(unmatched_output_path_parquet, index=False)
                    unmatched_lipids.to_csv(unmatched_output_path_csv, index=False)
                    log_messages.append(f"Unmatched data saved to: {unmatched_output_path_parquet} and {unmatched_output_path_csv}")
                except Exception as e:
                    log_messages.append(f"Failed to save unmatched data, Error: {str(e)}")

        except Exception as e:
            log_messages.append(f"Error processing files with key {key}: {str(e)}")

        # Save log messages to a log file
        if output_log_dir:
            log_filename = f"{key}_notpossibleLOG.txt"
            log_output_path = os.path.join(output_log_dir, log_filename)
            try:
                with open(log_output_path, 'w') as log_file:
                    log_file.write("\n".join(log_messages))
                print(f"Log saved to: {log_output_path}")
            except Exception as e:
                print(f"Failed to save log to: {log_output_path}, Error: {str(e)}")
