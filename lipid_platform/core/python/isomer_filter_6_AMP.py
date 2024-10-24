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
        'Key': ["_".join(os.path.basename(f).split('_')[2:7]).replace('.parquet', '') for f in files2]
    })

    # Remove '5_' prefix from keys if present
    df1['Key'] = df1['Key'].str.replace('^5_', '', regex=True)
    df2['Key'] = df2['Key'].str.replace('^5_', '', regex=True)

    # Print all keys before merging
    print("\nOzOFF Keys:")
    for idx, row in df1.iterrows():
        print(f"{idx + 1}. File: {row['File']}")
        print(f"   Key: {row['Key']}\n")

    print("\nOzON Keys:")
    for idx, row in df2.iterrows():
        print(f"{idx + 1}. File: {row['File']}")
        print(f"   Key: {row['Key']}\n")

    # Merge the DataFrames on the Key column
    merged_df = pd.merge(
        df1, 
        df2, 
        on='Key', 
        how='outer',
        suffixes=('_OzOFF', '_OzON')
    )

    return merged_df

def filter_ozON_by_ozOFF(input_dir, off_possible_dir, retention_time_tolerance=0.3, isomer_filter_output='Projects/AMP/isomer_filter_6/'):
    """
    Iterates over all files in the input directory, filtering each df_analysis (OzON data) based on Retention_Time
    and Species from the corresponding off_possible file. Uses key-based matching logic.
    
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

    print("\nStarting file matching process...")
    print(f"Looking for parquet files in:")
    print(f"OzOFF directory: {off_possible_dir}")
    print(f"OzON directory: {input_dir}")

    # Get matched files using the new function
    file_matches = list_files_in_dirs(off_possible_dir, input_dir)
    
    # Process only rows where we have both OzON and OzOFF files
    matched_files = file_matches.dropna()
    
    print(f"\nFound {len(matched_files)} matched file pairs")
    
    # Print matched pairs
    print("\nMatched Pairs:")
    for idx, row in matched_files.iterrows():
        print(f"\nMatch {idx + 1}:")
        print(f"Key: {row['Key']}")
        print(f"OzOFF: {row['File_OzOFF']}")
        print(f"OzON: {row['File_OzON']}")
    
    for _, row in matched_files.iterrows():
        ozon_file = row['File_OzON']
        ozoff_file = row['File_OzOFF']
        key = row['Key']
        
        print(f"\nProcessing matched pair with key: {key}")
        print(f"OzON file: {ozon_file}")
        print(f"OzOFF file: {ozoff_file}")
        
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
                        
                        print(f"Matching OzON: Lipid={lipid} (n-{n_position}), Adjusted_RT={retention_time_on} "
                              f"with OzOFF: Species={species}, Isomer={isomer_off}, Retention_Time={retention_time_off}")
                        
                        if is_match:
                            print("Result: Match")
                        else:
                            print("Result: No Match")

                # Concatenate the filtered results
                filtered_df = pd.concat([filtered_df, filtered_rows], ignore_index=True)

            # Construct the output file path
            sample_full_name = df_analysis['Sample'].iloc[0] if 'Sample' in df_analysis.columns else ozon_file.replace('.parquet', '')
            output_file = os.path.join(isomer_filter_output, f'{sample_full_name}_isomer_filtered_6.parquet')

            # Save the filtered DataFrame
            filtered_df.to_parquet(output_file, index=False)
            print(f"Filtered data saved to: {output_file}")
            
        except Exception as e:
            print(f"Error processing files with key {key}: {str(e)}")
            continue

    # Report unmatched files
    unmatched = file_matches[file_matches.isnull().any(axis=1)]
    if not unmatched.empty:
        print("\nUnmatched files:")
        for _, row in unmatched.iterrows():
            if pd.isna(row['File_OzON']):
                print(f"OzOFF file without match: {row['File_OzOFF']} (Key: {row['Key']})")
            if pd.isna(row['File_OzOFF']):
                print(f"OzON file without match: {row['File_OzON']} (Key: {row['Key']})")