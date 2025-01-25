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

    # Debugging: Print file lists
    print(f"Files in {dir1}: {files1}")
    print(f"Files in {dir2}: {files2}")



    # Extract keys for OFF files
    def extract_key_from_off(filename):
        parts = os.path.basename(filename).replace('df_analysis_5_', '').replace('_OFF.parquet', '').split('_')
        return "_".join(parts[:5])  # Ensure matching format as ON files
    
    # Extract keys for ON files
    def extract_key_from_on(filename):
        parts = os.path.basename(filename).split('_')
        if len(parts) >= 4:
            return "_".join(parts[:4])  # Extract first 4 parts: region_condition_mouseID_sampleID
        return filename  # Fallback to full filename if parsing fails

 


    # Create DataFrames with extracted keys
    df1 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files1],
        'Key': [extract_key_from_on(f) for f in files1]
    })
    df2 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files2],
        'Key': [extract_key_from_off(f) for f in files2]
    })

    # Debugging: Print DataFrames before merging
    print(f"Processed DataFrame 1 (ON files):\n{df1}")
    print(f"Processed DataFrame 2 (OFF files):\n{df2}")

    # Merge the DataFrames on the Key column
    merged_df = pd.merge(
        df1, 
        df2, 
        on='Key', 
        how='inner',  # Only keep pairs with matching keys
        suffixes=('_OzON', '_OzOFF')
    )

    # Debugging: Print merged DataFrame
    print(f"Merged DataFrame:\n{merged_df}")

    return merged_df



def match_lipids_with_adjusted_rt_to_rt_from_dirs(dir1, dir2, output_dir, rt_window=0.05):
    """
    Matches lipids between files in two directories using retention time and keys.
    Appends results to master matched and unmatched DataFrames and saves them as Parquet files.

    Parameters:
    dir1 (str): Directory containing OzON files.
    dir2 (str): Directory containing OzOFF files.
    output_dir (str): Directory to save matched and unmatched results.
    rt_window (float): Retention time window for matching (default: 0.5).
    """
    file_pairs = list_files_in_dirs(dir1, dir2)

    # Initialize master DataFrames
    master_matched = pd.DataFrame()
    master_unmatched = pd.DataFrame()

    for _, row in file_pairs.iterrows():
        ozon_file = os.path.join(dir1, row['File_OzON'])
        ozoff_file = os.path.join(dir2, row['File_OzOFF'])

        ozon_test = pd.read_parquet(ozon_file)
        sorted_notpossible_lipids = pd.read_parquet(ozoff_file)

        print(f"OzON DataFrame for {row['File_OzON']}:\n{ozon_test.head()}")
        print(f"OzOFF DataFrame for {row['File_OzOFF']}:\n{sorted_notpossible_lipids.head()}")

        ozon_test['Matched_Lipid_OFF'] = None
        ozon_test['Intensity_OFF'] = None
        ozon_test['Retention_Time_OFF'] = None

        for index_on, row_on in ozon_test.iterrows():
            on_lipid = row_on['Lipid'].strip()
            on_adjusted_rt = row_on['Adjusted_RT']
            match_found = False

            for _, row_off in sorted_notpossible_lipids.iterrows():
                off_lipids = [lipid.strip() for lipid in row_off['Lipid'].split('|')]
                off_retention_time = row_off['Retention_Time']

                print(f"Comparing ON lipid {on_lipid} (RT: {on_adjusted_rt}) with OFF lipids {off_lipids} (RT: {off_retention_time})")

                if any(on_lipid.lower() == off_lipid.lower() for off_lipid in off_lipids) and abs(on_adjusted_rt - off_retention_time) <= rt_window:
                    ozon_test.at[index_on, 'Matched_Lipid_OFF'] = row_off['Lipid']
                    ozon_test.at[index_on, 'Intensity_OFF'] = row_off['OzESI_Intensity']
                    ozon_test.at[index_on, 'Retention_Time_OFF'] = off_retention_time
                    match_found = True
                    break

            if not match_found:
                print(f"No match found for ON lipid {on_lipid} with Adjusted RT {on_adjusted_rt}")

        matched_lipids_df = ozon_test.dropna(subset=['Matched_Lipid_OFF']).copy()
        unmatched_lipids_df = ozon_test[ozon_test['Matched_Lipid_OFF'].isna()].copy()

        os.makedirs(output_dir, exist_ok=True)
        matched_file = os.path.join(output_dir, f"matched_{row['Key']}.csv")
        unmatched_file = os.path.join(output_dir, f"unmatched_{row['Key']}.csv")

        matched_lipids_df.to_csv(matched_file, index=False)
        unmatched_lipids_df.to_csv(unmatched_file, index=False)

        # Append to master DataFrames
        master_matched = pd.concat([master_matched, matched_lipids_df], ignore_index=True)
        master_unmatched = pd.concat([master_unmatched, unmatched_lipids_df], ignore_index=True)

        print(f"Processed files: {row['File_OzON']} and {row['File_OzOFF']}")
        print(f"Saved matched lipids to {matched_file}")
        print(f"Saved unmatched lipids to {unmatched_file}")

    # Save master DataFrames as Parquet files
    master_matched_file = os.path.join(output_dir, "master_matched.parquet")
    master_unmatched_file = os.path.join(output_dir, "master_unmatched.parquet")

    master_matched.to_parquet(master_matched_file, index=False)
    master_unmatched.to_parquet(master_unmatched_file, index=False)

    print(f"Saved master matched DataFrame to {master_matched_file}")
    print(f"Saved master unmatched DataFrame to {master_unmatched_file}")



if __name__ == "__main__":
    # Example usage
    ozon_dir = 'Projects/AMP/isomer_filter_6/'
    ozoff_dir = 'Projects/AMP/analysis/OFF/notpossible/'
    output_dir = 'Projects/AMP/notpossible_9/'

    match_lipids_with_adjusted_rt_to_rt_from_dirs(ozon_dir, ozoff_dir, output_dir, rt_window=0.05)
