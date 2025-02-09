import pandas as pd
import os
import glob
import re
import argparse
import numpy as np

def list_files_in_dirs(dir1, dir2, extension='*.parquet'):
    """
    Lists file names from two directories and creates keys for matching.
    """
    # Get list of files in both directories
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))

    print(f"Files in {dir1}: {files1}")
    print(f"Files in {dir2}: {files2}")

    def extract_key_from_on(filename):
        basename = os.path.basename(filename)
        match = re.search(r'(NIST_n\d+)', basename)
        if match:
            return match.group(1)
        match = re.search(r'(Blank)', basename)
        if match:
            return match.group(1)
        return basename.replace('_isomer_filtered.parquet', '')

    def extract_key_from_off(filename):
        basename = os.path.basename(filename)
        match = re.search(r'(NIST_n\d+)', basename)
        if match:
            return match.group(1)
        match = re.search(r'(Blank)', basename)
        if match:
            return match.group(1)
        return basename.replace('df_analysis_5_', '').replace('_OFF.parquet', '')

    df1 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files1],
        'Key': [extract_key_from_on(f) for f in files1]
    })
    df2 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files2],
        'Key': [extract_key_from_off(f) for f in files2]
    })

    print(f"Processed DataFrame 1 (ON files):\n{df1}")
    print(f"Processed DataFrame 2 (OFF files):\n{df2}")

    merged_df = pd.merge(
        df1, 
        df2, 
        on='Key', 
        how='inner',  # Only keep pairs with matching keys
        suffixes=('_OzON', '_OzOFF')
    )

    print(f"Merged DataFrame:\n{merged_df}")

    return merged_df

def match_lipids_with_adjusted_rt_to_rt_from_dirs(dir1, dir2, output_dir, rt_window=0.05):
    """
    Matches lipids between files in two directories using species from FAME and retention time.
    First the FAME species are matched (extracted from the Lipid column) and then, for ON rows,
    the OFF retention time for that FAME species is compared to the ON Adjusted_RT.
    Additional debug output is provided.
    """
    file_pairs = list_files_in_dirs(dir1, dir2)

    master_matched = pd.DataFrame()
    master_unmatched = pd.DataFrame()

    # Regular expression to extract species from lipid strings, e.g., "FA(14:1)"
    species_regex = r'(FA\(\d+:\d+\))'

    for _, row in file_pairs.iterrows():
        ozon_file = os.path.join(dir1, row['File_OzON'])
        ozoff_file = os.path.join(dir2, row['File_OzOFF'])

        # Read the files
        ozon_df = pd.read_parquet(ozon_file)
        off_df = pd.read_parquet(ozoff_file)

        print(f"\n--- Processing File Pair ---")
        print(f"OzON file: {row['File_OzON']}")
        print(f"OzOFF file: {row['File_OzOFF']}")
        print(f"OzON DataFrame head:\n{ozon_df.head()}")
        print(f"OzOFF DataFrame head:\n{off_df.head()}")

        # Prepare ON data: create a lower-case lipid and extract species from the Lipid column
        ozon_df['lipid_lower'] = ozon_df['Lipid'].str.strip().str.lower()
        ozon_df['species'] = ozon_df['Lipid'].str.extract(species_regex, flags=re.IGNORECASE)
        ozon_df['species'] = ozon_df['species'].str.lower()

        # Prepare OFF data:
        # For OFF, the Lipid column might contain multiple candidates separated by '|'
        off_df = off_df.copy()
        off_df['off_lipid_list'] = off_df['Lipid'].str.split('|')
        off_exploded = off_df.explode('off_lipid_list').reset_index(drop=True)
        off_exploded['off_lipid'] = off_exploded['off_lipid_list'].str.strip().str.lower()
        # Extract species from the OFF candidate lipid
        off_exploded['species'] = off_exploded['off_lipid'].str.extract(species_regex, flags=re.IGNORECASE)
        off_exploded['species'] = off_exploded['species'].str.lower()
        off_exploded = off_exploded.rename(columns={
            'Lipid': 'Lipid_OFF',
            'Retention_Time': 'Retention_Time_OFF',
            'OzESI_Intensity': 'Intensity_OFF'
        })

        # Reset index in ozon_df so we can track rows for merging
        ozon_df = ozon_df.reset_index().rename(columns={'index': 'orig_index'})

        # Merge ON and OFF on the species field (i.e. FAME species must match)
        merged = pd.merge(ozon_df, off_exploded, on='species', how='left', suffixes=('', '_off'))

        # Debug: Print a sample of the merged DataFrame
        print(f"\nDEBUG: Merged DataFrame head (by species):\n{merged.head(10)}")

        # Compute absolute difference between ON Adjusted_RT and OFF Retention_Time_OFF
        merged['rt_diff'] = (merged['Adjusted_RT'] - merged['Retention_Time_OFF']).abs()

        print(f"DEBUG: Total merged rows before RT filtering: {len(merged)}")

        # Filter valid matches: retention time difference within the window
        valid_matches = merged[merged['rt_diff'] <= rt_window].copy()
        print(f"DEBUG: Valid matches within rt_window (<= {rt_window}) found: {len(valid_matches)}")
        if not valid_matches.empty:
            print(f"DEBUG: Sample valid matches:\n{valid_matches[['Lipid', 'Adjusted_RT', 'Lipid_OFF', 'Retention_Time_OFF', 'rt_diff', 'species']].head(10)}")
        else:
            print("DEBUG: No valid matches found after retention time filtering.")

        # For each original ON row, select the candidate (by species) with the smallest RT difference
        if not valid_matches.empty:
            best_matches = valid_matches.sort_values('rt_diff').groupby('orig_index', as_index=False).first()
        else:
            best_matches = pd.DataFrame(columns=merged.columns)

        if not best_matches.empty:
            print("\nDEBUG: Best matches per ON row (by species and RT difference):")
            for _, m in best_matches.iterrows():
                print(f"   ON row (index {m['orig_index']}): Lipid='{m['Lipid']}', species='{m['species']}', Adjusted_RT={m['Adjusted_RT']} -> OFF Lipid='{m['Lipid_OFF']}', Retention_Time_OFF={m['Retention_Time_OFF']}, rt_diff={m['rt_diff']}")
        else:
            print("DEBUG: No best matches selected for any ON row.")

        # Report ON rows that did not get any match and show best candidate regardless of RT window.
        all_on_indices = set(ozon_df['orig_index'])
        matched_on_indices = set(best_matches['orig_index']) if not best_matches.empty else set()
        unmatched_on_indices = all_on_indices - matched_on_indices
        if unmatched_on_indices:
            print("\nDEBUG: ON rows with no valid matches:")
            for idx in unmatched_on_indices:
                on_row = ozon_df.loc[ozon_df['orig_index'] == idx].iloc[0]
                lipid_val = on_row['Lipid']
                on_rt = on_row['Adjusted_RT']
                species_val = on_row['species']
                print(f"   ON row (index {idx}): Lipid='{lipid_val}', species='{species_val}', Adjusted_RT={on_rt}")

                # Look at all candidate OFF rows for this ON row (ignoring the rt_window filter)
                candidates = merged[(merged['orig_index'] == idx) & (merged['off_lipid'].notnull())]
                if not candidates.empty:
                    best_candidate = candidates.sort_values('rt_diff').iloc[0]
                    off_rt = best_candidate['Retention_Time_OFF']
                    diff = best_candidate['rt_diff']
                    print(f"      -> Best candidate OFF Retention_Time={off_rt} with difference {diff}")
                else:
                    print("      -> No candidate OFF row available for matching.")
        else:
            print("\nDEBUG: All ON rows found a valid match.")

        # Assign the best match information back to ozon_df for the rows that have valid matches.
        ozon_df['Matched_Lipid_OFF'] = np.nan
        ozon_df['Intensity_OFF'] = np.nan
        ozon_df['Retention_Time_OFF'] = np.nan

        if not best_matches.empty:
            for idx in best_matches['orig_index']:
                match_row = best_matches[best_matches['orig_index'] == idx].iloc[0]
                ozon_df.loc[ozon_df['orig_index'] == idx, 'Matched_Lipid_OFF'] = match_row['Lipid_OFF']
                ozon_df.loc[ozon_df['orig_index'] == idx, 'Intensity_OFF'] = match_row['Intensity_OFF']
                ozon_df.loc[ozon_df['orig_index'] == idx, 'Retention_Time_OFF'] = match_row['Retention_Time_OFF']

        # Drop temporary helper columns
        ozon_df = ozon_df.drop(columns=['orig_index', 'lipid_lower', 'species'])

        # Separate matched and unmatched rows
        matched_lipids_df = ozon_df.dropna(subset=['Matched_Lipid_OFF']).copy()
        unmatched_lipids_df = ozon_df[ozon_df['Matched_Lipid_OFF'].isna()].copy()

        os.makedirs(output_dir, exist_ok=True)
        matched_file = os.path.join(output_dir, f"matched_{row['Key']}.csv")
        unmatched_file = os.path.join(output_dir, f"unmatched_{row['Key']}.csv")

        matched_lipids_df.to_csv(matched_file, index=False)
        unmatched_lipids_df.to_csv(unmatched_file, index=False)

        print(f"\nProcessed files: {row['File_OzON']} and {row['File_OzOFF']}")
        print(f"Saved matched lipids to {matched_file}")
        print(f"Saved unmatched lipids to {unmatched_file}")

        # Append to master DataFrames
        master_matched = pd.concat([master_matched, matched_lipids_df], ignore_index=True)
        master_unmatched = pd.concat([master_unmatched, unmatched_lipids_df], ignore_index=True)

    # Save master DataFrames as Parquet files
    master_matched_file = os.path.join(output_dir, "master_matched.parquet")
    master_unmatched_file = os.path.join(output_dir, "master_unmatched.parquet")

    master_matched.to_parquet(master_matched_file, index=False)
    master_unmatched.to_parquet(master_unmatched_file, index=False)

    print(f"\nSaved master matched DataFrame to {master_matched_file}")
    print(f"Saved master unmatched DataFrame to {master_unmatched_file}")

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Match lipids between OzON and OzOFF directories based on FAME species and retention time."
    )
    parser.add_argument('--ozon_dir', type=str, required=True, help='Path to the OzON directory.')
    parser.add_argument('--ozoff_dir', type=str, required=True, help='Path to the OzOFF directory.')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to the output directory.')
    parser.add_argument('--rt_window', type=float, default=0.05, help='Retention time window for matching (default: 0.05).')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    match_lipids_with_adjusted_rt_to_rt_from_dirs(
        dir1=args.ozon_dir,
        dir2=args.ozoff_dir,
        output_dir=args.output_dir,
        rt_window=args.rt_window
    )
