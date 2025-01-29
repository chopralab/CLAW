import pandas as pd
import os
import glob
import logging
import argparse

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def generate_key(filename):
    """
    Generates a matching key from filename with proper normalization for 'CT_OFF_', '5_', and 'max_' cases.
    """
    basename = os.path.basename(filename).replace('df_analysis_', '').replace('.parquet', '')

    # Special handling for CisTrans files
    if 'CisTrans' in basename:
        return 'MAX_CISTRANS' if 'max_' in basename else 'CISTRANS'

    # Special handling for FAME files
    if 'FAME' in basename:
        return 'MAX_FAME' if 'max_' in basename else 'FAME'

    # Remove 'CT_OFF_' and '5_' prefixes
    basename = basename.replace('CT_OFF_', '').replace('5_', '')

    # Convert 'max_' to 'MAX_'
    if 'max_' in basename:
        basename = basename.replace('max_', 'MAX_')

    return basename.upper()

def list_files_in_dirs(dir1, dir2, extension='*.parquet', specific_files1=None, specific_files2=None):
    """
    Lists file names from two directories and matches them based on generated keys.
    Includes debug logging for matching process.
    """
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))

    logging.info(f"Found {len(files1)} files in {dir1}:")
    for f in files1:
        logging.info(f"  - {f}")

    logging.info(f"Found {len(files2)} files in {dir2}:")
    for f in files2:
        logging.info(f"  - {f}")

    df1 = pd.DataFrame({'File_OzOFF': [os.path.basename(f) for f in files1], 'Key': [generate_key(f) for f in files1]})
    df2 = pd.DataFrame({'File_OzON': [os.path.basename(f) for f in files2], 'Key': [generate_key(f) for f in files2]})

    merged_df = pd.merge(df1, df2, on='Key', how='outer')

    matched_files = merged_df.dropna(subset=['File_OzOFF', 'File_OzON'])
    unmatched_off = merged_df[pd.isna(merged_df['File_OzON'])]
    unmatched_on = merged_df[pd.isna(merged_df['File_OzOFF'])]

    logging.info(f"\nTotal matched file pairs: {len(matched_files)}")
    for _, row in matched_files.iterrows():
        logging.info(f"Matched Pair -> OzOFF: {row['File_OzOFF']} <-> OzON: {row['File_OzON']} (Key: {row['Key']})")

    if not unmatched_off.empty:
        logging.warning("\nUnmatched OzOFF files:")
        for _, row in unmatched_off.iterrows():
            logging.warning(f"OzOFF without match: {row['File_OzOFF']} (Key: {row['Key']})")

    if not unmatched_on.empty:
        logging.warning("\nUnmatched OzON files:")
        for _, row in unmatched_on.iterrows():
            logging.warning(f"OzON without match: {row['File_OzON']} (Key: {row['Key']})")

    return merged_df

def filter_ozON_by_ozOFF(input_dir, off_possible_dir, output_dir, retention_time_tolerance=0.3, 
                          specific_off_files=None, specific_on_files=None):
    """
    Filters OzON data based on retention time and species from corresponding OzOFF files.
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure OUTPUT_DIR exists

    logging.info(f"Scanning OzOFF directory: {off_possible_dir}")
    logging.info(f"Scanning OzON directory: {input_dir}")

    file_matches = list_files_in_dirs(off_possible_dir, input_dir, 
                                      specific_files1=specific_off_files, 
                                      specific_files2=specific_on_files)

    matched_files = file_matches.dropna(subset=['File_OzOFF', 'File_OzON'])
    if matched_files.empty:
        logging.info("No matched file pairs found.")
        return

    for _, row in matched_files.iterrows():
        ozon_file = row['File_OzON']
        ozoff_file = row['File_OzOFF']
        key = row['Key']

        logging.info(f"\nProcessing matched files for Key: {key}")
        logging.info(f"OzON File: {ozon_file}")
        logging.info(f"OzOFF File: {ozoff_file}")

        try:
            df_analysis = pd.read_parquet(os.path.join(input_dir, ozon_file))
            off_possible = pd.read_parquet(os.path.join(off_possible_dir, ozoff_file))

            df_analysis['Adjusted_RT'] = df_analysis.get('Retention_Time', 0) + df_analysis.get('STD_RT_Dif', 0)

            if 'Species' not in df_analysis.columns:
                df_analysis['Species'] = df_analysis['Lipid'].str.extract(r'\((\d+:\d+)\)')[0]

            df_analysis['n_position'] = df_analysis['Lipid'].str.extract(r'n-(\d+)')[0].astype(float)

            filtered_df = pd.DataFrame()

            for _, off_row in off_possible.iterrows():
                species = off_row['Species']
                isomer_off = off_row['Isomer']
                retention_time_off = off_row['Retention_Time']

                rt_start = retention_time_off - retention_time_tolerance
                rt_end = retention_time_off + retention_time_tolerance

                condition = (
                    (df_analysis['Species'] == species) &
                    (df_analysis['Adjusted_RT'] >= rt_start) &
                    (df_analysis['Adjusted_RT'] <= rt_end)
                )
                matches = df_analysis[condition].copy()
                
                if not matches.empty:
                    matches['OzOFF_Isomer'] = isomer_off
                    filtered_df = pd.concat([filtered_df, matches], ignore_index=True)

            if filtered_df.empty:
                logging.info(f"No retention time matches found for {ozon_file}")
                continue

            sample_full_name = df_analysis['Sample'].iloc[0] if 'Sample' in df_analysis.columns else os.path.splitext(ozon_file)[0]
            output_file = os.path.join(output_dir, f"{sample_full_name}_isomer_filtered.parquet")

            filtered_df.to_parquet(output_file, index=False)
            logging.info(f"Filtered data saved to: {output_file}")

        except Exception as e:
            logging.error(f"Error processing {ozon_file}: {str(e)}")
            continue

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter OzON data based on OzOFF data.")
    parser.add_argument('--input_dir', required=True, help='Directory containing OzON parquet files.')
    parser.add_argument('--off_possible_dir', required=True, help='Directory containing OzOFF parquet files.')
    parser.add_argument('--output_dir', required=True, help='Directory to save filtered parquet files.')
    parser.add_argument('--retention_time_tolerance', type=float, default=0.3, help='Retention time tolerance for filtering.')
    return parser.parse_args()

def main():
    args = parse_arguments()
    filter_ozON_by_ozOFF(
        input_dir=args.input_dir,
        off_possible_dir=args.off_possible_dir,
        output_dir=args.output_dir,
        retention_time_tolerance=args.retention_time_tolerance
    )

if __name__ == "__main__":
    main()
