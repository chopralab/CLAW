import pandas as pd
import os
import glob
import logging
import argparse

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
        pd.DataFrame: DataFrame containing file names from both directories and their keys.
    """
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))

    # Filter files1 if specific_files1 is provided
    if specific_files1:
        specific_files1_full = [os.path.join(dir1, f) for f in specific_files1]
        files1 = [f for f in specific_files1_full if f in files1]
        missing_files1 = set(specific_files1) - set([os.path.basename(f) for f in files1])
        if missing_files1:
            logging.warning(f"Specified files in dir1 not found: {missing_files1}")

    # Filter files2 if specific_files2 is provided
    if specific_files2:
        specific_files2_full = [os.path.join(dir2, f) for f in specific_files2]
        files2 = [f for f in specific_files2_full if f in files2]
        missing_files2 = set(specific_files2) - set([os.path.basename(f) for f in files2])
        if missing_files2:
            logging.warning(f"Specified files in dir2 not found: {missing_files2}")

    def generate_key(filename):
        parts = os.path.basename(filename).split('_')
        if 'FAME' in (part.upper() for part in parts):
            return 'FAME'
        if len(parts) >= 7:
            return "_".join(parts[2:7]).replace('.parquet', '')
        else:
            logging.warning(f"Filename '{os.path.basename(filename)}' does not have enough parts for key generation.")
            return "_".join(parts[2:]).replace('.parquet', '')

    df1 = pd.DataFrame({
        'File_OzOFF': [os.path.basename(f) for f in files1],
        'Key': [generate_key(f) for f in files1]
    })

    df2 = pd.DataFrame({
        'File_OzON': [os.path.basename(f) for f in files2],
        'Key': [generate_key(f) for f in files2]
    })

    df1['Key'] = df1['Key'].astype(str).str.replace('^5_', '', regex=True)
    df2['Key'] = df2['Key'].astype(str).str.replace('^5_', '', regex=True)

    logging.info("\nOzOFF Keys:")
    for idx, row in df1.iterrows():
        logging.info(f"{idx + 1}. File: {row['File_OzOFF']} | Key: {row['Key']}")

    logging.info("\nOzON Keys:")
    for idx, row in df2.iterrows():
        logging.info(f"{idx + 1}. File: {row['File_OzON']} | Key: {row['Key']}")

    merged_df = pd.merge(
        df1,
        df2,
        on='Key',
        how='outer'
    )

    return merged_df

def filter_ozON_by_ozOFF(input_dir, off_possible_dir, output_dir, retention_time_tolerance=0.3, 
                        isomer_filter_output='isomer_filter_output', specific_off_files=None, specific_on_files=None):
    """
    Filters OzON data based on retention time and species from corresponding OzOFF files.

    Parameters:
        input_dir (str): Directory containing OzON parquet files.
        off_possible_dir (str): Directory containing OzOFF parquet files.
        output_dir (str): Directory to save filtered parquet files.
        retention_time_tolerance (float): Tolerance for retention time filtering.
        isomer_filter_output (str): Subdirectory within output_dir for filtered files.
        specific_off_files (list, optional): Specific OzOFF file names to include.
        specific_on_files (list, optional): Specific OzON file names to include.

    Returns:
        None
    """
    os.makedirs(isomer_filter_output, exist_ok=True)

    logging.info("Starting file matching process...")
    logging.info(f"OzOFF directory: {off_possible_dir}")
    logging.info(f"OzON directory: {input_dir}")

    file_matches = list_files_in_dirs(
        off_possible_dir, input_dir, 
        specific_files1=specific_off_files, 
        specific_files2=specific_on_files
    )

    matched_files = file_matches.dropna(subset=['File_OzOFF', 'File_OzON'])
    logging.info(f"\nFound {len(matched_files)} matched file pairs")

    if not matched_files.empty:
        logging.info("\nMatched Pairs:")
        for idx, row in matched_files.iterrows():
            logging.info(f"Match {idx + 1}: Key={row['Key']} | OzOFF={row['File_OzOFF']} | OzON={row['File_OzON']}")
    else:
        logging.info("No matched file pairs found.")
        return

    for _, row in matched_files.iterrows():
        ozon_file = row['File_OzON']
        ozoff_file = row['File_OzOFF']
        key = row['Key']

        logging.info(f"\nProcessing pair Key: {key} | OzON: {ozon_file} | OzOFF: {ozoff_file}")

        try:
            ozon_path = os.path.join(input_dir, ozon_file)
            ozoff_path = os.path.join(off_possible_dir, ozoff_file)

            df_analysis = pd.read_parquet(ozon_path)
            off_possible = pd.read_parquet(ozoff_path)

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

                # Debugging logs
                for _, analysis_row in df_analysis.iterrows():
                    lipid = analysis_row['Lipid']
                    retention_time_on = analysis_row['Adjusted_RT']
                    n_position = analysis_row['n_position']
                    is_match = rt_start <= retention_time_on <= rt_end

                    logging.info(
                        f"Lipid={lipid} (n-{n_position}), Adjusted_RT={retention_time_on} | "
                        f"Species={species}, Isomer={isomer_off}, Retention_Time={retention_time_off} | "
                        f"Match: {'Yes' if is_match else 'No'}"
                    )

            if filtered_df.empty:
                logging.info("No matches found after filtering.")
                continue

            sample_full_name = df_analysis['Sample'].iloc[0] if 'Sample' in df_analysis.columns else os.path.splitext(ozon_file)[0]
            output_file = os.path.join(isomer_filter_output, f"{sample_full_name}_isomer_filtered.parquet")

            filtered_df.to_parquet(output_file, index=False)
            logging.info(f"Filtered data saved to: {output_file}")

        except Exception as e:
            logging.error(f"Error processing Key {key}: {e}")
            continue

    # Report unmatched files
    unmatched = file_matches[file_matches.isnull().any(axis=1)]
    if not unmatched.empty:
        logging.info("\nUnmatched files:")
        for _, row in unmatched.iterrows():
            if pd.isna(row['File_OzON']):
                logging.info(f"OzOFF without match: {row['File_OzOFF']} (Key: {row['Key']})")
            if pd.isna(row['File_OzOFF']):
                logging.info(f"OzON without match: {row['File_OzON']} (Key: {row['Key']})")

def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Filter OzON data based on OzOFF data.")

    parser.add_argument('--input_dir', required=True, help='Directory containing OzON parquet files.')
    parser.add_argument('--off_possible_dir', required=True, help='Directory containing OzOFF parquet files.')
    parser.add_argument('--output_dir', required=True, help='Directory to save filtered parquet files.')
    parser.add_argument('--retention_time_tolerance', type=float, default=0.3, help='Retention time tolerance for filtering.')
    parser.add_argument('--isomer_filter_output', default='isomer_filter_output', help='Subdirectory for filtered output.')
    parser.add_argument('--specific_off_files', nargs='*', help='List of specific OzOFF files to include.')
    parser.add_argument('--specific_on_files', nargs='*', help='List of specific OzON files to include.')

    return parser.parse_args()

def main():
    args = parse_arguments()
    filter_ozON_by_ozOFF(
        input_dir=args.input_dir,
        off_possible_dir=args.off_possible_dir,
        output_dir=args.output_dir,
        retention_time_tolerance=args.retention_time_tolerance,
        isomer_filter_output=args.isomer_filter_output,
        specific_off_files=args.specific_off_files,
        specific_on_files=args.specific_on_files
    )

if __name__ == "__main__":
    main()
