import pandas as pd
import os
import glob
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def list_files_in_dirs(dir1, dir2, extension='*.parquet'):
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))
    
    logging.debug(f"Files in {dir1}: {files1}")
    logging.debug(f"Files in {dir2}: {files2}")
    
    def get_key(filename):
        base_name = os.path.basename(filename).replace('.parquet', '')
        return "FAME" if "FAME" in base_name else ("CisTrans" if "CisTrans" in base_name else base_name)
    
    df1 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files1],
        'Key': [get_key(f) for f in files1]
    })
    df2 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files2],
        'Key': [get_key(f) for f in files2]
    })
    
    logging.debug(f"Initial DataFrame 1:\n{df1}")
    logging.debug(f"Initial DataFrame 2:\n{df2}")
    
    merged_df = pd.merge(df1, df2, on='Key', how='inner', suffixes=('_OzON', '_OzOFF'))
    
    logging.debug(f"Merged DataFrame:\n{merged_df}")
    
    return merged_df

def match_lipids_with_adjusted_rt_to_rt_from_dirs(dir1, dir2, output_dir, rt_window=0.5):
    file_pairs = list_files_in_dirs(dir1, dir2)
    
    for _, row in file_pairs.iterrows():
        ozon_file = os.path.join(dir1, row['File_OzON'])
        ozoff_file = os.path.join(dir2, row['File_OzOFF'])
        
        logging.debug(f"Processing OzON file: {ozon_file}")
        logging.debug(f"Processing OzOFF file: {ozoff_file}")
        
        ozon_test = pd.read_parquet(ozon_file)
        sorted_notpossible_lipids = pd.read_parquet(ozoff_file)
        
        logging.debug(f"OzON DataFrame (first 5 rows):\n{ozon_test.head()}")
        logging.debug(f"OzOFF DataFrame (first 5 rows):\n{sorted_notpossible_lipids.head()}")
        
        ozon_test['Matched_Lipid_OFF'] = None
        ozon_test['Intensity_OFF'] = None
        ozon_test['Retention_Time_OFF'] = None

        for index_on, row_on in ozon_test.iterrows():
            on_lipid = row_on['Lipid'].strip()
            on_adjusted_rt = row_on['Adjusted_RT']

            for _, row_off in sorted_notpossible_lipids.iterrows():
                off_lipids = [lipid.strip() for lipid in row_off['Lipid'].split('|')]
                off_retention_time = row_off['Retention_Time']

                if any(on_lipid.lower() == off_lipid.lower() for off_lipid in off_lipids) and abs(on_adjusted_rt - off_retention_time) <= rt_window:
                    ozon_test.at[index_on, 'Matched_Lipid_OFF'] = row_off['Lipid']
                    ozon_test.at[index_on, 'Intensity_OFF'] = row_off['OzESI_Intensity']
                    ozon_test.at[index_on, 'Retention_Time_OFF'] = off_retention_time
                    logging.debug(f"Matched {on_lipid} (Adjusted_RT: {on_adjusted_rt}) with {row_off['Lipid']} (Retention_Time: {off_retention_time})")
                    break

        matched_lipids_df = ozon_test.dropna(subset=['Matched_Lipid_OFF']).copy()
        unmatched_lipids_df = ozon_test[ozon_test['Matched_Lipid_OFF'].isna()].copy()
        
        os.makedirs(output_dir, exist_ok=True)
        matched_file = os.path.join(output_dir, f"matched_{row['Key']}.csv")
        unmatched_file = os.path.join(output_dir, f"unmatched_{row['Key']}.csv")
        
        matched_lipids_df.to_csv(matched_file, index=False)
        unmatched_lipids_df.to_csv(unmatched_file, index=False)
        
        logging.debug(f"Saved matched results to {matched_file}")
        logging.debug(f"Saved unmatched results to {unmatched_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Match lipids between OzON and OzOFF datasets.")
    parser.add_argument("--ozon_dir", type=str, required=True, help="Directory containing OzON files.")
    parser.add_argument("--ozoff_dir", type=str, required=True, help="Directory containing OzOFF files.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument("--rt_window", type=float, default=0.5, help="Retention time window for matching (default: 0.5).")
    
    args = parser.parse_args()
    match_lipids_with_adjusted_rt_to_rt_from_dirs(args.ozon_dir, args.ozoff_dir, args.output_dir, args.rt_window)
