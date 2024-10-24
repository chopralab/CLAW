import sys
import os
import time
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.signal import find_peaks
import logging

class SampleIDExtract:
    def __init__(self, new_columns=None):
        if new_columns is None:
            new_columns = {
                'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
                'Genotype': ['5xFAD', 'WT'],
                'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
                'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
                'Other': ['Blank', 'blank']
            }
        self.new_columns = new_columns

    def extract_sample_parts(self, sample_id, std):
        parts = sample_id.replace('-', '_').split('_')
        matched_parts = {key: None for key in ['Biology', 'Genotype', 'Mouse', 'Cage', 'STD']}
        
        if std == 'yes':
            if 'FAME' in sample_id:
                matched_parts['STD'] = 'FAME'
            return 'FAME', matched_parts['STD']
        
        for part in parts:
            for key in matched_parts.keys():
                if matched_parts[key] is None:
                    for value in self.new_columns.get(key, []):
                        if value in part:
                            matched_parts[key] = value
                            break

        if matched_parts['Biology'] and matched_parts['Genotype'] and matched_parts['Mouse'] and matched_parts['Cage']:
            sample_name = '_'.join([matched_parts['Biology'], matched_parts['Genotype'], matched_parts['Mouse'], matched_parts['Cage']])
        elif matched_parts['STD']:
            sample_name = matched_parts['STD']
        else:
            sample_name = 'Unknown'
        # Debugging print statement
        logging.debug(f"Sample ID: {sample_id}, Extracted Sample Name: {sample_name}")

        std_name = matched_parts['STD'] if matched_parts['STD'] else 'None'
        return sample_name, std_name

    def find_std_rt_on(self, df, std, parent_ion, product_ion, tolerance):
        logging.info(f"Filtering DataFrame for standard based on Parent Ion: {parent_ion}, Product Ion: {product_ion}, with tolerance: {tolerance}")
        condition = (abs(df['Parent_Ion'] - parent_ion) <= tolerance) & (abs(df['Product_Ion'] - product_ion) <= tolerance)
        filtered_df = df[condition].copy()
        logging.info(f"Number of rows matching the ion filter condition: {len(filtered_df)}")
        df['STD_RT_ON'] = np.nan

        def get_highest_intensity_peak(group):
            peaks, _ = find_peaks(group['OzESI_Intensity'])
            if len(peaks) > 0:
                peak_idx = group.iloc[peaks]['OzESI_Intensity'].idxmax()
                return group.loc[peak_idx, 'Retention_Time']
            else:
                return np.nan

        rt_on_map = filtered_df.groupby('Sample').apply(get_highest_intensity_peak).to_dict()
        df['STD_RT_ON'] = df['Sample'].map(rt_on_map)
        return df

    def apply_extraction(self, df, std, parent_ion, product_ion, tolerance):
        tqdm.pandas(desc="Extracting Sample Parts")
        df[['Sample', 'STD']] = df['Sample_ID'].progress_apply(lambda x: self.extract_sample_parts(x, std)).apply(pd.Series)
        df = self.find_std_rt_on(df, std, parent_ion, product_ion, tolerance)
        return df

def main():
    logging.basicConfig(filename='sample_extraction_debug.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(console_handler)

    start_time = time.time()

    if len(sys.argv) > 1:
        std = sys.argv[1].strip().lower()
    else:
        logging.error("Error: Please provide STD as a command-line argument (yes/no)")
        sys.exit(1)

    new_columns = {
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
        'Other': ['Blank', 'blank']
    }

    output_dir = 'Projects/AMP/samples/ON/'
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Loading Parquet file...")
    OzESI_df = pd.read_parquet("Projects/AMP/mzml_parsed/ON/df_mzml_parser_1_AMP.parquet")
    
    logging.info("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(new_columns)

    # Defining parameters for standard ion extraction
    parent_ion = 425.40  # Mass-to-charge ratio of the parent ion for the standard
    product_ion = 183    # Mass-to-charge ratio of the product ion for the standard
    tolerance = 0.3      # Tolerance for matching ion values

    logging.info(f"Applying extraction with STD: {std} and finding STD_RT_ON...")
    OzON_Data = sample_extractor.apply_extraction(OzESI_df, std, parent_ion, product_ion, tolerance)
    logging.info("Extraction and STD_RT_ON calculation applied successfully.")

    unique_samples = OzON_Data['Sample'].unique()
    file_names = []
    for sample in tqdm(unique_samples, desc="Processing Samples"):
        sample_start_time = time.time()
        logging.info(f"Processing sample: {sample}")
        sample_df = OzON_Data[OzON_Data['Sample'] == sample]
        filename = os.path.join(output_dir, f"df_sample_2_{sample}_ON.parquet")
        logging.info(f"Saving {filename}...")
        sample_df.to_parquet(filename, index=False, compression="brotli")
        file_names.append(filename)
        sample_end_time = time.time()
        elapsed_time = sample_end_time - sample_start_time
        logging.info(f"File {filename} saved successfully. Processing took {elapsed_time:.2f} seconds.")

    end_time = time.time()
    total_time = end_time - start_time

    logging.info("Summary of Execution:")
    summary = (
        f"Number of unique Sample_ID values in input DataFrame: {len(OzESI_df['Sample_ID'].unique())}\n"
        f"Number of files created in output directory: {len(file_names)}\n"
        f"All output file names sorted: {sorted(file_names)}\n"
        f"Total script execution time: {total_time:.2f} seconds"
    )
    logging.info(summary)
    print(summary)

if __name__ == "__main__":
    main()
