import sys
import os
import time
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.signal import find_peaks

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
        print(f"Sample ID: {sample_id}, Extracted Sample Name: {sample_name}")

        std_name = matched_parts['STD'] if matched_parts['STD'] else 'None'
        return sample_name, std_name

    def find_std_rt_on(self, df, std, parent_ion, product_ion, tolerance):
        condition = (abs(df['Parent_Ion'] - parent_ion) <= tolerance) & (abs(df['Product_Ion'] - product_ion) <= tolerance)
        filtered_df = df[condition].copy()
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
    if len(sys.argv) > 1:
        std = sys.argv[1].strip().lower()
    else:
        print("Error: Please provide STD as a command-line argument (yes/no)")
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

    print("Loading Parquet file...")
    OzESI_df = pd.read_parquet("Projects/AMP/mzml_parsed/ON/df_mzml_parser_1_AMP.parquet")
    
    print("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(new_columns)

    parent_ion = 425.40
    product_ion = 183
    tolerance = 0.3

    print(f"Applying extraction with STD: {std} and finding STD_RT_ON...")
    OzON_Data = sample_extractor.apply_extraction(OzESI_df, std, parent_ion, product_ion, tolerance)
    print("Extraction and STD_RT_ON calculation applied successfully.")

    unique_samples = OzON_Data['Sample'].unique()
    for sample in tqdm(unique_samples, desc="Processing Samples"):
        start_time = time.time()
        print(f"Processing sample: {sample}")
        sample_df = OzON_Data[OzON_Data['Sample'] == sample]
        filename = os.path.join(output_dir, f"df_sample_2_{sample}_ON.parquet")
        print(f"Saving {filename}...")
        sample_df.to_parquet(filename, index=False, compression="brotli")
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"File {filename} saved successfully. Processing took {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
