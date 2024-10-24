import pandas as pd
from tqdm import tqdm
import os
import numpy as np
from scipy.signal import find_peaks
import time

class SampleIDExtract:
    def __init__(self, new_columns=None):
        if new_columns is None:
            new_columns = {
                'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
                'Genotype': ['5xFAD', 'WT'],
                'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
                'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
                'Other': ['Blank', 'blank'],
                'STD': ['FAME']  # Added standard identifier
            }
        self.new_columns = new_columns

    def extract_sample_parts(self, sample_id):
        parts = sample_id.replace('-', '_').split('_')
        matched_parts = {key: None for key in ['Biology', 'Genotype', 'Mouse', 'Cage', 'STD']}

        for part in parts:
            for key in matched_parts.keys():
                if matched_parts[key] is None:
                    for value in self.new_columns.get(key, []):
                        if value in part:
                            matched_parts[key] = value
                            break

        # Construct the Sample name in the desired order
        if matched_parts['Biology'] and matched_parts['Genotype'] and matched_parts['Mouse'] and matched_parts['Cage']:
            sample_name = '_'.join([
                matched_parts['Biology'], 
                matched_parts['Genotype'], 
                matched_parts['Mouse'], 
                matched_parts['Cage']
            ])
        elif matched_parts['STD']:
            sample_name = matched_parts['STD']
        else:
            sample_name = 'Unknown'  # Fallback if neither criteria are met

        # Determine if 'FAME' or other STD is present in the sample ID
        std_name = matched_parts['STD'] if matched_parts['STD'] else 'None'

        return sample_name, std_name

    def find_std_rt_off(self, df, std, parent_ion, product_ion, tolerance):
        condition = (abs(df['Parent_Ion'] - parent_ion) <= tolerance) & \
                    (abs(df['Product_Ion'] - product_ion) <= tolerance)

        filtered_df = df[condition].copy()

        df['STD_RT_OFF'] = np.nan  # Change from 'STD_RT_ON' to 'STD_RT_OFF'

        def get_highest_intensity_peak(group):
            peaks, _ = find_peaks(group['OzESI_Intensity'])

            if len(peaks) > 0:
                peak_idx = group.iloc[peaks]['OzESI_Intensity'].idxmax()
                return group.loc[peak_idx, 'Retention_Time']
            else:
                return np.nan

        rt_off_map = filtered_df.groupby('Sample').apply(get_highest_intensity_peak).to_dict()

        df['STD_RT_OFF'] = df['Sample'].map(rt_off_map)

        return df

    def apply_extraction(self, df, std, parent_ion, product_ion, tolerance):
        tqdm.pandas(desc="Extracting Sample Parts")
        extracted = df['Sample_ID'].progress_apply(self.extract_sample_parts)
        df[['Sample', 'STD']] = pd.DataFrame(extracted.tolist(), index=df.index)
        df = self.find_std_rt_off(df, std, parent_ion, product_ion, tolerance)  # Change method call to 'find_std_rt_off'
        return df

def main():
    new_columns = {
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
        'Other': ['Blank', 'blank'],
        'STD': ['FAME']  # Standard identifier added here
    }

    output_dir = 'Projects/AMP/samples/OFF/'
    os.makedirs(output_dir, exist_ok=True)

    print("Loading Parquet file...")
    OzESI_df = pd.read_parquet("Projects/AMP/mzml_parsed/OFF/df_mzml_parser_1_OFF.parquet")
    print("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(new_columns)

    std = 'd2-16:0'
    parent_ion = 425.40
    product_ion = 183
    tolerance = 0.3

    print("Applying extraction and finding STD_RT_OFF...")
    OzON_Data = sample_extractor.apply_extraction(OzESI_df, std, parent_ion, product_ion, tolerance)
    print("Extraction and STD_RT_OFF calculation applied successfully.")

    unique_samples = OzON_Data['Sample'].unique()
    for sample in tqdm(unique_samples, desc="Processing Samples"):
        start_time = time.time()
        
        print(f"Processing sample: {sample}")
        sample_df = OzON_Data[OzON_Data['Sample'] == sample]
        filename = os.path.join(output_dir, f"df_sample_2_{sample}.parquet")
        
        print(f"Saving {filename}...")
        sample_df.to_parquet(filename, index=False, compression="brotli")
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        print(f"File {filename} saved successfully.")
        print(f"Processing and saving took {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
