import pandas as pd
from tqdm import tqdm
import os
import numpy as np
from scipy.signal import find_peaks
import time
import sys

class SampleIDExtract:
    def __init__(self, sample_id_switch):
        self.sample_id_switch = sample_id_switch

    def extract_sample_id(self, df):
        if self.sample_id_switch == 'yes':
            # Set Sample column to the values in Sample_ID column
            df['Sample'] = df['Sample_ID']
        return df

    def find_std_rt_off(self, df, std, parent_ion, product_ion, tolerance):
        condition = (abs(df['Parent_Ion'] - parent_ion) <= tolerance) & \
                    (abs(df['Product_Ion'] - product_ion) <= tolerance)

        filtered_df = df[condition].copy()

        df['STD_RT_OFF'] = np.nan

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
        # Apply the Sample_ID extraction logic if requested
        df = self.extract_sample_id(df)
        df = self.find_std_rt_off(df, std, parent_ion, product_ion, tolerance)
        return df

def main():
    # Receive user input for SampleID (yes or no)
    sample_id_switch = sys.argv[1] if len(sys.argv) > 1 else 'no'

    output_dir = 'Projects/saniya/samples/OFF/'
    os.makedirs(output_dir, exist_ok=True)

    print("Loading Parquet file...")
    OzESI_df = pd.read_parquet("Projects/saniya/mzml_parsed/OFF/df_mzml_parser_1_OFF.parquet")
    print("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(sample_id_switch)

    std = 'd2-16:0'
    parent_ion = 425.40
    product_ion = 183
    tolerance = 0.3

    print(f"SampleID switch set to '{sample_id_switch}'.")
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
