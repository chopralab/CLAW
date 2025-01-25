import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import time
from scipy.signal import find_peaks
from scipy.integrate import trapz
from concurrent.futures import ProcessPoolExecutor
import argparse
import json

class SampleIDExtract:
    def __init__(self, new_columns=None):
        self.new_columns = new_columns or {}

    def extract_sample_parts(self, sample_id):
        # Split the Sample_ID into parts
        parts = sample_id.replace('-', '_').split('_')
        matched_parts = {key: None for key in ['Biology', 'Genotype', 'Mouse', 'Cage', 'STD']}
        
        # Check if "Blank" is in the Sample_ID
        if "Blank" in parts:
            sample_name = "Blank"  # Assign Sample as "Blank"
            std_name = None  # No STD for Blank
            return sample_name, std_name

        # Check if "CisTrans" is in the Sample_ID
        if "CisTrans" in parts:
            sample_name = "CT"  # Assign Sample as "CT"
            std_name = None  # No STD for CisTrans
            return sample_name, std_name

        # Match other parts based on the columns configuration
        for part in parts:
            for key in matched_parts.keys():
                if matched_parts[key] is None:
                    for value in self.new_columns.get(key, []):
                        if value in part:
                            matched_parts[key] = value
                            break

        # Determine the sample_name and std_name based on the matches
        if matched_parts['Biology'] and matched_parts['Genotype'] and matched_parts['Mouse'] and matched_parts['Cage']:
            sample_name = '_'.join([
                matched_parts['Biology'], 
                matched_parts['Genotype'], 
                matched_parts['Mouse'], 
                matched_parts['Cage']
            ])
        elif matched_parts['STD']:
            sample_name = matched_parts['STD']  # If only STD matches, use it as the sample name
        else:
            sample_name = 'Unknown'  # Default to "Unknown" if no matches

        std_name = matched_parts['STD'] if matched_parts['STD'] else 'None'

        return sample_name, std_name





    def get_highest_intensity_peak_and_area(self, group):
        peaks, _ = find_peaks(group['OzESI_Intensity'])
        if len(peaks) > 0:
            peak_idx = group.iloc[peaks]['OzESI_Intensity'].idxmax()
            rt_peak = group.loc[peak_idx, 'Retention_Time']
            peak_data = group.iloc[peaks]
            peak_area = trapz(peak_data['OzESI_Intensity'], peak_data['Retention_Time'])
            return rt_peak, group.loc[peak_idx, 'OzESI_Intensity'], peak_area
        else:
            return np.nan, np.nan, np.nan

    def find_std_rt_off(self, df, std, parent_ion, product_ion, tolerance):
        condition = (abs(df['Parent_Ion'] - parent_ion) <= tolerance) & \
                    (abs(df['Product_Ion'] - product_ion) <= tolerance)

        filtered_df = df[condition].copy()

        def apply_peak_info(group):
            rt_peak, peak_intensity, peak_area = self.get_highest_intensity_peak_and_area(group)
            return pd.Series([rt_peak, peak_intensity, peak_area], 
                           index=['STD_RT_OFF', 'STD_Peak_Intensity', 'STD_Peak_Area'])

        peak_info_map = filtered_df.groupby('Sample').apply(apply_peak_info)
        df['STD_RT_OFF'] = df['Sample'].map(lambda x: peak_info_map.loc[x, 'STD_RT_OFF'] 
                                           if x in peak_info_map.index else np.nan)
        df['STD_Peak_Intensity'] = df['Sample'].map(lambda x: peak_info_map.loc[x, 'STD_Peak_Intensity'] 
                                                   if x in peak_info_map.index else np.nan)
        df['STD_Peak_Area'] = df['Sample'].map(lambda x: peak_info_map.loc[x, 'STD_Peak_Area'] 
                                              if x in peak_info_map.index else np.nan)

        return df

    def apply_extraction(self, df, std, parent_ion, product_ion, tolerance):
        print("[DEBUG] Extracting sample parts for each row...")
        tqdm.pandas(desc="Extracting Sample Parts")
        
        # Apply `extract_sample_parts` to each row of the DataFrame
        extracted = df['Sample_ID'].progress_apply(self.extract_sample_parts)
        
        # Assign the extracted parts back to the DataFrame
        df[['Sample', 'STD']] = pd.DataFrame(extracted.tolist(), index=df.index)
        print("[DEBUG] Sample and STD columns added to DataFrame.")

        # Apply `find_std_rt_off` to the DataFrame
        print("[DEBUG] Finding STD_RT_OFF for the entire DataFrame...")
        df = self.find_std_rt_off(df, std, parent_ion, product_ion, tolerance)
        print("[DEBUG] STD_RT_OFF calculation completed.")
        
        return df

def process_sample(sample, OzON_Data, output_dir):
    print(f"[DEBUG] Processing sample: {sample}")
    sample_df = OzON_Data[OzON_Data['Sample'] == sample]
    filename = os.path.join(output_dir, f"df_sample_2_{sample}.parquet")
    sample_df.to_parquet(filename, index=False, compression="brotli")
    print(f"[DEBUG] Sample {sample} written to file: {filename}")
    return f"Sample {sample} processed."

def main():
    parser = argparse.ArgumentParser(description='Process OzESI data with configurable parameters')
    parser.add_argument('--input_file', required=True, help='Input parquet file path')
    parser.add_argument('--output_dir', required=True, help='Output directory path')
    parser.add_argument('--std', default='d2-16:0', help='Standard identifier')
    parser.add_argument('--parent_ion', type=float, default=425.40, help='Parent ion mass')
    parser.add_argument('--product_ion', type=float, default=183, help='Product ion mass')
    parser.add_argument('--tolerance', type=float, default=0.3, help='Mass tolerance')
    parser.add_argument('--max_workers', type=int, default=64, help='Maximum number of parallel workers')
    parser.add_argument('--columns_config', type=str, required=True, help='Path to columns configuration JSON file')
    
    args = parser.parse_args()
    start_time = time.time()

    # Load columns configuration from JSON file
    with open(args.columns_config, 'r') as f:
        new_columns = json.load(f)

    os.makedirs(args.output_dir, exist_ok=True)
    
    print("Loading Parquet file...")
    OzESI_df = pd.read_parquet(args.input_file)
    print("Parquet file loaded successfully.")
    
    sample_extractor = SampleIDExtract(new_columns)
    
    print("Applying extraction and finding STD_RT_OFF...")
    OzON_Data = sample_extractor.apply_extraction(
        OzESI_df, 
        args.std, 
        args.parent_ion, 
        args.product_ion, 
        args.tolerance
    )
    print("Extraction and STD_RT_OFF calculation applied successfully.")
    
    unique_samples = OzON_Data['Sample'].unique()
    
    print("Processing samples in parallel...")
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        results = list(tqdm(
            executor.map(
                process_sample, 
                unique_samples, 
                [OzON_Data] * len(unique_samples), 
                [args.output_dir] * len(unique_samples)
            ),
            total=len(unique_samples), 
            desc="Processing Samples"
        ))
        print("[DEBUG] Parallel sample processing complete.")
        for result in results:
            print(result)
    
    print(f"Total execution time: {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()