import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import time
from scipy.signal import find_peaks
from scipy.integrate import trapz
from concurrent.futures import ProcessPoolExecutor
import dask.dataframe as dd

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
            sample_name = 'Unknown'

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
            return pd.Series([rt_peak, peak_intensity, peak_area], index=['STD_RT_OFF', 'STD_Peak_Intensity', 'STD_Peak_Area'])

        peak_info_map = filtered_df.groupby('Sample').apply(apply_peak_info)
        df['STD_RT_OFF'] = df['Sample'].map(lambda x: peak_info_map.loc[x, 'STD_RT_OFF'] if x in peak_info_map.index else np.nan)
        df['STD_Peak_Intensity'] = df['Sample'].map(lambda x: peak_info_map.loc[x, 'STD_Peak_Intensity'] if x in peak_info_map.index else np.nan)
        df['STD_Peak_Area'] = df['Sample'].map(lambda x: peak_info_map.loc[x, 'STD_Peak_Area'] if x in peak_info_map.index else np.nan)

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
    start_time = time.time()  # Start tracking the execution time
    
    new_columns = {
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
        'Other': ['Blank', 'blank'],
        'STD': ['FAME']
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

    print("Processing samples in parallel...")
    with ProcessPoolExecutor(max_workers=64) as executor:
        results = list(tqdm(executor.map(process_sample, unique_samples, [OzON_Data] * len(unique_samples), [output_dir] * len(unique_samples)),
                            total=len(unique_samples), desc="Processing Samples"))
        print("[DEBUG] Parallel sample processing complete.")
        for result in results:
            print(result)

    # Log the total execution time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total execution time: {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    main()
