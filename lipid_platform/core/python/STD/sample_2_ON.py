from scipy.signal import find_peaks
import pandas as pd
import numpy as np
import os
import sys
import time
from tqdm import tqdm
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
        
        std_name = matched_parts['STD'] if matched_parts['STD'] else 'None'
        return sample_name, std_name

    def calculate_peak_metrics(self, data, intensity_column, rt_column, peak_window=0.2):
        """
        Identify peaks and calculate peak retention time, intensity, and area.
        Ensures the data is sorted by Retention_Time before analysis.
        """
        # Ensure the data is sorted by Retention_Time
        data = data.sort_values(by=rt_column)

        # Identify peaks in the intensity data
        peaks, properties = find_peaks(data[intensity_column], prominence=1)  # Adjust 'prominence' as needed
        if len(peaks) > 0:
            # Find the highest peak
            peak_idx = data.iloc[peaks][intensity_column].idxmax()
            peak_rt = data.loc[peak_idx, rt_column]
            peak_intensity = data.loc[peak_idx, intensity_column]

            # Calculate area around the peak within the specified window
            area_window = data[(data[rt_column] >= peak_rt - peak_window) &
                               (data[rt_column] <= peak_rt + peak_window)]
            peak_area = np.trapz(area_window[intensity_column], area_window[rt_column])

            return pd.Series({'STD_RT_ON': peak_rt, 'STD_Peak_Intensity': peak_intensity, 'STD_Peak_Area': peak_area})
        else:
            return pd.Series({'STD_RT_ON': np.nan, 'STD_Peak_Intensity': np.nan, 'STD_Peak_Area': np.nan})

    def find_peak_and_area(self, df, parent_ion, product_ion, tolerance):
        """
        Find peak retention time, intensity, and area for each sample.
        Ensures each group is sorted by Retention_Time before calculating metrics.
        """
        # Filter for the standard
        condition = (abs(df['Parent_Ion'] - parent_ion) <= tolerance) & (abs(df['Product_Ion'] - product_ion) <= tolerance)
        filtered_df = df[condition].copy()

        # Group by Sample and calculate peak metrics
        peak_data = filtered_df.groupby('Sample').apply(
            lambda group: self.calculate_peak_metrics(group, 'OzESI_Intensity', 'Retention_Time')
        ).reset_index()
        # Merge back with the original DataFrame
        df = df.merge(peak_data, on='Sample', how='left')
        return df

    def apply_extraction(self, df, std, parent_ion, product_ion, tolerance):
        """
        Extract and calculate peaks and areas, adding new columns to the DataFrame.
        """
        df[['Sample', 'STD']] = df['Sample_ID'].apply(lambda x: self.extract_sample_parts(x, std)).apply(pd.Series)
        df = self.find_peak_and_area(df, parent_ion, product_ion, tolerance)
        return df

    @staticmethod
    def flat_baseline_per_parent_ion(df, window_start, window_end, filter_column='Retention_Time', average_column='OzESI_Intensity', group_column='Parent_Ion'):
        logging.info("Starting flat_baseline_per_parent_ion calculation.")
        logging.debug("DataFrame columns: {}".format(df.columns.tolist()))
        logging.debug("Window: {} to {}, Filter column: '{}', Average column: '{}', Group column: '{}'".format(
            window_start, window_end, filter_column, average_column, group_column
        ))

        # Define a function to calculate the baseline for each group
        def calculate_group_baseline(group):
            group_name = group[group_column].iloc[0]
            logging.debug("Processing group '{}', number of rows: {}".format(group_name, len(group)))
            
            df_window = group[(group[filter_column] >= window_start) & (group[filter_column] <= window_end)]
            logging.debug("Group '{}': Number of rows in window: {}".format(group_name, len(df_window)))

            if df_window.empty:
                average_value = np.nan
                logging.warning("Group '{}': No data in window {} to {} for baseline calculation. Setting baseline to NaN.".format(group_name, window_start, window_end))
            else:
                average_value = df_window[average_column].mean()
                logging.debug("Group '{}': Calculated baseline value: {}".format(group_name, average_value))

            group['flat_baseline'] = average_value
            return group

        # Apply the function to each group
        df = df.groupby(group_column).apply(calculate_group_baseline).reset_index(drop=True)
        logging.info("Completed flat_baseline_per_parent_ion calculation.")
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

    output_dir = 'Projects/STD/samples/ON/'
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Loading Parquet file...")
    OzESI_df = pd.read_parquet("Projects/STD/mzml_parsed/ON/df_mzml_parser_1_FAME.parquet")
    
    logging.info("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(new_columns)

    # Defining parameters for standard ion extraction
    parent_ion = 425.40  # Mass-to-charge ratio of the parent ion for the standard
    product_ion = 183    # Mass-to-charge ratio of the product ion for the standard
    tolerance = 0.3      # Tolerance for matching ion values

    logging.info(f"Applying extraction with STD: {std} and finding STD_RT_ON...")

    # Apply extraction and add the new columns for peak metrics
    OzON_Data = sample_extractor.apply_extraction(OzESI_df, std, parent_ion, product_ion, tolerance)
    logging.info("Extraction and STD_RT_ON calculation applied successfully.")

    # Define retention time window for flat baseline calculation
    window_start = 24  # Adjust as needed based on your data
    window_end = 25    # Adjust as needed based on your data

    logging.info("Applying flat_baseline_per_parent_ion function.")
    OzON_Data = sample_extractor.flat_baseline_per_parent_ion(
        OzON_Data,
        window_start=window_start,
        window_end=window_end,
        filter_column='Retention_Time',
        average_column='OzESI_Intensity',
        group_column='Parent_Ion'
    )
    logging.info("flat_baseline_per_parent_ion applied successfully.")

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
