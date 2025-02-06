#!/usr/bin/env python

import argparse
import logging
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from tqdm import tqdm


class SampleIDExtract:
    def __init__(self, new_columns):
        self.new_columns = new_columns

    def extract_sample_parts(self, sample_id):
        """
        Extracts sample name and STD name from the Sample_ID.

        Parameters:
            sample_id (str): The Sample_ID string.

        Returns:
            tuple: (sample_name, std_name)
        """
        # Split the Sample_ID into parts
        parts = sample_id.replace('-', '_').split('_')
        matched_parts = {key: None for key in ['Biology', 'Genotype', 'Mouse', 'Cage', 'STD']}

        # Convert all parts to lowercase for case-insensitive matching
        lower_parts = [part.lower() for part in parts]

        # Check if any part contains "blank"
        if any("blank" in part for part in lower_parts):
            sample_name = "Blank"  # Assign Sample as "Blank"
            std_name = None  # No STD for Blank
            return sample_name, std_name

        # Unified check for any form of 'cistrans'
        if any('cistrans' in part.lower() for part in parts):
            sample_name = "CT"
            std_name = None
            return sample_name, std_name


        # Match other parts based on the columns configuration
        for part in parts:
            lower_part = part.lower()
            for key in matched_parts.keys():
                if matched_parts[key] is None:
                    for value in self.new_columns.get(key, []):
                        if value.lower() in lower_part:  # Case-insensitive matching
                            matched_parts[key] = value
                            break

        # Determine the sample_name based on the matches
        if all(matched_parts[key] for key in ['Biology', 'Genotype', 'Mouse', 'Cage']):
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

    def calculate_peak_metrics(self, data, intensity_column, rt_column, peak_window=0.2):
        """
        Calculates peak metrics for a given dataset.

        Parameters:
            data (DataFrame): The input data.
            intensity_column (str): Column name for intensity.
            rt_column (str): Column name for retention time.
            peak_window (float): Window around the peak to calculate area.

        Returns:
            Series: Peak metrics.
        """
        data = data.sort_values(by=rt_column)
        peaks, properties = find_peaks(data[intensity_column], prominence=1)

        if peaks.size > 0:
            peak_idx = peaks[np.argmax(properties["prominences"])]
            peak_rt = data.iloc[peak_idx][rt_column]
            peak_intensity = data.iloc[peak_idx][intensity_column]

            area_window = data[(data[rt_column] >= peak_rt - peak_window) &
                               (data[rt_column] <= peak_rt + peak_window)]
            peak_area = np.trapz(area_window[intensity_column], area_window[rt_column])

            return pd.Series({
                'STD_RT_ON': peak_rt,
                'STD_Peak_Intensity': peak_intensity,
                'STD_Peak_Area': peak_area
            })
        else:
            return pd.Series({
                'STD_RT_ON': np.nan,
                'STD_Peak_Intensity': np.nan,
                'STD_Peak_Area': np.nan
            })

    def find_peak_and_area(self, df, parent_ion, product_ion, tolerance):
        """
        Filters the dataframe based on ion parameters and calculates peak metrics.

        Parameters:
            df (DataFrame): The input dataframe.
            parent_ion (float): Parent ion m/z ratio.
            product_ion (float): Product ion m/z ratio.
            tolerance (float): Tolerance for ion matching.

        Returns:
            DataFrame: Merged dataframe with peak metrics.
        """
        condition = (
            (df['Parent_Ion'].between(parent_ion - tolerance, parent_ion + tolerance)) &
            (df['Product_Ion'].between(product_ion - tolerance, product_ion + tolerance))
        )
        filtered_df = df[condition].copy()

        peak_data = filtered_df.groupby('Sample').apply(
            lambda group: self.calculate_peak_metrics(group, 'OzESI_Intensity', 'Retention_Time')
        ).reset_index()

        return df.merge(peak_data, on='Sample', how='left')

    def apply_extraction(self, df, std, parent_ion, product_ion, tolerance):
        """
        Applies sample extraction and peak metric calculations.

        Parameters:
            df (DataFrame): The input dataframe.
            std (str): STD usage flag.
            parent_ion (float): Parent ion m/z ratio.
            product_ion (float): Product ion m/z ratio.
            tolerance (float): Tolerance for ion matching.

        Returns:
            DataFrame: Processed dataframe.
        """
        df[['Sample', 'STD']] = df['Sample_ID'].apply(
            self.extract_sample_parts
        ).apply(pd.Series)
        df = self.find_peak_and_area(df, parent_ion, product_ion, tolerance)
        return df


def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
        Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Extract sample information and calculate peak metrics.")
    parser.add_argument('--std', choices=['yes', 'no'], required=True, help="Specify if STD is used ('yes' or 'no').")
    parser.add_argument('--input_parquet', required=True, help="Path to the input Parquet file.")
    parser.add_argument('--output_dir', required=True, help="Directory to save the output Parquet files.")
    parser.add_argument('--parent_ion', type=float, default=425.40, help="Parent ion m/z ratio.")
    parser.add_argument('--product_ion', type=float, default=183, help="Product ion m/z ratio.")
    parser.add_argument('--tolerance', type=float, default=0.3, help="Tolerance for ion matching.")
    return parser.parse_args()


def setup_logging():
    """
    Sets up logging configuration.
    """
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def main():
    args = parse_arguments()
    setup_logging()

    start_time = time.time()

    # Define new_columns
    new_columns = {
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
        'STD': ['STD1', 'STD2'],  # Example STD values; update as needed
        # 'Other': ['Blank', 'blank']  # Removed if handled separately
    }

    os.makedirs(args.output_dir, exist_ok=True)

    logging.info("Loading Parquet file...")
    try:
        OzESI_df = pd.read_parquet(args.input_parquet)
        logging.info("Parquet file loaded successfully.")
    except Exception as e:
        logging.error(f"Failed to load Parquet file: {e}")
        sys.exit(1)

    sample_extractor = SampleIDExtract(new_columns)

    logging.info(f"Applying extraction with STD: {args.std} and finding STD_RT_ON...")

    # Apply extraction and calculate peak metrics
    OzON_Data = sample_extractor.apply_extraction(
        OzESI_df,
        args.std,
        args.parent_ion,
        args.product_ion,
        args.tolerance
    )
    logging.info("Extraction and STD_RT_ON calculation applied successfully.")

    unique_samples = OzON_Data['Sample'].unique()
    file_names = []

    for sample in tqdm(unique_samples, desc="Processing Samples"):
        sample_start_time = time.time()
        logging.info(f"Processing sample: {sample}")
        sample_df = OzON_Data[OzON_Data['Sample'] == sample]
        filename = os.path.join(args.output_dir, f"df_sample_2_{sample}_ON.parquet")
        try:
            sample_df.to_parquet(filename, index=False, compression="brotli")
            file_names.append(filename)
            logging.info(f"Saved {filename}.")
        except Exception as e:
            logging.error(f"Failed to save {filename}: {e}")
        sample_end_time = time.time()
        logging.info(f"Processing took {sample_end_time - sample_start_time:.2f} seconds.")

    total_time = time.time() - start_time

    summary = (
        f"Number of unique Sample_ID values in input DataFrame: {OzESI_df['Sample_ID'].nunique()}\n"
        f"Number of files created in output directory: {len(file_names)}\n"
        f"All output file names sorted: {sorted(file_names)}\n"
        f"Total script execution time: {total_time:.2f} seconds"
    )
    logging.info("Summary of Execution:")
    logging.info(summary)
    print(summary)


if __name__ == "__main__":
    main()
