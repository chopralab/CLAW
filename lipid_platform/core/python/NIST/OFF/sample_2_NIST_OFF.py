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
        1) If 'blank' appears in sample_id (case-insensitive), return ('Blank', None, None).
        2) Otherwise:
           - Force STD='NIST'
           - Parse replicate from sample_id (n1, n2, n3, etc.) if found
           - Otherwise replicate='n0'
           - Then Sample='NIST_{replicate}'
        """
        parts = sample_id.replace('-', '_').split('_')
        lower_parts = [p.lower() for p in parts]

        # 1) Detect if 'blank'
        if any("blank" in lp for lp in lower_parts):
            # Return 'Blank' sample, no STD, no replicate
            return "Blank", None, None

        # 2) Not a blank → parse replicate
        matched_replicate = None
        for p in parts:
            lower_p = p.lower()
            for val in self.new_columns.get('Replicate', []):
                if val.lower() in lower_p:
                    matched_replicate = val
                    break
            if matched_replicate:
                break

        if matched_replicate is None:
            matched_replicate = "n0"

        # Hardcode
        sample_name = f"NIST_{matched_replicate}"
        std_name = "NIST"

        return sample_name, std_name, matched_replicate

    def calculate_peak_metrics(self, data, intensity_column, rt_column, peak_window=0.2):
        """
        Calculates peak metrics for a given dataset.
        """
        data = data.sort_values(by=rt_column)
        peaks, properties = find_peaks(data[intensity_column], prominence=1)

        if peaks.size > 0:
            peak_idx = peaks[np.argmax(properties["prominences"])]
            peak_rt = data.iloc[peak_idx][rt_column]
            peak_intensity = data.iloc[peak_idx][intensity_column]

            area_window = data[
                (data[rt_column] >= peak_rt - peak_window) &
                (data[rt_column] <= peak_rt + peak_window)
            ]
            peak_area = np.trapz(area_window[intensity_column], area_window[rt_column])

            return pd.Series({
                'STD_RT_OFF': peak_rt,
                'STD_Peak_Intensity': peak_intensity,
                'STD_Peak_Area': peak_area
            })
        else:
            return pd.Series({
                'STD_RT_OFF': np.nan,
                'STD_Peak_Intensity': np.nan,
                'STD_Peak_Area': np.nan
            })

    def find_peak_and_area(self, df, parent_ion, product_ion, tolerance):
        """
        Filters the dataframe based on ion parameters and calculates peak metrics.
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
        1) Extract sample info (including replicate).
        2) Handle blanks so that each blank Sample_ID stays at replicate = n0.
        3) Find peak/area data.
        """
        # Extract (Sample, STD, Replicate) from Sample_ID
        df[['Sample', 'STD', 'Replicate']] = df['Sample_ID'].apply(
            self.extract_sample_parts
        ).apply(pd.Series)

        # --- UPDATED: Assign replicate='n0' for all blank rows at once ---
        blank_mask = (df['Sample'] == 'Blank')
        df.loc[blank_mask, 'Replicate'] = 'n0'
        df.loc[blank_mask, 'STD'] = None  # if you want STD cleared for blanks

        # Now do peak-finding
        df = self.find_peak_and_area(df, parent_ion, product_ion, tolerance)
        return df


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract sample info and calculate peak metrics.")
    parser.add_argument('--std', choices=['yes', 'no'], required=True,
                        help="Specify if STD is used ('yes' or 'no').")
    parser.add_argument('--input_parquet', required=True,
                        help="Path to the input Parquet file.")
    parser.add_argument('--output_dir', required=True,
                        help="Directory to save the output Parquet files.")
    parser.add_argument('--parent_ion', type=float, default=425.40,
                        help="Parent ion m/z ratio.")
    parser.add_argument('--product_ion', type=float, default=183,
                        help="Product ion m/z ratio.")
    parser.add_argument('--tolerance', type=float, default=0.3,
                        help="Tolerance for ion matching.")
    return parser.parse_args()


def setup_logging():
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

    # We only need replicates here
    new_columns = {
        'Replicate': ['n1', 'n2', 'n3', 'n4', 'n5', 'n6']
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

    logging.info(f"Applying extraction with STD: {args.std} and finding STD_RT_OFF...")

    # Run extraction & peak metrics
    OzON_Data = sample_extractor.apply_extraction(
        OzESI_df,
        args.std,
        args.parent_ion,
        args.product_ion,
        args.tolerance
    )
    logging.info("Extraction and STD_RT_OFF calculation applied successfully.")

    # 1) Separate blanks vs. non-blanks
    blank_df = OzON_Data[OzON_Data['Sample'] == 'Blank'].copy()
    non_blanks = OzON_Data[OzON_Data['Sample'] != 'Blank'].copy()

    # 2) Save non-blanks by grouping (Sample, Replicate)
    group_iter = non_blanks.groupby(['Sample', 'Replicate'])
    file_names = []

    for (sample, replicate), sample_df in tqdm(group_iter, desc="Processing (Sample, Replicate)"):
        sample_start_time = time.time()
        logging.info(f"Processing sample: {sample}, replicate: {replicate}")

        filename = os.path.join(
            args.output_dir,
            f"df_sample_2_{sample}_ON.parquet"
        )

        try:
            sample_df.to_parquet(filename, index=False, compression="brotli")
            file_names.append(filename)
            logging.info(f"Saved {filename}.")
        except Exception as e:
            logging.error(f"Failed to save {filename}: {e}")

        sample_end_time = time.time()
        logging.info(f"Processing took {sample_end_time - sample_start_time:.2f} seconds.")

    # 3) Save all blanks in ONE file
    if not blank_df.empty:
        blank_file = os.path.join(args.output_dir, "df_sample_2_Blank_ALL_ON.parquet")
        try:
            blank_df.to_parquet(blank_file, index=False, compression="brotli")
            file_names.append(blank_file)
            logging.info(f"Saved all Blanks to {blank_file}")
        except Exception as e:
            logging.error(f"Failed to save {blank_file}: {e}")

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
