#!/usr/bin/env python3

import argparse
import os
import sys
import time
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.integrate import trapz
from tqdm import tqdm

class LipidAnalysis:
    def __init__(self, data, height=1000, width=None, rel_height=0.5):
        """
        Initialize the LipidAnalysis class.

        :param data: DataFrame containing lipid analysis data.
        :param height: Minimum height of peaks.
        :param width: Minimum width of peaks.
        :param rel_height: Relative height for peak width calculation.
        """
        self.data = data
        self.height = height
        self.width = width
        self.rel_height = rel_height

    def extract_species_info(self, species):
        """
        Extract carbon number and double bond information from species string.

        :param species: String containing species information.
        :return: Tuple containing carbon number and double bond count for sorting purposes.
        """
        parts = species.split(':')
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            carbon_number = int(parts[0])
            double_bond = int(parts[1])
        else:
            carbon_number = float('inf')  # Use a large number for unknown formats
            double_bond = float('inf')
        return carbon_number, double_bond

    def find_lipid_peaks(self, ignore_columns=False):
        """
        Find peaks in lipid data.

        :param ignore_columns: Boolean to determine whether to drop 'Biology', 'Genotype', 'Cage', 'Mouse' columns.
        :return: DataFrame containing peak data.
        """
        peak_data = []

        filter_col = 'group_by_lipid'
        unique_groups = self.data[filter_col].unique()

        for group in tqdm(unique_groups, desc=f"Processing {filter_col} groups"):
            group_data = self.data[self.data[filter_col] == group]
            peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=self.height, width=self.width)
            if len(peaks) == 0:
                continue  # Skip if no peaks found
            results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=self.rel_height)

            retention_times = group_data['Retention_Time'].values
            if len(retention_times) > 1:
                sampling_interval = retention_times[1] - retention_times[0]
            else:
                sampling_interval = 1  # Fallback value in case there's only one retention time

            for i, peak in enumerate(peaks):
                # Extract metadata
                metadata = group_data.iloc[peak][[
                    'Parent_Ion', 'Product_Ion', 'Sample', 'Species', 
                    'group_by_lipid', 'group_by_ion', 'Lipid', 'STD', 
                    'STD_RT_OFF', 'STD_Peak_Intensity', 'STD_Peak_Area',
                    'Sample_ID', 'Transition', 'Class'
                ]].to_dict()

                # Include additional columns if not ignored
                if not ignore_columns:
                    metadata.update({
                        'Biology': group_data.iloc[peak].get('Biology', None),
                        'Genotype': group_data.iloc[peak].get('Genotype', None),
                        'Cage': group_data.iloc[peak].get('Cage', None),
                        'Mouse': group_data.iloc[peak].get('Mouse', None)
                    })

                left_ip = results_half[2][i]
                right_ip = results_half[3][i]
                try:
                    left_time = group_data['Retention_Time'].iloc[int(left_ip)]
                    right_time = group_data['Retention_Time'].iloc[int(right_ip)]
                except IndexError:
                    # Handle cases where left_ip or right_ip are out of bounds
                    left_time = group_data['Retention_Time'].iloc[0]
                    right_time = group_data['Retention_Time'].iloc[-1]
                width_in_time = right_time - left_time

                # Calculate FWHM
                fwhm = results_half[0][i] * sampling_interval

                # Calculate peak area using trapezoidal integration
                peak_indices = np.arange(int(left_ip), int(right_ip) + 1)
                # Ensure indices are within the data range
                peak_indices = peak_indices[(peak_indices >= 0) & (peak_indices < len(group_data))]
                intensity_values = group_data['OzESI_Intensity'].iloc[peak_indices].values
                time_values = group_data['Retention_Time'].iloc[peak_indices].values
                peak_area = trapz(intensity_values, time_values)  # Integrate intensity over time

                peak_data.append({
                    'Lipid': metadata['Lipid'],
                    'Retention_Time': group_data.iloc[peak]['Retention_Time'],
                    'OzESI_Intensity': group_data.iloc[peak]['OzESI_Intensity'],
                    'group_by_ion': metadata['group_by_ion'],
                    'group_by_lipid': metadata['group_by_lipid'],
                    'Sample_ID': metadata['Sample_ID'],
                    'Transition': metadata['Transition'],
                    'Sample': metadata['Sample'],
                    'Parent_Ion': metadata['Parent_Ion'],
                    'Product_Ion': metadata['Product_Ion'],
                    'Species': metadata['Species'],
                    'Class': metadata['Class'],
                    'STD': metadata['STD'],
                    'STD_RT_OFF': metadata['STD_RT_OFF'],
                    'STD_Peak_Intensity': metadata['STD_Peak_Intensity'],
                    'STD_Peak_Area': metadata['STD_Peak_Area'],
                    'Peak_Height': properties['peak_heights'][i],
                    'FWHM': fwhm,
                    'Peak_Width': width_in_time,
                    'Peak_Area': peak_area,
                    'Filter_Column': filter_col
                })

        peaks_df = pd.DataFrame(peak_data)

        if peaks_df.empty:
            print("No peaks were found with the given parameters.")
            return peaks_df

        # Sort the DataFrame by species information
        peaks_df['species_sort'] = peaks_df['Species'].apply(self.extract_species_info)
        peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

        return peaks_df

    def isomer_selection(self, peaks_df):
        """
        Identify cis/trans isomers for lipid species with double bonds (e.g., 18:1 or 16:1).

        :param peaks_df: DataFrame with peak information.
        :return: Updated DataFrame with a new column 'Isomer' containing 'cis' or 'trans' values.
        """
        # Filter for lipids with species ending in ':1' (e.g., 18:1, 16:1)
        isomer_lipids = peaks_df[peaks_df['Species'].str.endswith(':1')]

        # Initialize a column for Isomer labeling
        peaks_df['Isomer'] = None

        for lipid in isomer_lipids['Species'].unique():
            lipid_data = peaks_df[peaks_df['Species'] == lipid]

            # Find the two largest peaks based on intensity
            if len(lipid_data) >= 2:
                # Sort by OzESI_Intensity to get the two largest peaks
                top_peaks = lipid_data.nlargest(2, 'OzESI_Intensity').sort_values(by='Retention_Time')

                # Label the first peak (based on retention time) as cis and the second as trans
                peaks_df.loc[top_peaks.index[0], 'Isomer'] = 'cis'
                peaks_df.loc[top_peaks.index[1], 'Isomer'] = 'trans'

        return peaks_df

def parse_arguments():
    parser = argparse.ArgumentParser(description="Lipid Peak Analysis Script")
    parser.add_argument("input_file", type=str, help="Path to the input parquet file.")
    parser.add_argument("height", type=float, help="Minimum height of peaks.")
    parser.add_argument("width", type=float, help="Minimum width of peaks.")
    parser.add_argument("rel_height", type=float, help="Relative height for peak width calculation.")
    parser.add_argument("ignore_columns_flag", type=str, choices=["ignore", "keep"], help="Flag to ignore specific columns ('ignore' or 'keep').")
    parser.add_argument("--output_dir", type=str, default=None, help="Directory to save the output file. Defaults to the input file's directory if not specified.")
    parser.add_argument("--max_peaks", action='store_true', help="Flag to calculate and save maximum peaks.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Validate input file
    if not os.path.isfile(args.input_file):
        print(f"Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Determine output directory
    output_dir = args.output_dir if args.output_dir else os.path.dirname(args.input_file)
    os.makedirs(output_dir, exist_ok=True)

    # Load the input DataFrame
    print("Loading Parquet file...")
    df_grouped = pd.read_parquet(args.input_file)
    print("Parquet file loaded successfully.")

    # Determine whether to ignore specific columns
    ignore_columns = args.ignore_columns_flag.lower() == "ignore"

    # Initialize LipidAnalysis class and set parameters
    analysis = LipidAnalysis(
        data=df_grouped,
        height=args.height,
        width=args.width,
        rel_height=args.rel_height
    )

    # Apply peak extraction
    print("Applying peak extraction...")
    peaks_df = analysis.find_lipid_peaks(ignore_columns=ignore_columns)
    print("Peak extraction applied successfully.")

    if peaks_df.empty:
        print("No peaks to process. Exiting.")
        sys.exit(0)

    # Perform isomer selection
    print("Applying isomer selection...")
    isomer_df = analysis.isomer_selection(peaks_df)
    print("Isomer selection applied successfully.")

    # Construct output file name
    sample_id = df_grouped['Sample'].iloc[0] if 'Sample' in df_grouped.columns else "unknown_sample"
    output_filename = f"df_analysis_CT_OFF_{sample_id}.parquet"
    output_file = os.path.join(output_dir, output_filename)

    # Save the results
    isomer_df.to_parquet(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")

