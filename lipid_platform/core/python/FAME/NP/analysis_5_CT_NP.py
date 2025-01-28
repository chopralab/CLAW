#!/usr/bin/env python3

import os
import sys
import time
import argparse
import pandas as pd
from tqdm import tqdm
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt

class LipidAnalysis:
    def __init__(self, data):
        """
        Initialize the LipidAnalysis class.

        :param data: DataFrame containing lipid analysis data.
        """
        self.data = data
        self.height = 1000
        self.width = None
        self.rel_height = 0.5

    def set_parameters(self, height=1000, width=None, rel_height=0.5):
        """
        Set parameters for peak finding.

        :param height: Minimum height of peaks.
        :param width: Minimum width of peaks.
        :param rel_height: Relative height for peak width calculation.
        """
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

    def isomer_selection(self, peaks_df):
        """
        Identify cis/trans isomers for lipid species with double bonds (e.g., 18:1 or 16:1).
        The function finds the two largest peaks by intensity for species with a double bond of 1
        and labels the first peak (based on retention time) as 'cis' and the second as 'trans'.

        :param peaks_df: DataFrame with peak information, including Peak_Height, FWHM, etc.
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

        # Return the DataFrame with isomer information and peak data
        return peaks_df

    def find_lipid_peaks(self, max_peaks=False, ignore_columns=False):
        """
        Find peaks in lipid data.

        :param max_peaks: Boolean to determine if maximum peaks should be calculated.
        :param ignore_columns: Boolean to determine whether to drop 'Biology', 'Genotype', 'Cage', 'Mouse' columns.
        :return: DataFrame containing peak data.
        """
        peak_data = []

        for filter_col in ['group_by_lipid']:
            unique_groups = self.data[filter_col].unique()

            for group in tqdm(unique_groups, desc=f"Processing {filter_col} groups"):
                group_data = self.data[self.data[filter_col] == group]
                peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=self.height, width=self.width)
                results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=self.rel_height)

                retention_times = group_data['Retention_Time'].values
                if len(retention_times) > 1:
                    sampling_interval = retention_times[1] - retention_times[0]
                else:
                    sampling_interval = 1  # Fallback value in case there's only one retention time

                for i, peak in enumerate(peaks):
                    metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 
                                  'group_by_lipid', 'group_by_ion', 'Lipid', 
                                  'STD', 'STD_RT_OFF']].to_dict()

                    # Ignore specific columns if flagged
                    if not ignore_columns:
                        metadata.update({
                            'Biology': group_data.iloc[peak]['Biology'],
                            'Genotype': group_data.iloc[peak]['Genotype'],
                            'Cage': group_data.iloc[peak]['Cage'],
                            'Mouse': group_data.iloc[peak]['Mouse']
                        })

                    left_ip = results_half[2][i]
                    right_ip = results_half[3][i]
                    left_time = group_data['Retention_Time'].iloc[int(left_ip)]
                    right_time = group_data['Retention_Time'].iloc[int(right_ip)]
                    width_in_time = right_time - left_time

                    fwhm = results_half[0][i] * sampling_interval

                    peak_data.append({
                        'Lipid': metadata['Lipid'],
                        'Retention_Time': group_data.iloc[peak]['Retention_Time'],
                        'OzESI_Intensity': group_data.iloc[peak]['OzESI_Intensity'],
                        'group_by_ion': metadata['group_by_ion'],
                        'group_by_lipid': metadata['group_by_lipid'],
                        'Sample_ID': group_data.iloc[peak]['Sample_ID'],
                        'Transition': group_data.iloc[peak]['Transition'],
                        'Sample': metadata['Sample'],
                        'Parent_Ion': metadata['Parent_Ion'],
                        'Product_Ion': metadata['Product_Ion'],
                        'Species': metadata['Species'],
                        'Class': group_data.iloc[peak]['Class'],
                        'STD': metadata['STD'],  # Ensure 'STD' is included
                        'STD_RT_OFF': metadata['STD_RT_OFF'],  # Ensure 'STD_RT_OFF' is included
                        'Peak_Height': properties['peak_heights'][i],
                        'FWHM': fwhm,
                        'Peak_Width': width_in_time,
                        'Peak_Area': properties['peak_heights'][i] * width_in_time,
                        'Filter_Column': filter_col  # Track which column was used for filtering
                    })

        peaks_df = pd.DataFrame(peak_data)

        # Sort the DataFrame by species information and return
        peaks_df['species_sort'] = peaks_df['Species'].apply(self.extract_species_info)
        peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

        if max_peaks:
            peaks_df = self.create_max_peaks_df(peaks_df)
            return peaks_df
        else:
            return peaks_df

    def create_max_peaks_df(self, peaks_df):
        """
        Create a DataFrame containing the maximum peaks per group.

        :param peaks_df: DataFrame containing all peaks.
        :return: DataFrame with maximum peaks.
        """
        # Placeholder for actual implementation
        # You need to define how to select max peaks
        return peaks_df

    def add_suffix_to_filename(self, filename, suffix):
        """
        Add a suffix to the filename before the file extension.

        :param filename: Original filename.
        :param suffix: Suffix to add.
        :return: Modified filename with suffix.
        """
        base, ext = os.path.splitext(filename)
        return f"{base}{suffix}{ext}"

def parse_arguments():
    parser = argparse.ArgumentParser(description="Lipid Analysis Peak Detection Script")

    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the input Parquet file containing lipid data."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory where the output Parquet file will be saved."
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1000,
        help="Minimum height of peaks. Default is 1000."
    )
    parser.add_argument(
        "--width",
        type=float,
        default=None,
        help="Minimum width of peaks. Default is None."
    )
    parser.add_argument(
        "--rel_height",
        type=float,
        default=0.5,
        help="Relative height for peak width calculation. Default is 0.5."
    )
    parser.add_argument(
        "--ignore_columns",
        action='store_true',
        help="Flag to determine whether to ignore specific columns ('Biology', 'Genotype', 'Cage', 'Mouse')."
    )
    parser.add_argument(
        "--max_peaks",
        action='store_true',
        help="Flag to determine if maximum peaks should be calculated."
    )

    return parser.parse_args()

def main():
    args = parse_arguments()

    # Start timing
    start_time = time.time()
    print("Loading Parquet file...")

    # Load the input DataFrame
    try:
        df_grouped = pd.read_parquet(args.input_file)
        print("Parquet file loaded successfully.")
    except Exception as e:
        print(f"Error loading Parquet file: {e}", file=sys.stderr)
        sys.exit(1)

    # Initialize LipidAnalysis class and set parameters
    analysis = LipidAnalysis(df_grouped)
    analysis.set_parameters(height=args.height, width=args.width, rel_height=args.rel_height)

    # Apply peak extraction
    print("Applying extraction...")
    peaks_df = analysis.find_lipid_peaks(max_peaks=args.max_peaks, ignore_columns=args.ignore_columns)
    print("Extraction applied successfully.")

    # Perform isomer selection
    print("Applying isomer selection...")
    isomer_df = analysis.isomer_selection(peaks_df)
    print("Isomer selection applied successfully.")

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Construct output file name
    sample_id = df_grouped['Sample'].iloc[0] if 'Sample' in df_grouped.columns else 'unknown_sample'
    output_file = os.path.join(args.output_dir, f"df_analysis_5_{sample_id}_OFF.parquet")

    # Save the results
    try:
        isomer_df.to_parquet(output_file, index=False)
        print(f"Results saved to {output_file}")
    except Exception as e:
        print(f"Error saving output file: {e}", file=sys.stderr)
        sys.exit(1)

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
