import os
import sys
import time
import argparse
import pandas as pd
from tqdm import tqdm
from scipy.signal import find_peaks, peak_widths
from scipy.stats import linregress
import numpy as np
import re
import matplotlib.pyplot as plt

class LipidAnalysis:
    def __init__(self, data, height=1000, width=None, rel_height=0.5):
        """
        Initialize the LipidAnalysis class.

        :param data: DataFrame containing lipid analysis data.
        :param height: Minimum height of peaks to find.
        :param width: Required width of peaks.
        :param rel_height: Relative height at which to calculate peak widths.
        """
        self.data = data
        self.height = height
        self.width = width
        self.rel_height = rel_height

    @staticmethod
    def extract_species_info(species):
        """
        Extract carbon number and double bond information from species string.

        :param species: Species string (e.g., 'd2-16:1').
        :return: Tuple of (carbon_number, double_bond).
        """
        parts = species.split(':')
        carbon_number = float(parts[0].replace('d2-', '').replace('inf', '0')) if parts[0].replace('d2-', '').replace('inf', '0').isdigit() else float('inf')
        double_bond = float(parts[1]) if len(parts) > 1 and parts[1].replace('inf', '0').isdigit() else float('inf')
        return carbon_number, double_bond

    @staticmethod
    def extract_lipid_info(lipid):
        """
        Extract carbon number and double bond information from lipid string.

        :param lipid: Lipid string (e.g., 'FA(16:1)').
        :return: Tuple of (carbon_number, double_bond).
        """
        match = re.match(r'FA\((\d+):(\d+)\)', lipid)
        if match:
            return int(match.group(1)), int(match.group(2))
        else:
            return float('inf'), float('inf')  # For unknown formats

    def find_lipid_peaks(self):
        """
        Find peaks in the lipid data.

        :return: DataFrame containing peak information.
        """
        peak_data = []
        filter_col = 'group_by_lipid'
        unique_groups = self.data[filter_col].unique()

        for group in tqdm(unique_groups, desc=f"Processing {filter_col} groups"):
            group_data = self.data[self.data[filter_col] == group].reset_index(drop=True)
            peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=self.height, width=self.width)
            results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=self.rel_height)

            retention_times = group_data['Retention_Time'].values
            ozesi_intensity = group_data['OzESI_Intensity'].values
            sampling_interval = retention_times[1] - retention_times[0] if len(retention_times) > 1 else 1

            for i, peak in enumerate(peaks):
                metadata_columns = ['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 'group_by_lipid', 'group_by_ion', 'Lipid', 
                                    'STD_RT_ON', 'STD_RT_OFF', 'STD_RT_Dif', 'STD_Peak_Intensity', 
                                    'STD_Peak_Area', 'Adjusted_RT']
                if group_data.at[peak, 'Sample'] != 'FAME':
                    metadata_columns += ['Biology', 'Genotype', 'Cage', 'Mouse']

                metadata = group_data.loc[peak, metadata_columns]

                left_ip = int(results_half[2][i])
                right_ip = int(results_half[3][i])

                left_time = group_data.at[left_ip, 'Retention_Time']
                right_time = group_data.at[right_ip, 'Retention_Time']
                width_in_time = right_time - left_time
                fwhm = results_half[0][i] * sampling_interval

                # Linear regression for left and right of the peak
                slope_left = self.calculate_slope(retention_times, ozesi_intensity, peak, direction='left')
                slope_right = self.calculate_slope(retention_times, ozesi_intensity, peak, direction='right')

                # Calculate peak area using the trapezoidal rule
                peak_area = np.trapz(ozesi_intensity[left_ip:right_ip + 1], retention_times[left_ip:right_ip + 1])

                # Calculate asymmetry factor
                asymmetry_factor = self.calculate_peak_asymmetry(peak, left_ip, right_ip, ozesi_intensity)

                # Calculate baseline drift
                baseline_drift = self.calculate_baseline_drift(retention_times, ozesi_intensity, left_ip, right_ip)

                # Detect peak shoulders
                shoulder_count = self.calculate_peak_shouldering(ozesi_intensity, peak, left_ip, right_ip)

                peak_info = {
                    'Lipid': metadata['Lipid'],
                    'Retention_Time': group_data.at[peak, 'Retention_Time'],
                    'OzESI_Intensity': group_data.at[peak, 'OzESI_Intensity'],
                    'group_by_ion': metadata['group_by_ion'],
                    'group_by_lipid': metadata['group_by_lipid'],
                    'Sample_ID': group_data.at[peak, 'Sample_ID'],
                    'Transition': group_data.at[peak, 'Transition'],
                    'Sample': metadata['Sample'],
                    'Parent_Ion': metadata['Parent_Ion'],
                    'Product_Ion': metadata['Product_Ion'],
                    'Species': metadata['Species'],
                    'Class': group_data.at[peak, 'Class'],
                    'Possible_Lipids': group_data.at[peak, 'Possible_Lipids'],
                    'Peak_Height': properties['peak_heights'][i],
                    'Prominence': properties['prominences'][i],
                    'FWHM': fwhm,
                    'Peak_Width': width_in_time,
                    'Peak_Area': peak_area,
                    'Filter_Column': filter_col,
                    'STD_RT_ON': metadata['STD_RT_ON'],
                    'STD_RT_OFF': metadata['STD_RT_OFF'],
                    'STD_RT_Dif': metadata['STD_RT_Dif'],
                    'STD_Peak_Intensity': metadata['STD_Peak_Intensity'],
                    'STD_Peak_Area': metadata['STD_Peak_Area'],
                    'Adjusted_RT': metadata['Adjusted_RT'],
                    'Slope_Left': slope_left,
                    'Slope_Right': slope_right,
                    'Asymmetry_Factor': asymmetry_factor,
                    'Baseline_Drift': baseline_drift,
                    'Shoulder_Count': shoulder_count
                }

                if group_data.at[peak, 'Sample'] != 'FAME':
                    peak_info.update({
                        'Biology': metadata['Biology'],
                        'Genotype': metadata['Genotype'],
                        'Cage': metadata['Cage'],
                        'Mouse': metadata['Mouse']
                    })

                peak_data.append(peak_info)

        peaks_df = pd.DataFrame(peak_data)

        # Normalize Peak_Area to STD_Peak_Area
        peaks_df = self.normalize_peak_area(peaks_df)

        # Sort DataFrame
        peaks_df['species_sort'] = peaks_df['Species'].apply(self.extract_species_info)
        peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

        return peaks_df

    @staticmethod
    def calculate_slope(retention_times, intensity, peak, direction='left', window=5):
        """
        Calculate the slope of the peak on either the left or right side.

        :param retention_times: Array of retention times.
        :param intensity: Array of intensity values.
        :param peak: Index of the peak apex.
        :param direction: 'left' or 'right'.
        :param window: Number of points to consider on each side.
        :return: Slope value.
        """
        if direction == 'left':
            indices = range(max(0, peak - window), peak)
        elif direction == 'right':
            indices = range(peak + 1, min(peak + 1 + window, len(intensity)))
        else:
            raise ValueError("Direction must be 'left' or 'right'.")

        if len(indices) > 1:
            slope, _, _, _, _ = linregress(retention_times[list(indices)], intensity[list(indices)])
            return slope
        else:
            return float('nan')

    @staticmethod
    def calculate_peak_asymmetry(peak, left_ip, right_ip, intensity, percentage=0.2):
        """
        Calculate the asymmetry factor of a peak at a specified height percentage.

        :param peak: Index of the peak apex.
        :param left_ip: Index of the left intersection point.
        :param right_ip: Index of the right intersection point.
        :param intensity: Array of intensity values.
        :param percentage: The percentage height at which to calculate asymmetry (default is 20%).
        :return: Asymmetry factor.
        """
        height_at_percentage = intensity[peak] * percentage
        left_percentage = left_ip
        right_percentage = right_ip

        for i in range(left_ip, peak):
            if intensity[i] >= height_at_percentage:
                left_percentage = i
                break

        for i in range(right_ip, peak, -1):
            if intensity[i] >= height_at_percentage:
                right_percentage = i
                break

        A = peak - left_percentage
        B = right_percentage - peak
        return B / A if A > 0 else float('inf')

    @staticmethod
    def calculate_baseline_drift(retention_times, intensity, left_ip, right_ip):
        """
        Calculate the baseline drift of a peak.

        :param retention_times: Array of retention times.
        :param intensity: Array of intensity values.
        :param left_ip: Index of the left intersection point.
        :param right_ip: Index of the right intersection point.
        :return: Slope of the baseline drift.
        """
        baseline_region_times = retention_times[left_ip:right_ip + 1]
        baseline_region_intensity = intensity[left_ip:right_ip + 1]
        slope, _, _, _, _ = linregress(baseline_region_times, baseline_region_intensity)
        return slope

    @staticmethod
    def calculate_peak_shouldering(intensity, peak, left_ip, right_ip):
        """
        Detect the presence of shoulders around a peak.

        :param intensity: Array of intensity values.
        :param peak: Index of the peak apex.
        :param left_ip: Index of the left intersection point.
        :param right_ip: Index of the right intersection point.
        :return: Number of shoulders detected.
        """
        region_intensity = intensity[left_ip:right_ip]
        small_peaks, _ = find_peaks(region_intensity, height=intensity[peak] * 0.5, prominence=intensity[peak] * 0.1)
        return len(small_peaks)

    @staticmethod
    def normalize_peak_area(df):
        """
        Normalize the Peak_Area to the STD_Peak_Area.

        :param df: DataFrame containing 'Peak_Area' and 'STD_Peak_Area' columns.
        :return: DataFrame with an added 'Normalized_Peak_Area' column.
        """
        df['Normalized_Peak_Area'] = df['Peak_Area'] / df['STD_Peak_Area']
        return df

    def create_max_peaks_df(self, peaks_df):
        """
        Create a DataFrame containing only the maximum peaks.

        :param peaks_df: DataFrame containing peak data.
        :return: DataFrame containing maximum peaks.
        """
        sorted_df = peaks_df.sort_values(by='OzESI_Intensity', ascending=False)
        max_peaks_df = sorted_df.groupby(['group_by_lipid', 'Sample']).first().reset_index()
        max_peaks_df['species_sort'] = max_peaks_df['Species'].apply(self.extract_species_info)
        max_peaks_df = max_peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

        columns = ['Lipid'] + [col for col in max_peaks_df.columns if col != 'Lipid']
        max_peaks_df = max_peaks_df[columns]

        return max_peaks_df

def extract_n_value(lipid):
    """
    Extract the n-value from the lipid string.

    :param lipid: String containing lipid information.
    :return: Integer n-value extracted from the lipid string.
    """
    match = re.search(r'n-(\d+)', lipid)
    n_value = int(match.group(1)) if match else float('inf')
    return n_value

def sort_dataframe(df):
    """
    Sort the DataFrame first by Species column and then by n_value column.

    :param df: DataFrame to sort.
    :return: Sorted DataFrame.
    """
    df_sorted = df.sort_values(by=['Species', 'n_value'])
    return df_sorted

def parse_arguments():
    """
    Parse command-line arguments.

    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Lipid Analysis Script")
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input Parquet file.')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the output files.')
    parser.add_argument('--height', type=int, default=1000, help='Minimum height of peaks to find.')
    parser.add_argument('--width', type=float, default=None, help='Required width of peaks.')
    parser.add_argument('--rel_height', type=float, default=0.5, help='Relative height at which to calculate peak widths.')
    parser.add_argument('--max_peaks', action='store_true', help='Flag to create a DataFrame with only maximum peaks.')
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
        print(f"Error loading Parquet file: {e}")
        sys.exit(1)

    # Initialize LipidAnalysis class and set parameters
    analysis = LipidAnalysis(df_grouped, height=args.height, width=args.width, rel_height=args.rel_height)

    # Apply peak extraction
    print("Applying extraction...")
    peaks_df = analysis.find_lipid_peaks()
    print("Extraction applied successfully.")

    # Normalize and sort DataFrame
    peaks_df['Normalized_Peak_Area'] = peaks_df['Peak_Area'] / peaks_df['STD_Peak_Area']
    peaks_df['species_sort'] = peaks_df['Species'].apply(LipidAnalysis.extract_species_info)
    peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Determine output file name
    sample_value = peaks_df['Sample'].iloc[0] if 'Sample' in peaks_df.columns else 'unknown_sample'
    output_file = os.path.join(args.output_dir, f"df_analysis_5_{sample_value}.parquet")

    # Add n_value column and sort DataFrame
    peaks_df['n_value'] = peaks_df['Lipid'].apply(extract_n_value)
    peaks_df = sort_dataframe(peaks_df)

    # Save the results
    print(f"Saving the result to {output_file}...")
    peaks_df.to_parquet(output_file, index=False)
    print("File saved successfully.")
    print("Added n_value column and sorted DataFrame by Species and n_value.")

    # Optionally create max peaks DataFrame
    if args.max_peaks:
        max_peaks_df = analysis.create_max_peaks_df(peaks_df)
        max_output_file = os.path.join(args.output_dir, f"df_analysis_5_max_{sample_value}.parquet")
        max_peaks_df.to_parquet(max_output_file, index=False)
        print(f"Max peaks DataFrame saved to {max_output_file}.")

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
