import os
import sys
import time
import pandas as pd
from tqdm import tqdm
from scipy.signal import find_peaks, peak_widths
from scipy.stats import linregress
import numpy as np
import re
import matplotlib.pyplot as plt

import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, filename='lipid_analysis.log', filemode='w',
                    format='%(levelname)s:%(message)s')


class LipidAnalysis:
    def __init__(self, data):
        """
        Initialize the LipidAnalysis class.

        :param data: DataFrame containing lipid analysis data with 'flat_baseline'.
        """
        self.data = data
        self.height = 1000
        self.width = None
        self.rel_height = 0.5

    @staticmethod
    def extract_species_info(species):
        """
        Extract carbon number and double bond information from species string.
        """
        parts = species.split(':')
        carbon_number = float(parts[0].replace('d2-', '').replace('inf', '0')) if parts[0].isdigit() else float('inf')
        double_bond = float(parts[1]) if len(parts) > 1 and parts[1].isdigit() else float('inf')
        return carbon_number, double_bond

    @staticmethod
    def extract_lipid_info(lipid):
        """
        Extract carbon number and double bond information from lipid string.
        """
        match = re.match(r'FA\((\d+):(\d+)\)', lipid)
        if match:
            return int(match.group(1)), int(match.group(2))
        else:
            return float('inf'), float('inf')  # For unknown formats

    def set_parameters(self, height=1000, width=None, rel_height=0.5):
        """
        Set parameters for peak finding.
        """
        self.height = height
        self.width = width
        self.rel_height = rel_height

    def find_lipid_peaks(self, output_file, user_input="OFF", max_peaks=False):
        """
        Find peaks in lipid data.

        :param output_file: Name of the output file.
        :param user_input: Whether to use user input (ON/OFF).
        :param max_peaks: Whether to return only maximum peaks.
        """
        if user_input not in ["ON", "OFF"]:
            raise ValueError("user_input must be 'ON' or 'OFF'")

        peak_data = []

        for filter_col in ['group_by_lipid']:
            unique_groups = self.data[filter_col].unique()

            for group in tqdm(unique_groups, desc=f"Processing {filter_col} groups"):
                group_data = self.data[self.data[filter_col] == group]
                peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=self.height, width=self.width)
                results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=self.rel_height)

                # Calculate slope values for each peak
                retention_times = group_data['Retention_Time'].values
                ozesi_intensity = group_data['OzESI_Intensity'].values
                if len(retention_times) > 1:
                    sampling_interval = retention_times[1] - retention_times[0]
                else:
                    sampling_interval = 1  # Fallback

                for i, peak in enumerate(peaks):
                    metadata_columns = ['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 'group_by_lipid', 'group_by_ion', 'Lipid', 'STD_RT_ON', 'STD_RT_OFF', 'STD_RT_Dif', 'STD_Peak_Intensity', 'STD_Peak_Area', 'Adjusted_RT']
                    if group_data.iloc[peak]['Sample'] != 'FAME':
                        metadata_columns += ['Biology', 'Genotype', 'Cage', 'Mouse']

                    metadata = group_data.iloc[peak][metadata_columns]

                    left_ip = int(results_half[2][i])
                    right_ip = int(results_half[3][i])

                    left_time = group_data['Retention_Time'].iloc[left_ip]
                    right_time = group_data['Retention_Time'].iloc[right_ip]
                    width_in_time = right_time - left_time

                    fwhm = results_half[0][i] * sampling_interval

                    # Linear regression for left and right of the peak
                    left_indices = range(max(0, peak - 5), peak)  # Left side points
                    right_indices = range(peak + 1, min(peak + 6, len(ozesi_intensity)))  # Right side points

                    if len(left_indices) > 1:
                        slope_left = linregress(retention_times[left_indices], ozesi_intensity[left_indices]).slope
                    else:
                        slope_left = float('nan')

                    if len(right_indices) > 1:
                        slope_right = linregress(retention_times[right_indices], ozesi_intensity[right_indices]).slope
                    else:
                        slope_right = float('nan')

                    # Calculate peak area using the trapezoidal rule
                    peak_region_times = retention_times[left_ip:right_ip + 1]
                    peak_region_intensity = ozesi_intensity[left_ip:right_ip + 1]
                    peak_area = np.trapz(peak_region_intensity, peak_region_times)

                    # Calculate asymmetry factor
                    asymmetry_factor = self.calculate_peak_asymmetry(peak, left_ip, right_ip, ozesi_intensity)

                    # Calculate baseline drift
                    baseline_drift = self.calculate_baseline_drift(retention_times, ozesi_intensity, left_ip, right_ip)

                    # Detect peak shoulders
                    shoulder_count = self.calculate_peak_shouldering(ozesi_intensity, peak, left_ip, right_ip)

                    # Calculate SN_Sides ratio
                    sn_sides_ratio = self.SN_sides(group_data, peak, left_ip, right_ip)

                    peak_entry = {
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
                        'Possible_Lipids': group_data.iloc[peak]['Possible_Lipids'],
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
                        'Shoulder_Count': shoulder_count,
                        'SN_Sides': sn_sides_ratio,  # Include SN_Sides here
                    }

                    if group_data.iloc[peak]['Sample'] != 'FAME':
                        peak_entry.update({
                            'Biology': metadata['Biology'],
                            'Genotype': metadata['Genotype'],
                            'Cage': metadata['Cage'],
                            'Mouse': metadata['Mouse']
                        })

                    peak_data.append(peak_entry)

        peaks_df = pd.DataFrame(peak_data)

        # Merge flat_baseline into peaks_df
        peaks_df = peaks_df.merge(self.data[['group_by_lipid', 'Retention_Time', 'flat_baseline']],
                                  on=['group_by_lipid', 'Retention_Time'], how='left')

        # Apply SN_flat function to compute S/N ratio
        peaks_df['SN_Flat'] = peaks_df.apply(self.SN_flat, axis=1)

        # Reorder columns to place SN_Sides next to SN_Flat
        cols = list(peaks_df.columns)
        # Find the positions of SN_Flat and SN_Sides
        sn_flat_index = cols.index('SN_Flat')
        sn_sides_index = cols.index('SN_Sides')
        # Remove SN_Sides from its current position
        cols.pop(sn_sides_index)
        # Insert SN_Sides right after SN_Flat
        cols.insert(sn_flat_index + 1, 'SN_Sides')
        # Reorder the DataFrame columns
        peaks_df = peaks_df[cols]

        # Normalize Peak_Area to STD_Peak_Area
        peaks_df = self.normalize_peak_area(peaks_df)

        peaks_df['species_sort'] = peaks_df['Species'].apply(self.extract_species_info)
        peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

        if output_file:
            peaks_df.to_parquet(output_file, index=False)

        if max_peaks:
            max_peaks_df = self.create_max_peaks_df(peaks_df)
            max_csv_filename = self.add_suffix_to_filename(output_file, '_max')
            max_peaks_df.to_parquet(max_csv_filename, index=False)
            return max_peaks_df
        else:
            return peaks_df

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
        """
        baseline_region_times = retention_times[left_ip:right_ip + 1]
        baseline_region_intensity = intensity[left_ip:right_ip + 1]
        slope, _, _, _, _ = linregress(baseline_region_times, baseline_region_intensity)
        return slope

    @staticmethod
    def calculate_peak_shouldering(intensity, peak, left_ip, right_ip):
        """
        Detect the presence of shoulders around a peak.
        """
        small_peaks, _ = find_peaks(intensity[left_ip:right_ip], height=intensity[peak] * 0.5, prominence=intensity[peak] * 0.1)
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

    @staticmethod
    def add_suffix_to_filename(filename, suffix):
        """
        Add a suffix to the filename before the file extension.

        :param filename: Original filename.
        :param suffix: Suffix to add.
        :return: Modified filename with suffix.
        """
        base, ext = os.path.splitext(filename)
        return f"{base}{suffix}{ext}"

    def SN_flat(self, row):
        """
        Calculate the S/N ratio using the flat_baseline value from the row.

        :param row: A pandas Series representing a row in the DataFrame.
        :return: S/N ratio calculated using the flat_baseline value.
        """
        # Extract the peak intensity from the row
        peak_intensity = row['OzESI_Intensity']

        # Extract the flat_baseline value from the row
        noise_level = row['flat_baseline']

        # Avoid division by zero
        if noise_level == 0 or np.isnan(noise_level):
            logging.warning(f"Noise level is zero or NaN for row with index {row.name}. Cannot calculate SN_Flat.")
            return float('nan')
        else:
            # Calculate S/N ratio
            sn_ratio = peak_intensity / noise_level
            logging.debug(f"Calculated SN_Flat ratio for row {row.name}: {sn_ratio}")

        return sn_ratio

    def SN_sides(self, group_data, peak_idx, left_ip, right_ip, window_size=50):
        """
        Calculate the S/N ratio using the noise regions before and after the peak.

        :param group_data: DataFrame containing data for a specific group.
        :param peak_idx: Index of the peak in group_data.
        :param left_ip: Left interpolation point from peak_widths.
        :param right_ip: Right interpolation point from peak_widths.
        :param window_size: Number of points to include in noise regions.
        :return: S/N ratio calculated using noise regions adjacent to the peak.
        """
        # Extract the peak intensity
        peak_intensity = group_data['OzESI_Intensity'].iloc[peak_idx]

        # Define noise regions before and after the peak
        # Avoid including the peak itself
        noise_start_before = max(0, left_ip - window_size)
        noise_end_before = left_ip - 1  # Up to the start of the peak

        noise_start_after = right_ip + 1  # After the end of the peak
        noise_end_after = min(len(group_data) - 1, right_ip + window_size)

        noise_region_before = group_data['OzESI_Intensity'].iloc[noise_start_before:noise_end_before]
        noise_region_after = group_data['OzESI_Intensity'].iloc[noise_start_after:noise_end_after]

        # Combine noise regions
        noise_region = pd.concat([noise_region_before, noise_region_after])

        # Ensure the noise region is not empty
        if noise_region.empty:
            # Handle empty noise region
            sn_ratio = float('nan')
        else:
            # Calculate noise level (standard deviation)
            noise_level = np.std(noise_region)

            # Avoid division by zero
            if noise_level == 0 or np.isnan(noise_level):
                sn_ratio = float('nan')
            else:
                # Calculate S/N ratio
                sn_ratio = peak_intensity / noise_level

        return sn_ratio

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

if __name__ == "__main__":
    # Start timing
    start_time = time.time()
    print("Loading Parquet file...")

    # Parse command-line arguments
    if len(sys.argv) < 5:
        print("Usage: script.py input_file height width rel_height")
        sys.exit(1)

    input_file = sys.argv[1]
    height_arg = sys.argv[2]
    width_arg = sys.argv[3]
    rel_height = float(sys.argv[4])

    # Handle 'None' for width
    if width_arg.lower() == 'none':
        width = None
    else:
        width = int(width_arg)

    height = int(height_arg)

    plot_results_folder = "./results/"

    # Load the input DataFrame with flat_baseline
    df_grouped = pd.read_parquet(input_file)
    print("Parquet file loaded successfully.")

    # Initialize LipidAnalysis class and set parameters
    analysis = LipidAnalysis(df_grouped)
    analysis.set_parameters(height=height, width=width, rel_height=rel_height)

    # Apply peak extraction
    print("Applying extraction...")
    peaks_df = analysis.find_lipid_peaks(output_file=None)
    print("Extraction applied successfully.")

    # Get the sample value from the DataFrame for naming the output file
    sample_value = peaks_df['Sample'].iloc[0] if 'Sample' in peaks_df.columns else 'unknown_sample'
    output_file = f"Projects/STD/analysis/ON/df_analysis_5_{sample_value}.parquet"

    # Save the results
    print(f"Saving the result to {output_file}...")
    peaks_df['n_value'] = peaks_df['Lipid'].apply(extract_n_value)
    peaks_df = sort_dataframe(peaks_df)
    peaks_df.to_parquet(output_file, index=False)
    print("File saved successfully.")
    print("Added n_value column and sorted DataFrame by Species and n_value.")

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")
