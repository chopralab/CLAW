import os
import sys
import time
import pandas as pd
from tqdm import tqdm
from scipy.signal import find_peaks, peak_widths
import re
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

    @staticmethod
    def extract_species_info(species):
        """
        Extract carbon number and double bond information from species string.

        :param species: String containing species information.
        :return: Tuple containing carbon number and double bond count.
        """
        parts = species.split(':')
        carbon_number = float(parts[0].replace('d2-', '').replace('inf', '0')) if parts[0].isdigit() else float('inf')
        double_bond = float(parts[1]) if len(parts) > 1 and parts[1].isdigit() else float('inf')
        return carbon_number, double_bond

    @staticmethod
    def extract_lipid_info(lipid):
        """
        Extract carbon number and double bond information from lipid string.

        :param lipid: String containing lipid information.
        :return: Tuple containing carbon number and double bond count.
        """
        match = re.match(r'FA\((\d+):(\d+)\)', lipid)
        if match:
            return int(match.group(1)), int(match.group(2))
        else:
            return float('inf'), float('inf')  # Return a large number to push unknown formats to the end

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

    def find_lipid_peaks(self, output_file, user_input="OFF", max_peaks=False):
        """
        Find peaks in lipid data.

        :param output_file: Path to save the output DataFrame.
        :param user_input: String to determine if user input is ON or OFF.
        :param max_peaks: Boolean to determine if maximum peaks should be calculated.
        :return: DataFrame containing peak data.
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

                retention_times = group_data['Retention_Time'].values
                if len(retention_times) > 1:
                    sampling_interval = retention_times[1] - retention_times[0]
                else:
                    sampling_interval = 1  # Fallback value in case there's only one retention time

                for i, peak in enumerate(peaks):
                    metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 'group_by_lipid', 'group_by_ion', 'Biology', 'Genotype', 'Cage', 'Mouse', 'Lipid']]
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
                        'Possible_Lipids': group_data.iloc[peak]['Possible_Lipids'],
                        'Biology': metadata['Biology'],
                        'Genotype': metadata['Genotype'],
                        'Cage': metadata['Cage'],
                        'Mouse': metadata['Mouse'],
                        'Peak_Height': properties['peak_heights'][i],
                        'FWHM': fwhm,
                        'Peak_Width': width_in_time,
                        'Peak_Area': properties['peak_heights'][i] * width_in_time,
                        'Filter_Column': filter_col  # Track which column was used for filtering
                    })

        peaks_df = pd.DataFrame(peak_data)

        peaks_df['species_sort'] = peaks_df['Species'].apply(self.extract_species_info)
        peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

        peaks_df.to_parquet('df_analysis_5_sorted.parquet', index=False)
        
        if max_peaks:
            max_peaks_df = self.create_max_peaks_df(peaks_df)
            max_csv_filename = self.add_suffix_to_filename(output_file, '_max')
            max_peaks_df.to_parquet(max_csv_filename, index=False)
            return max_peaks_df
        else:
            return peaks_df

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
    input_file = sys.argv[1]
    height = int(sys.argv[2])
    width = int(sys.argv[3])
    rel_height = float(sys.argv[4])
    plot_results_folder = "./results/"

    # Load the input DataFrame
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
    output_file = f"Projects/AMP/analysis/test/df_analysis_5_{sample_value}.parquet"
    
    # Save the results
    print(f"Saving the result to {output_file}...")
    peaks_df['n_value'] = peaks_df['Lipid'].apply(extract_n_value)
    peaks_df = sort_dataframe(peaks_df)
    peaks_df.to_parquet(output_file, index=False)
    print("File saved successfully.")
    print("Added n_value column and sorted DataFrame by Species and n_value.")

    # Plot the peaks
    analysis.plot_peaks(peaks_df, project_results=plot_results_folder, file_name_to_save="lipid_analysis")

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")
