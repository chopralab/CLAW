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

    def find_lipid_peaks(self, output_file, user_input="OFF", max_peaks=False, ignore_columns=False):
        """
        Find peaks in lipid data.

        :param output_file: Path to save the output DataFrame.
        :param user_input: String to determine if user input is ON or OFF.
        :param max_peaks: Boolean to determine if maximum peaks should be calculated.
        :param ignore_columns: Boolean to determine whether to drop 'Biology', 'Genotype', 'Cage', 'Mouse' columns.
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
                    # metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 'group_by_lipid', 'group_by_ion', 'Lipid']]
                    metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 
                                  'group_by_lipid', 'group_by_ion', 'Lipid', 
                                  'STD', 'STD_RT_OFF']]  # Add 'STD_RT_OFF' explicitly
                    # Ignore specific columns if flagged
                    if not ignore_columns:
                        metadata['Biology'] = group_data.iloc[peak]['Biology']
                        metadata['Genotype'] = group_data.iloc[peak]['Genotype']
                        metadata['Cage'] = group_data.iloc[peak]['Cage']
                        metadata['Mouse'] = group_data.iloc[peak]['Mouse']

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
            max_peaks_df = self.create_max_peaks_df(peaks_df)
            max_csv_filename = self.add_suffix_to_filename(output_file, '_max')
            max_peaks_df.to_parquet(max_csv_filename, index=False)
            return max_peaks_df
        else:
            return peaks_df



if __name__ == "__main__":
    # Start timing
    start_time = time.time()
    print("Loading Parquet file...")

    # Parse command-line arguments
    input_file = sys.argv[1]
    height = int(sys.argv[2])
    width = int(sys.argv[3])
    rel_height = float(sys.argv[4])
    ignore_columns_flag = sys.argv[5]

    # Determine whether to ignore the specific columns
    ignore_columns = ignore_columns_flag == "ignore"

    # Load the input DataFrame
    df_grouped = pd.read_parquet(input_file)
    print("Parquet file loaded successfully.")

    # Initialize LipidAnalysis class and set parameters
    analysis = LipidAnalysis(df_grouped)
    analysis.set_parameters(height=height, width=width, rel_height=rel_height)

    # Apply peak extraction
    print("Applying extraction...")
    peaks_df = analysis.find_lipid_peaks(output_file=None, ignore_columns=ignore_columns)
    print("Extraction applied successfully.")

    # Save the results
    output_file = f"Projects/AMP/analysis/OFF/df_analysis_5_{df_grouped['Sample'].iloc[0]}_OFF.parquet"
    peaks_df.to_parquet(output_file, index=False)
    print(f"Results saved to {output_file}")

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script execution time: {elapsed_time:.2f} seconds")



# ########### update code
# import os
# import sys
# import time
# import pandas as pd
# from tqdm import tqdm
# from scipy.signal import find_peaks, peak_widths
# import re
# import matplotlib.pyplot as plt

# class LipidAnalysis:
#     def __init__(self, data):
#         """
#         Initialize the LipidAnalysis class.

#         :param data: DataFrame containing lipid analysis data.
#         """
#         self.data = data
#         self.height = 1000
#         self.width = None
#         self.rel_height = 0.5

#     def set_parameters(self, height=1000, width=None, rel_height=0.5):
#         """
#         Set parameters for peak finding.

#         :param height: Minimum height of peaks.
#         :param width: Minimum width of peaks.
#         :param rel_height: Relative height for peak width calculation.
#         """
#         self.height = height
#         self.width = width
#         self.rel_height = rel_height

#     def extract_species_info(self, species):
#         """
#         Extract carbon number and double bond information from species string.

#         :param species: String containing species information.
#         :return: Tuple containing carbon number and double bond count for sorting purposes.
#         """
#         parts = species.split(':')
#         if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
#             carbon_number = int(parts[0])
#             double_bond = int(parts[1])
#         else:
#             carbon_number = float('inf')  # Use a large number for unknown formats
#             double_bond = float('inf')
#         return carbon_number, double_bond

#     def isomer_selection(self, peaks_df):
#         """
#         Identify cis/trans isomers for lipid species with double bonds (e.g., 18:1 or 16:1).
#         The function finds the two largest peaks by intensity for species with a double bond of 1
#         and labels the first peak (based on retention time) as 'cis' and the second as 'trans'.

#         :param peaks_df: DataFrame with peak information, including Peak_Height, FWHM, etc.
#         :return: Updated DataFrame with a new column 'Isomer' containing 'cis' or 'trans' values.
#         """
#         # Filter for lipids with species ending in ':1' (e.g., 18:1, 16:1)
#         isomer_lipids = peaks_df[peaks_df['Species'].str.endswith(':1')]

#         # Initialize a column for Isomer labeling
#         peaks_df['Isomer'] = None

#         for lipid in isomer_lipids['Species'].unique():
#             lipid_data = peaks_df[peaks_df['Species'] == lipid]

#             # Find the two largest peaks based on intensity
#             if len(lipid_data) >= 2:
#                 # Sort by OzESI_Intensity to get the two largest peaks
#                 top_peaks = lipid_data.nlargest(2, 'OzESI_Intensity').sort_values(by='Retention_Time')

#                 # Label the first peak (based on retention time) as cis and the second as trans
#                 peaks_df.loc[top_peaks.index[0], 'Isomer'] = 'cis'
#                 peaks_df.loc[top_peaks.index[1], 'Isomer'] = 'trans'

#         # Return the DataFrame with isomer information and peak data
#         return peaks_df

#     def find_lipid_peaks(self, output_file, user_input="OFF", max_peaks=False, ignore_columns=False):
#         """
#         Find peaks in lipid data.

#         :param output_file: Path to save the output DataFrame.
#         :param user_input: String to determine if user input is ON or OFF.
#         :param max_peaks: Boolean to determine if maximum peaks should be calculated.
#         :param ignore_columns: Boolean to determine whether to drop 'Biology', 'Genotype', 'Cage', 'Mouse' columns.
#         :return: DataFrame containing peak data.
#         """
#         if user_input not in ["ON", "OFF"]:
#             raise ValueError("user_input must be 'ON' or 'OFF'")

#         peak_data = []

#         for filter_col in ['group_by_lipid']:
#             unique_groups = self.data[filter_col].unique()

#             for group in tqdm(unique_groups, desc=f"Processing {filter_col} groups"):
#                 group_data = self.data[self.data[filter_col] == group]
#                 peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=self.height, width=self.width)
#                 results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=self.rel_height)

#                 retention_times = group_data['Retention_Time'].values
#                 if len(retention_times) > 1:
#                     sampling_interval = retention_times[1] - retention_times[0]
#                 else:
#                     sampling_interval = 1  # Fallback value in case there's only one retention time

#                 for i, peak in enumerate(peaks):
#                     metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'Sample', 'Species', 'group_by_lipid', 'group_by_ion', 'Lipid', 'STD', 'STD_RT_OFF']]
                    
#                     # Ignore specific columns if flagged
#                     if not ignore_columns:
#                         metadata['Biology'] = group_data.iloc[peak]['Biology']
#                         metadata['Genotype'] = group_data.iloc[peak]['Genotype']
#                         metadata['Cage'] = group_data.iloc[peak]['Cage']
#                         metadata['Mouse'] = group_data.iloc[peak]['Mouse']

#                     left_ip = results_half[2][i]
#                     right_ip = results_half[3][i]
#                     left_time = group_data['Retention_Time'].iloc[int(left_ip)]
#                     right_time = group_data['Retention_Time'].iloc[int(right_ip)]
#                     width_in_time = right_time - left_time

#                     fwhm = results_half[0][i] * sampling_interval

#                     peak_data.append({
#                         'Lipid': metadata['Lipid'],
#                         'Retention_Time': group_data.iloc[peak]['Retention_Time'],
#                         'OzESI_Intensity': group_data.iloc[peak]['OzESI_Intensity'],
#                         'group_by_ion': metadata['group_by_ion'],
#                         'group_by_lipid': metadata['group_by_lipid'],
#                         'Sample_ID': group_data.iloc[peak]['Sample_ID'],
#                         'Transition': group_data.iloc[peak]['Transition'],
#                         'Sample': metadata['Sample'],
#                         'Parent_Ion': metadata['Parent_Ion'],
#                         'Product_Ion': metadata['Product_Ion'],
#                         'Species': metadata['Species'],
#                         'Class': group_data.iloc[peak]['Class'],
#                         'STD': metadata['STD'],  # Keep STD column
#                         'STD_RT_OFF': metadata['STD_RT_OFF'],  # Keep STD_RT_OFF column
#                         'Peak_Height': properties['peak_heights'][i],
#                         'FWHM': fwhm,
#                         'Peak_Width': width_in_time,
#                         'Peak_Area': properties['peak_heights'][i] * width_in_time,
#                         'Filter_Column': filter_col  # Track which column was used for filtering
#                     })

#         peaks_df = pd.DataFrame(peak_data)

#         # Sort the DataFrame by species information and return
#         peaks_df['species_sort'] = peaks_df['Species'].apply(self.extract_species_info)
#         peaks_df = peaks_df.sort_values(by=['species_sort', 'Parent_Ion', 'Sample'], ascending=[True, False, True]).drop(columns='species_sort')

#         if max_peaks:
#             max_peaks_df = self.create_max_peaks_df(peaks_df)
#             max_csv_filename = self.add_suffix_to_filename(output_file, '_max')
#             max_peaks_df.to_parquet(max_csv_filename, index=False)
#             return max_peaks_df
#         else:
#             return peaks_df




# if __name__ == "__main__":
#     # Start timing
#     start_time = time.time()
#     print("Loading Parquet file...")

#     # Parse command-line arguments
#     input_file = sys.argv[1]
#     height = int(sys.argv[2])
#     width = int(sys.argv[3])
#     rel_height = float(sys.argv[4])
#     ignore_columns_flag = sys.argv[5]

#     # Determine whether to ignore the specific columns
#     ignore_columns = ignore_columns_flag == "ignore"

#     # Load the input DataFrame
#     df_grouped = pd.read_parquet(input_file)
#     print("Parquet file loaded successfully.")

#     # Initialize LipidAnalysis class and set parameters
#     analysis = LipidAnalysis(df_grouped)
#     analysis.set_parameters(height=height, width=width, rel_height=rel_height)

#     # Apply peak extraction
#     print("Applying extraction...")
#     peaks_df = analysis.find_lipid_peaks(output_file=None, ignore_columns=ignore_columns)
#     print("Extraction applied successfully.")

#     # Perform isomer selection
#     print("Applying isomer selection...")
#     isomer_df = analysis.isomer_selection(peaks_df)  # peaks_df is now passed as an argument
#     print(isomer_df.columns)  # Debug to check for STD_RT_OFF
#     print("Isomer selection applied successfully.")

#     # Save the results
#     # Save the results
#     print(peaks_df.head())  # Print first few rows to confirm Peak_Height, FWHM, Peak_Width, Peak_Area
#     output_file = f"Projects/AMP/analysis/OFF/df_analysis_5_{df_grouped['Sample'].iloc[0]}_OFF.parquet"
#     isomer_df.to_parquet(output_file, index=False)
#     print(f"Results saved to {output_file}")

#     # End timing
#     end_time = time.time()
#     elapsed_time = end_time - start_time
#     print(f"Script execution time: {elapsed_time:.2f} seconds")
