import pandas as pd
import re
import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks, peak_widths
from OzESI_utils import create_folder



class RawDataParser:
    def __init__(self, mrm_csv_path, oze_esi_csv_path):
        self.df_MRM = pd.read_csv(mrm_csv_path)
        self.df_OzESI = pd.read_csv(oze_esi_csv_path)

    def create_match_group(self, df):
        df['Match_Group'] = df.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()
        return df

    def filter_retention_time(self, df, retention_time_range):
        df_filtered = df[(df['Retention_Time'] >= retention_time_range[0]) & (df['Retention_Time'] <= retention_time_range[1])]
        return df_filtered

    def is_within_tolerance(self, ion1, ion2, tolerance=0.3):
        return abs(ion1 - ion2) <= tolerance

    def match_lipids(self, d1a, d1b):
        if 'Match_Group' not in d1a.columns:
            d1a['Match_Group'] = d1a.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()

        d1b['Match_Group'] = d1a['Match_Group']
        d1b['Lipid'] = None

        for group in d1a['Match_Group'].unique():
            group_row = d1a[d1a['Match_Group'] == group].iloc[0]

            for _, mrm_row in self.df_MRM.iterrows():
                if self.is_within_tolerance(group_row['Parent_Ion'], mrm_row['Parent_Ion']) and self.is_within_tolerance(group_row['Product_Ion'], mrm_row['Product_Ion']):
                    d1b.loc[d1b['Match_Group'] == group, 'Lipid'] = mrm_row['Lipid']
                    break

        return d1b

    def extract_details_from_sample_id(self, df, column_name='Sample_ID', new_columns=None):
        if new_columns is None or not isinstance(new_columns, dict):
            raise ValueError("new_columns must be a dictionary with column names as keys and list of strings as values")

        for col, values in new_columns.items():
            pattern = f"(?P<{col}>{'|'.join(values)})"
            df_extracted = df[column_name].str.extract(pattern)
            df = pd.concat([df, df_extracted[[col]]], axis=1)

        return df

    def extract_fac_values(self, df):
        def extract_fac(lipid):
            if pd.isna(lipid):
                return []
            return re.findall(r'\d+:\d+', lipid)

        df['FAC'] = df['Lipid'].apply(extract_fac)
        return df

    def group(self, df, group_columns=None):
        if group_columns is None:
            group_columns = ['Lipid', 'Sample_ID', 'Biology', 'Genotype', 'Mouse', 'Cage']

        if not isinstance(group_columns, list):
            raise ValueError("group_columns must be a list of column names")

        invalid_columns = set(group_columns) - set(df.columns)
        if invalid_columns:
            raise ValueError(f"The following columns are not present in the DataFrame: {', '.join(invalid_columns)}")

        df['Group_Sample'] = df.groupby(group_columns).ngroup()
        return df

    def parse_data(self, retention_time_range, new_columns, group_columns=None, project_results=None, file_name_to_save=None, mode=None):
        self.df_MRM = self.create_match_group(self.df_MRM)

        d1 = self.df_OzESI.iloc[:, 1:]
        d1a = self.filter_retention_time(d1, retention_time_range)
        d1a = self.create_match_group(d1a)

        d1b = d1a.copy()
        d1b = self.match_lipids(d1a, d1b)

        d1c = d1b.copy()
        d1c = self.extract_details_from_sample_id(d1c, new_columns=new_columns)
        d1c = self.extract_fac_values(d1c)

        d1d = d1c.copy()
        d1d = self.group(d1d, group_columns)

        csv_data_folder = f'{project_results}csv_data/'
        create_folder(csv_data_folder)

        output_csv = f"{csv_data_folder}{file_name_to_save}_RawDataParser_{mode}.csv"
        d1d.to_csv(output_csv, index=False)
        return f"RawDataParser complete, output saved to {output_csv}"

    def save_raw_data(self, project_results, file_name_to_save, mode):
        csv_data_folder = f'{project_results}csv_data/'
        create_folder(csv_data_folder)
        self.df_MRM.to_csv(f'{csv_data_folder}{file_name_to_save}_df_MRM_{mode}.csv', index=False)
        self.df_OzESI.to_csv(f'{csv_data_folder}{file_name_to_save}_df_OzESI_{mode}.csv', index=False)
        return f"df_MRM and df_OzESI saved to {csv_data_folder}"

    def plot_full_spectrum(self, csv_file_path):
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file_path)
        
        # Plot the data
        plt.figure(figsize=(10, 6))
        plt.scatter(df['Retention_Time'], df['OzESI_Intensity'])
        plt.xlabel('Retention Time')
        plt.ylabel('OzESI Intensity')
        plt.title('Retention Time vs OzESI Intensity')
        plt.show()



class PeakAnalysis:
    def __init__(self, raw_data_csv, mode):
        self.data = pd.read_csv(raw_data_csv)
        self.mode = mode

    def find_lipid_peaks(self, use_match_group=True, height=1000, width=None, rel_height=0.5, project_results=None, file_name_to_save=None, user_input="OFF"):
        if user_input not in ["ON", "OFF"]:
            raise ValueError("user_input must be 'ON' or 'OFF'")
        
        filter_col = 'Match_Group' if use_match_group else 'Group_Sample'
        unique_groups = self.data[filter_col].unique()
        peak_data = []

        for group in unique_groups:
            group_data = self.data[self.data[filter_col] == group]
            peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=height, width=width)
            results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=rel_height)

            # Calculate sampling interval
            retention_times = group_data['Retention_Time'].values
            if len(retention_times) > 1:
                sampling_interval = retention_times[1] - retention_times[0]
            else:
                sampling_interval = 1  # Fallback value in case there's only one retention time

            print(f"Sampling Interval: {sampling_interval}")

            for i, peak in enumerate(peaks):
                metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'FAC', 'Group_Sample', 'Match_Group', 'Biology', 'Genotype', 'Cage', 'Mouse', 'Lipid']]
                left_ip = results_half[2][i]
                right_ip = results_half[3][i]
                left_time = group_data['Retention_Time'].iloc[int(left_ip)]
                right_time = group_data['Retention_Time'].iloc[int(right_ip)]
                width_in_time = right_time - left_time

                fwhm = results_half[0][i] * sampling_interval

                print(f"Group_Sample: {metadata['Group_Sample']} - Peak {i}: Start Time = {left_time}, End Time = {right_time}, Width = {width_in_time}, FWHM = {fwhm}")

                peak_data.append({
                    'Lipid': metadata['Lipid'],
                    'Retention_Time': group_data.iloc[peak]['Retention_Time'],
                    'OzESI_Intensity': group_data.iloc[peak]['OzESI_Intensity'],
                    'Match_Group': metadata['Match_Group'],
                    'Group_Sample': metadata['Group_Sample'],
                    'Sample_ID': group_data.iloc[peak]['Sample_ID'],
                    'Parent_Ion': metadata['Parent_Ion'],
                    'Product_Ion': metadata['Product_Ion'],
                    'FAC': metadata['FAC'],
                    'Biology': metadata['Biology'],
                    'Genotype': metadata['Genotype'],
                    'Cage': metadata['Cage'],
                    'Mouse': metadata['Mouse'],
                    'Peak_Height': properties['peak_heights'][i],
                    'FWHM': fwhm,
                    'Peak_Width': width_in_time,
                    'Peak_Area': properties['peak_heights'][i] * width_in_time
                })

        peaks_df = pd.DataFrame(peak_data)

        csv_data_folder = f'{project_results}csv_data/'
        create_folder(csv_data_folder)

        output_csv = f"{csv_data_folder}{file_name_to_save}_PeakAnalysis_{user_input}.csv"
        peaks_df.to_csv(output_csv, index=False)
        return peaks_df

    def plot_data_and_peaks(self, raw_data_csv, peak_analysis_csv, group_value, group_type='Match_Group', height=50000, width=None, rel_height=0.5):
        # Read the raw data CSV file into a DataFrame
        raw_df = pd.read_csv(raw_data_csv)

        # Read the peak analysis CSV file into a DataFrame
        peak_df = pd.read_csv(peak_analysis_csv)

        if group_type not in ['Match_Group', 'Group_Sample']:
            raise ValueError("group_type must be 'Match_Group' or 'Group_Sample'")

        # Filter the DataFrame based on the specified group value
        filtered_raw_df = raw_df[raw_df[group_type] == group_value]
        filtered_peak_df = peak_df[peak_df[group_type] == group_value]

        if filtered_raw_df.empty:
            print(f"No data found for group {group_value} ({group_type})")
            return

        # Calculate the sampling interval
        retention_times = filtered_raw_df['Retention_Time'].values
        sampling_interval = retention_times[1] - retention_times[0]

        # Plot the full spectrum data
        plt.figure(figsize=(10, 6))
        plt.plot(filtered_raw_df['Retention_Time'], filtered_raw_df['OzESI_Intensity'], color='blue', label='Intensity vs. Retention Time')
        if not filtered_raw_df.empty:
            first_peak = filtered_raw_df.iloc[0]
            title_info = f"{first_peak['Lipid']}, MG {first_peak['Match_Group']}, GS {first_peak['Group_Sample']}, {first_peak['Parent_Ion']}, {first_peak['FAC']}, {first_peak['Biology']}, {first_peak['Genotype']}, {first_peak['Cage']}, {first_peak['Mouse']}"
            plt.title(title_info)
        plt.xlabel('Retention Time')
        plt.ylabel('OzESI Intensity')
        plt.legend()
        plt.grid(True)
        plt.show()

        # Plot the peak data
        peaks, properties = find_peaks(filtered_raw_df['OzESI_Intensity'], height=height, width=width)
        results_half = peak_widths(filtered_raw_df['OzESI_Intensity'], peaks, rel_height=rel_height)

        # Convert peak widths to retention time units
        left_ips = results_half[2]
        right_ips = results_half[3]
        peak_widths_in_time = (right_ips - left_ips) * sampling_interval

        plt.figure(figsize=(10, 6))
        plt.plot(filtered_raw_df['Retention_Time'], filtered_raw_df['OzESI_Intensity'], color='gray', alpha=0.5, label='Full Data')
        plt.scatter(filtered_raw_df.iloc[peaks]['Retention_Time'], filtered_raw_df.iloc[peaks]['OzESI_Intensity'], color='red', label='Peaks')
        
        # Plot peak widths using the converted time units and vertical lines for start and stop times
        for left_ip, right_ip, height in zip(left_ips, right_ips, results_half[1]):
            left_time = filtered_raw_df['Retention_Time'].iloc[int(left_ip)]
            right_time = filtered_raw_df['Retention_Time'].iloc[int(right_ip)]
            plt.hlines(y=height, xmin=left_time, xmax=right_time, color='blue', linestyle='--', label='Peak Width')
            plt.axvline(x=left_time, color='green', linestyle=':', label='Peak Start' if 'Peak Start' not in plt.gca().get_legend_handles_labels()[1] else "")
            plt.axvline(x=right_time, color='purple', linestyle=':', label='Peak Stop' if 'Peak Stop' not in plt.gca().get_legend_handles_labels()[1] else "")

        plt.xlabel('Retention Time')
        plt.ylabel('OzESI Intensity')
        plt.legend()
        plt.title(f"Peaks for Group {group_value} ({group_type})")
        plt.grid(True)
        plt.show()


# # Create the directory if it doesn't exist
# def create_base_directory(base_plot_directory):
#     if not os.path.exists(base_plot_directory):
#         os.makedirs(base_plot_directory)
#         print(f"Directory created at {base_plot_directory}")
#     else:
#         print(f"Directory already exists at {base_plot_directory}")

# # Define a function to generate filenames based on lipid names
# def generate_filename(base_plot_directory, lipid_name):
#     safe_lipid_name = lipid_name.replace("/", "-").replace(" ", "_").replace(":", "-")  # Ensure filename is safe for filesystems
#     return f"{base_plot_directory}{safe_lipid_name}_OzON.png"

# #### for ozone compare

# # Function to save a DataFrame to an Excel file
# def save_for_ozone_compare(peaks_df, project_results_directory, save_df_name):
#     # Ensure the directory exists
#     if not os.path.exists(project_results_directory):
#         os.makedirs(project_results_directory)
#         print(f"Directory created at {project_results_directory}")
#     else:
#         print(f"Directory already exists at {project_results_directory}")
    
#     # Save the DataFrame to an Excel file
#     file_path = os.path.join(project_results_directory, save_df_name)
#     peaks_df.csv(file_path, index=False)
#     print(f'peaks_df saved to excel in results folder {file_path}')
