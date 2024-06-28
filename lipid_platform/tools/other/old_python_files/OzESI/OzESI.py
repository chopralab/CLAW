import pandas as pd
import re
import matplotlib.pyplot as plt


class RawDataParser:
    def __init__(self, df_MRM, df_OzESI):
        self.df_MRM = df_MRM
        self.df_OzESI = df_OzESI
    
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
            group_columns = ['Lipid', 'Sample_ID','Biology','Genotype','Mouse','Cage']
        
        if not isinstance(group_columns, list):
            raise ValueError("group_columns must be a list of column names")
        
        invalid_columns = set(group_columns) - set(df.columns)
        if invalid_columns:
            raise ValueError(f"The following columns are not present in the DataFrame: {', '.join(invalid_columns)}")
        
        df['Group_Sample'] = df.groupby(group_columns).ngroup()
        return df
    
    def parse_data(self, retention_time_range, new_columns, group_columns=None):
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
        
        return d1d
    
    def plot_full_spectrum(self, d1a):
        plt.figure(figsize=(10, 6))
        plt.scatter(d1a['Retention_Time'], d1a['OzESI_Intensity'])
        plt.xlabel('Retention Time')
        plt.ylabel('OzESI Intensity')
        plt.title('Retention Time vs OzESI Intensity')
        plt.show()


#######################
        
import pandas as pd
from scipy.signal import find_peaks, peak_widths

class PeakAnalysis:
    def __init__(self, data):
        self.data = data

    def find_lipid_peaks(self, use_match_group=True, height=None, width=None):
        """
        Find peaks in lipid data based on retention time and intensity, and extract relevant metadata.

        Parameters:
        - use_match_group: boolean flag to determine filtering by 'Match_Group' or 'Group_Sample'.
        - height: Minimum height of peaks. Used as a threshold for peak intensity.
        - width: Minimum width of peaks in number of samples. Helps in distinguishing real peaks from noise.

        Returns:
        - peaks_df: DataFrame containing the peak data with additional metadata and calculated metrics.
        """
        filter_col = 'Match_Group' if use_match_group else 'Group_Sample'
        unique_groups = self.data[filter_col].unique()
        peak_data = []

        for group in unique_groups:
            group_data = self.data[self.data[filter_col] == group]
            peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=height, width=width)
            results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=0.5)

            for i, peak in enumerate(peaks):
                metadata = group_data.iloc[peak][['Parent_Ion', 'Product_Ion', 'FAC', 'Group_Sample', 'Match_Group', 'Biology', 'Genotype', 'Cage', 'Mouse', 'Lipid']]
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
                    'FWHM': results_half[1][i],  # Full width at half maximum
                    'Peak_Width': results_half[0][i],  # Width of the peak in samples
                    'Peak_Area': properties['peak_heights'][i] * results_half[0][i]  # Approximation of area
                })

        return pd.DataFrame(peak_data)



    
    def plot_data_and_peaks(self, group_value, group_type='Match_Group'):
        """
        Plot data and peaks for specified group value and group type, with detailed titles based on metadata.

        Parameters:
        - group_value: The specific value to filter and plot data for within the chosen group type.
        - group_type: 'Match_Group' or 'Group_Sample' to specify which group type to use for filtering.
        """
        # Validate group_type input
        if group_type not in ['Match_Group', 'Group_Sample']:
            raise ValueError("group_type must be 'Match_Group' or 'Group_Sample'")
        
        # Filter the DataFrame based on the specified group type and value
        filtered_df = self.data[self.data[group_type] == group_value]
        
        plt.figure(figsize=(10, 6))
        plt.plot(filtered_df['Retention_Time'], filtered_df['OzESI_Intensity'], color='blue', label='Intensity vs. Retention Time')
        if not filtered_df.empty:
            first_peak = filtered_df.iloc[0]
            title_info = f"{first_peak['Lipid']}, MG {first_peak['Match_Group']}, GS {first_peak['Group_Sample']}, {first_peak['Parent_Ion']}, {first_peak['FAC']}, {first_peak['Biology']}, {first_peak['Genotype']}, {first_peak['Cage']}, {first_peak['Mouse']}"
            plt.title(title_info)
        plt.xlabel('Retention Time')
        plt.ylabel('OzESI Intensity')
        plt.legend()
        plt.grid(True)
        plt.show()

        # Peak finding with predefined height threshold (modify as needed)
        peaks, properties = find_peaks(filtered_df['OzESI_Intensity'], height=50000)
        results_half = peak_widths(filtered_df['OzESI_Intensity'], peaks, rel_height=0.5)

        plt.figure(figsize=(10, 6))
        plt.plot(filtered_df['Retention_Time'], filtered_df['OzESI_Intensity'], color='gray', alpha=0.5, label='Full Data')
        plt.scatter(filtered_df.iloc[peaks]['Retention_Time'], filtered_df.iloc[peaks]['OzESI_Intensity'], color='red', zorder=5)
        for peak, left_ips, right_ips in zip(peaks, results_half[2], results_half[3]):
            plt.axvline(x=filtered_df.iloc[int(left_ips)]['Retention_Time'], color='green', linestyle=':', linewidth=1, label='Start of Peak')
            plt.axvline(x=filtered_df.iloc[int(right_ips)]['Retention_Time'], color='blue', linestyle=':', linewidth=1, label='End of Peak')
        if not filtered_df.empty:
            plt.title(title_info)
        plt.xlabel('Retention Time')
        plt.ylabel('Intensity')
        plt.legend()
        plt.grid(True)
        plt.show()

# Example usage:
# Assuming 'data' is a pandas DataFrame containing the relevant data
# peak_analysis = PeakAnalysis(data)
# peak_analysis.plot_data_and_peaks('Specific_Group_Value', group_type='Match_Group')
# peak_analysis.plot_data_and_peaks('Specific_Group_Value', group_type='Group_Sample')


import os

import os

# Create the directory if it doesn't exist
def create_base_directory(base_plot_directory):
    if not os.path.exists(base_plot_directory):
        os.makedirs(base_plot_directory)
        print(f"Directory created at {base_plot_directory}")
    else:
        print(f"Directory already exists at {base_plot_directory}")

# Define a function to generate filenames based on lipid names
def generate_filename(base_plot_directory, lipid_name):
    safe_lipid_name = lipid_name.replace("/", "-").replace(" ", "_").replace(":", "-")  # Ensure filename is safe for filesystems
    return f"{base_plot_directory}{safe_lipid_name}_OzON.png"

#### for ozone compare

# Function to save a DataFrame to an Excel file
def save_for_ozone_compare(peaks_df, project_results_directory, save_df_name):
    # Ensure the directory exists
    if not os.path.exists(project_results_directory):
        os.makedirs(project_results_directory)
        print(f"Directory created at {project_results_directory}")
    else:
        print(f"Directory already exists at {project_results_directory}")
    
    # Save the DataFrame to an Excel file
    file_path = os.path.join(project_results_directory, save_df_name)
    peaks_df.to_excel(file_path, index=False)
    print(f'peaks_df saved to excel in results folder {file_path}')
