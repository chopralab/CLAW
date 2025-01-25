import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
import re

class Plot:
    def __init__(self, raw_data_file, analyzed_data_file):
        self.raw_data = pd.read_parquet(raw_data_file)
        self.analyzed_data = pd.read_parquet(analyzed_data_file)

    @staticmethod
    def extract_n_value(lipid):
        match = re.search(r'n-(\d+)', lipid)
        db_match = re.search(r'FA\((\d+:\d+)\)', lipid)
        n_value = int(match.group(1)) if match else float('inf')
        db_value = db_match.group(1) if db_match else 'Unknown'
        return n_value, db_value

    def sort_lipids(self):
        self.analyzed_data[['n_value', 'DB']] = self.analyzed_data['Lipid'].apply(lambda x: pd.Series(self.extract_n_value(x)))
        self.analyzed_data.sort_values(by=['Species', 'DB', 'n_value'], inplace=True)
        self.analyzed_data.drop(columns=['n_value'], inplace=True)

    def plot_peaks(self, project_results=None, file_name_to_save=None):
        self.sort_lipids()
        unique_lipids = self.analyzed_data['Lipid'].unique()

        for lipid in unique_lipids:
            lipid_data = self.analyzed_data[self.analyzed_data['Lipid'] == lipid]
            plt.figure(figsize=(10, 6))

            plt.plot(lipid_data['Retention_Time'], lipid_data['OzESI_Intensity'], label='Intensity')
            plt.scatter(lipid_data['Retention_Time'], lipid_data['OzESI_Intensity'], color='red')

            for _, row in lipid_data.iterrows():
                plt.annotate(f"Peak Height: {row['Peak_Height']:.2f}\nFWHM: {row['FWHM']:.2f}\nArea: {row['Peak_Area']:.2f}",
                             (row['Retention_Time'], row['OzESI_Intensity']),
                             textcoords="offset points",
                             xytext=(10, -10),
                             ha='left',
                             fontsize=8)

            plt.title(f"Peaks for Lipid: {lipid}")
            plt.xlabel('Retention Time')
            plt.ylabel('OzESI Intensity')
            plt.legend()
            plt.grid(True)

            if project_results and file_name_to_save:
                plot_folder = f'{project_results}plots/'
                os.makedirs(plot_folder, exist_ok=True)
                sample_name = lipid_data['Sample'].iloc[0]
                plt.savefig(f"{plot_folder}{file_name_to_save}_{sample_name}_{lipid}_peaks.png")
            plt.show()
            plt.close()

    def plot_data_and_peaks(self, group_type, group_value, height=1000, width=None, rel_height=0.5):
        if group_type not in ['group_by_ion', 'group_by_lipid']:
            raise ValueError(f"group_type must be 'group_by_ion' or 'group_by_lipid'")

        if group_type not in self.raw_data.columns:
            raise ValueError(f"group_type '{group_type}' is not a valid column in the data")

        group_data = self.raw_data[self.raw_data[group_type] == group_value]
        peaks, properties = find_peaks(group_data['OzESI_Intensity'], height=height, width=width)
        results_half = peak_widths(group_data['OzESI_Intensity'], peaks, rel_height=rel_height)

        plt.figure(figsize=(10, 6))
        plt.plot(group_data['Retention_Time'], group_data['OzESI_Intensity'], label='Intensity')
        plt.scatter(group_data['Retention_Time'].iloc[peaks], group_data['OzESI_Intensity'].iloc[peaks], color='red')

        for i, peak in enumerate(peaks):
            left_ip = results_half[2][i]
            right_ip = results_half[3][i]
            left_time = group_data['Retention_Time'].iloc[int(left_ip)]
            right_time = group_data['Retention_Time'].iloc[int(right_ip)]
            width_in_time = right_time - left_time

            fwhm = results_half[0][i] * (group_data['Retention_Time'].values[1] - group_data['Retention_Time'].values[0])

            plt.annotate(
                f"Peak Height: {properties['peak_heights'][i]:.2f}\nFWHM: {fwhm:.2f}\nArea: {properties['peak_heights'][i] * width_in_time:.2f}",
                (group_data['Retention_Time'].iloc[peak], group_data['OzESI_Intensity'].iloc[peak]),
                textcoords="offset points",
                xytext=(10, -10),
                ha='left',
                fontsize=8
            )

        if group_type == 'group_by_lipid':
            lipid_value = group_data['Lipid'].iloc[0] if 'Lipid' in group_data.columns else 'Unknown'
            sample_name = group_data['Sample'].iloc[0] if 'Sample' in group_data.columns else 'Unknown'
            plot_title = f"{lipid_value} {sample_name}"
        else:  # group_type == 'group_by_ion'
            parent_ion_value = group_data['Parent_Ion'].iloc[0] if 'Parent_Ion' in group_data.columns else 'Unknown'
            plot_title = f"{parent_ion_value} {sample_name}"

        plt.title(f"{plot_title}")
        plt.xlabel('Retention Time')
        plt.ylabel('OzESI Intensity')
        plt.legend()
        plt.grid(True)
        plt.show()

# Usage example
# plot = Plot('raw_data.parquet', 'analyzed_data.parquet')
# plot.plot_peaks(project_results='path_to_project_results', file_name_to_save='results_file_name')
# plot.plot_data_and_peaks(group_type='group_by_lipid', group_value='some_lipid_value')
