import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import pymzml
import sys  # Import the sys module
import pyarrow as pq

class MzMLParser:
    def __init__(self):
        self.transition_summed_df = pd.DataFrame(columns=[
            'Parent_Ion', 
            'Product_Ion', 
            'Intensity', 
            'Transition', 
            'Sample_ID'
        ])
        self.OzESI_df = pd.DataFrame(columns=[
            'Lipid',
            'Parent_Ion', 
            'Product_Ion', 
            'Retention_Time', 
            'OzESI_Intensity', 
            'Sample_ID', 
            'Transition'
        ])
        self.time_and_intensity_df = pd.DataFrame(columns=['Time', 'Intensity'])

    def mzml_parser(self, file_path, plot_chromatogram=False):
        rows = []
        ozesi_rows = []

        run = pymzml.run.Reader(file_path, skip_chromatogram=False)
        q1_mz = 0
        q3_mz = 0

        for spectrum in run:
            for element in spectrum.ID.split(' '):
                if 'Q1' in element:
                    q1 = element.split('=')
                    q1_mz = np.round(float(q1[1]), 1)

                if 'Q3' in element:
                    q3 = element.split('=')
                    q3_mz = np.round(float(q3[1]), 1)
                    intensity_store = np.array([intensity for _, intensity in spectrum.peaks()])
                    intensity_sum = np.sum(intensity_store)
                    transition = f"{q1_mz} -> {q3_mz}"
                    sample_id = os.path.basename(file_path)[:-5]

                    rows.append({
                        'Parent_Ion': q1_mz,
                        'Product_Ion': q3_mz,
                        'Intensity': intensity_sum,
                        'Transition': transition,
                        'Sample_ID': sample_id
                    })

                    for time, intensity in spectrum.peaks():
                        ozesi_rows.append({
                            'Parent_Ion': q1_mz,
                            'Product_Ion': q3_mz,
                            'Retention_Time': time,
                            'OzESI_Intensity': intensity,
                            'Sample_ID': sample_id,
                            'Transition': transition
                        })

        df = pd.DataFrame(rows)
        self.OzESI_df = self.OzESI_df.append(pd.DataFrame(ozesi_rows), ignore_index=True)
        self.transition_summed_df = self.transition_summed_df.append(df, ignore_index=True)
        print(f'Finished parsing mzML file: {file_path}\n')

    def mzml_parser_batch(self, folder_name, plot_chromatogram=False):
        print(f"Attempting to list files in folder: {folder_name}", flush=True)

        data_folder = os.listdir(folder_name)
        data_folder.sort()

        for file in tqdm(data_folder, desc="Parsing mzML files"):
            if file.endswith('.mzML'):
                file_path = os.path.join(folder_name, file)
                self.mzml_parser(file_path, plot_chromatogram=plot_chromatogram)

        print('Finished parsing all mzML files\n')

    def get_transition_summed_df(self):
        return self.transition_summed_df

    def get_OzESI_df(self):
        return self.OzESI_df

    def save_and_measure_size(self, df, file_prefix):
        file_sizes = {}
        
        file_path = f"{file_prefix}.parquet"
        df.to_parquet(file_path, index=False, compression='brotli')
        file_sizes['Parquet'] = os.path.getsize(file_path)

        return file_sizes

# Example usage:
if __name__ == "__main__":
    print(f"Current working directory: {os.getcwd()}", flush=True)

    if len(sys.argv) < 2:
        print("Usage: python parse_mzml_files.py <path_to_mzml_data>")
        sys.exit(1)

    mzml_data = sys.argv[1]
    print(f"mzML data folder path: {mzml_data}", flush=True)


    parser = MzMLParser()
    parser.mzml_parser_batch(mzml_data)

    transition_summed_df = parser.get_transition_summed_df()
    OzESI_df = parser.get_OzESI_df()

    # Save the DataFrames in multiple formats and measure file sizes
    transition_summed_sizes = parser.save_and_measure_size(transition_summed_df, "Projects/STD/mzml_parsed/ON/sum/df_transition_summed_1_ON_test")
    OzESI_sizes = parser.save_and_measure_size(OzESI_df, "Projects/STD/mzml_parsed/ON/df_mzml_parser_1_ON_test")

    print("File sizes for df_transition_summed_1:")
    for format, size in transition_summed_sizes.items():
        print(f"{format}: {size / (1024 * 1024):.2f} MB")
    
    print("\nFile sizes for df_mzml_parser_1:")
    for format, size in OzESI_sizes.items():
        print(f"{format}: {size / (1024 * 1024):.2f} MB")
