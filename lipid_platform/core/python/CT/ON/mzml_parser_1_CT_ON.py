import os
import sys
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
import pymzml

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
        self.created_files = []  # Track created files

    def mzml_parser(self, file_path):
        transition_rows = []
        ozesi_rows = []

        run = pymzml.run.Reader(file_path, skip_chromatogram=False)

        for spectrum in run:
            q1_mz = q3_mz = 0.0
            for element in spectrum.ID.split(' '):
                if 'Q1=' in element:
                    q1_mz = round(float(element.split('=')[1]), 1)
                elif 'Q3=' in element:
                    q3_mz = round(float(element.split('=')[1]), 1)

            if q1_mz and q3_mz:
                intensity_store = np.array([intensity for _, intensity in spectrum.peaks()])
                intensity_sum = intensity_store.sum()
                transition = f"{q1_mz} -> {q3_mz}"
                sample_id = os.path.basename(file_path).replace('.mzML', '')

                transition_rows.append({
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

        # Append data to DataFrames using pd.concat
        if transition_rows:
            df_transition = pd.DataFrame(transition_rows)
            self.transition_summed_df = pd.concat([self.transition_summed_df, df_transition], ignore_index=True)

        if ozesi_rows:
            df_ozesi = pd.DataFrame(ozesi_rows)
            self.OzESI_df = pd.concat([self.OzESI_df, df_ozesi], ignore_index=True)

        print(f'Finished parsing mzML file: {file_path}\n')

    def mzml_parser_batch(self, folder_path):
        print(f"Listing files in folder: {folder_path}", flush=True)

        try:
            data_files = sorted([
                os.path.join(folder_path, file) 
                for file in os.listdir(folder_path) 
                if file.endswith('.mzML')
            ])
        except FileNotFoundError:
            print(f"Error: The folder '{folder_path}' does not exist.", file=sys.stderr)
            sys.exit(1)

        for file_path in tqdm(data_files, desc="Parsing mzML files"):
            self.mzml_parser(file_path)

        print('Finished parsing all mzML files\n')

    def get_transition_summed_df(self):
        return self.transition_summed_df

    def get_OzESI_df(self):
        return self.OzESI_df

    def save_and_measure_size(self, df, output_prefix, compression='brotli'):
        file_sizes = {}
        
        parquet_path = f"{output_prefix}.parquet"
        df.to_parquet(parquet_path, index=False, compression=compression)
        file_sizes['Parquet'] = os.path.getsize(parquet_path)
        self.created_files.append(parquet_path)

        return file_sizes

    def print_debug_summary(self):
        unique_samples = self.OzESI_df['Sample_ID'].nunique()
        print(f"Number of unique Sample_ID values: {unique_samples}")

        print("\nFiles created:")
        for file in self.created_files:
            print(file)

        print(f"\nTotal number of files created: {len(self.created_files)}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse mzML files and generate summarized data.")
    parser.add_argument(
        '--input_folder', 
        type=str, 
        required=True, 
        help="Path to the folder containing mzML files."
    )
    parser.add_argument(
        '--output_transition', 
        type=str, 
        required=True, 
        help="Output prefix path for the transition summed DataFrame."
    )
    parser.add_argument(
        '--output_ozesi', 
        type=str, 
        required=True, 
        help="Output prefix path for the OzESI DataFrame."
    )
    return parser.parse_args()

def main():
    args = parse_arguments()

    print(f"Current working directory: {os.getcwd()}", flush=True)
    print(f"mzML data folder path: {args.input_folder}", flush=True)

    parser = MzMLParser()
    parser.mzml_parser_batch(args.input_folder)

    transition_summed_df = parser.get_transition_summed_df()
    ozesi_df = parser.get_OzESI_df()

    # Save the DataFrames and measure file sizes
    transition_sizes = parser.save_and_measure_size(
        transition_summed_df, 
        args.output_transition
    )
    ozesi_sizes = parser.save_and_measure_size(
        ozesi_df, 
        args.output_ozesi
    )

    print("\nFile sizes for transition summed DataFrame:")
    for fmt, size in transition_sizes.items():
        print(f"{fmt}: {size / (1024 * 1024):.2f} MB")

    print("\nFile sizes for OzESI DataFrame:")
    for fmt, size in ozesi_sizes.items():
        print(f"{fmt}: {size / (1024 * 1024):.2f} MB")

    # Print the debugging summary
    parser.print_debug_summary()

if __name__ == "__main__":
    main()
