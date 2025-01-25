#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import pymzml

class MzMLParser:
    def __init__(self, output_dir: Path):
        """Initialize parser with output directory configuration."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize DataFrames with consistent column names
        self._initialize_dataframes()

    def _initialize_dataframes(self):
        """Initialize empty DataFrames with predefined schemas."""
        self.transition_df = pd.DataFrame(columns=[
            'Parent_Ion', 'Product_Ion', 'Intensity', 
            'Transition', 'Sample_ID'
        ])
        
        self.intensity_df = pd.DataFrame(columns=[
            'Parent_Ion', 'Product_Ion', 'Retention_Time',
            'OzESI_Intensity', 'Sample_ID', 'Transition'
        ])

    def parse_spectrum(self, spectrum, file_path: Path):
        """Extract and process spectrum data."""
        q1_mz = q3_mz = 0
        
        # Parse Q1 and Q3 values from spectrum ID
        for element in spectrum.ID.split():
            if 'Q1=' in element:
                q1_mz = np.round(float(element.split('=')[1]), 1)
            elif 'Q3=' in element:
                q3_mz = np.round(float(element.split('=')[1]), 1)

        if q1_mz and q3_mz:
            intensities = np.array([intensity for _, intensity in spectrum.peaks()])
            total_intensity = np.sum(intensities)
            
            if total_intensity > 0:
                sample_id = file_path.stem  # Get filename without extension
                transition = f"{q1_mz} -> {q3_mz}"
                
                return {
                    'transition_data': {
                        'Parent_Ion': q1_mz,
                        'Product_Ion': q3_mz,
                        'Intensity': total_intensity,
                        'Transition': transition,
                        'Sample_ID': sample_id
                    },
                    'intensity_data': [
                        {
                            'Parent_Ion': q1_mz,
                            'Product_Ion': q3_mz,
                            'Retention_Time': time,
                            'OzESI_Intensity': intensity,
                            'Sample_ID': sample_id,
                            'Transition': transition
                        }
                        for time, intensity in spectrum.peaks()
                    ]
                }
        return None

    def parse_file(self, file_path: Path):
        """Parse single mzML file."""
        print(f"Parsing: {file_path}")
        
        run = pymzml.run.Reader(str(file_path), skip_chromatogram=False)
        transition_rows = []
        intensity_rows = []
        
        for spectrum in run:
            result = self.parse_spectrum(spectrum, file_path)
            if result:
                transition_rows.append(result['transition_data'])
                intensity_rows.extend(result['intensity_data'])
        
        # Update DataFrames
        if transition_rows:
            self.transition_df = pd.concat([
                self.transition_df, 
                pd.DataFrame(transition_rows)
            ], ignore_index=True)
            
        if intensity_rows:
            self.intensity_df = pd.concat([
                self.intensity_df, 
                pd.DataFrame(intensity_rows)
            ], ignore_index=True)

    def parse_directory(self, input_dir: Path):
        """Parse all mzML files in directory."""
        mzml_files = list(input_dir.glob('*.mzML'))
        print(f"Found {len(mzml_files)} mzML files in {input_dir}")
        
        for file_path in tqdm(mzml_files, desc="Parsing files"):
            self.parse_file(file_path)

    def save_results(self, prefix: str):
        """Save processed data to parquet files."""
        output_files = {
            'transition': self.output_dir / f"{prefix}_transition_summed.parquet",
            'intensity': self.output_dir / f"{prefix}_intensity.parquet"
        }
        
        for name, df in [('transition', self.transition_df), 
                        ('intensity', self.intensity_df)]:
            file_path = output_files[name]
            df.to_parquet(file_path, index=False, compression='brotli')
            size_mb = file_path.stat().st_size / (1024 * 1024)
            print(f"Saved {name} data: {file_path} ({size_mb:.2f} MB)")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Parse mzML files for mass spectrometry data')
    parser.add_argument('input_dir', type=Path, help='Directory containing mzML files')
    parser.add_argument('output_dir', type=Path, help='Directory for output files')
    parser.add_argument('--prefix', type=str, default='mzml_parsed',
                      help='Prefix for output files')
    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_arguments()
    
    parser = MzMLParser(args.output_dir)
    parser.parse_directory(args.input_dir)
    parser.save_results(args.prefix)

if __name__ == "__main__":
    main()