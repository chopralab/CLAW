import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import logging

class MatchLipids:
    def __init__(self, mrm_database, tolerance=0.3):
        self.mrm_database = mrm_database
        self.tolerance = tolerance
        self._round_database_ions()

    def _round_database_ions(self):
        for col in ['Parent_Ion', 'Product_Ion']:
            self.mrm_database[col] = np.round(self.mrm_database[col], 1)

    def match_lipids_parser(self, df):
        for col in ['Parent_Ion', 'Product_Ion']:
            df[col] = np.round(df[col], 1)
        df['OzESI_Intensity'] = np.round(df['OzESI_Intensity'], 0)
        
        matched_data = []
        for _, row in tqdm(df.iterrows(), total=len(df)):
            matches = (self.within_tolerance(self.mrm_database['Parent_Ion'], row['Parent_Ion']) & 
                      self.within_tolerance(self.mrm_database['Product_Ion'], row['Product_Ion']))
            
            matched_data.append({
                'Lipid': ' | '.join(self.mrm_database.loc[matches, 'Lipid']) if matches.any() else '',
                'Class': ' | '.join(self.mrm_database.loc[matches, 'Class']) if matches.any() else ''
            })
            
        for key in ['Lipid', 'Class']:
            df[key] = [d[key] for d in matched_data]
        return df

    @staticmethod
    def within_tolerance(values1, values2, tolerance=None):
        return np.abs(values1 - values2) <= (tolerance or 0.3)

    @staticmethod
    def save_results(df, output_dir, sample_name='unknown'):
        output_file = os.path.join(output_dir, f"df_match_3_{sample_name}.parquet")
        df.to_parquet(output_file, index=False)
        return output_file

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--database', required=True, help='Path to OzOFF database')
    parser.add_argument('--input', required=True, help='Path to sample file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--tolerance', type=float, default=0.3, help='Matching tolerance')
    parser.add_argument('--log-level', default='INFO', help='Logging level')
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)
    
    database = pd.read_parquet(args.database)
    sample_data = pd.read_parquet(args.input)
    
    matcher = MatchLipids(database, args.tolerance)
    results = matcher.match_lipids_parser(sample_data)
    
    sample_name = sample_data.get('Sample', ['unknown']).iloc[0]
    output_file = matcher.save_results(results, args.output, sample_name)
    logging.info(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()