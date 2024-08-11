import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MatchLipids:
    def __init__(self, mrm_database, tolerance=0.3):
        self.mrm_database = mrm_database.copy()  # Ensure we work with a copy
        self.tolerance = tolerance
        
        # Round the ion values in the MRM database for easier comparison
        self.mrm_database['Parent_Ion'] = np.round(self.mrm_database['Parent_Ion'], 1)
        self.mrm_database['Product_Ion'] = np.round(self.mrm_database['Product_Ion'], 1)

    def within_tolerance(self, values1, values2):
        return np.abs(values1 - values2) <= self.tolerance

    def match_lipids_parser(self, df):
        logging.info("Starting lipid matching parser.")
        
        # Round the ion values in the input DataFrame for easier comparison
        df.loc[:, 'Parent_Ion'] = np.round(df['Parent_Ion'], 1)
        df.loc[:, 'Product_Ion'] = np.round(df['Product_Ion'], 1)
        df.loc[:, 'OzESI_Intensity'] = np.round(df['OzESI_Intensity'], 0)
        
        matched_lipids = []
        matched_classes = []

        for _, row in tqdm(df.iterrows(), total=len(df), desc="Matching Lipids"):
            parent_ion_matches = self.within_tolerance(self.mrm_database['Parent_Ion'].values, row['Parent_Ion'])
            product_ion_matches = self.within_tolerance(self.mrm_database['Product_Ion'].values, row['Product_Ion'])
            matches = parent_ion_matches & product_ion_matches
            
            if np.any(matches):
                matched_lipids.append(' | '.join(self.mrm_database.loc[matches, 'Lipid'].values))
                matched_classes.append(' | '.join(self.mrm_database.loc[matches, 'Class'].values))
            else:
                matched_lipids.append('')
                matched_classes.append('')
        
        df['Lipid'] = matched_lipids
        df['Class'] = matched_classes
        
        logging.info("Lipid matching parser completed.")
        return df

    @staticmethod
    def run_match_lipids(OzOFF_database, OzON_database, df_sample_2, tolerance=0.3, retention_time_window=0.5, output=None):
        logging.info("Starting lipid matching process.")
        
        OzON_results_df = pd.DataFrame()

        species_list = OzOFF_database['Species'].tolist()
        print('SPECIES LIST:', species_list)
        temp_OzON_database = OzON_database[OzON_database['Species'].isin(species_list)].copy()
        print('TEMP_OzON DATABASE:', temp_OzON_database)
        
        logging.info("Filtered OzON database based on OzOFF species.")

        matcher = MatchLipids(temp_OzON_database, tolerance=tolerance)

        for index, row in tqdm(OzOFF_database.iterrows(), total=len(OzOFF_database), desc="Processing OzOFF Database"):
            species = row['Species']
            retention_time_value = row['Retention_Time']
            sample_value = row['Sample']
            
            filtered_df_sample_2 = df_sample_2.loc[
                (df_sample_2['Sample'] == sample_value) &
                (df_sample_2['Retention_Time'] >= (retention_time_value - retention_time_window)) &
                (df_sample_2['Retention_Time'] <= (retention_time_value + retention_time_window))
            ].copy()

            if not filtered_df_sample_2.empty:
                filtered_df_sample_2.loc[:, 'Species'] = species
                temp_OzON_data = filtered_df_sample_2.copy()
                matched_temp_OzON_data = matcher.match_lipids_parser(temp_OzON_data)
                OzON_results_df = pd.concat([OzON_results_df, matched_temp_OzON_data], ignore_index=True)

        if output:
            sample_value = df_sample_2['Sample'].iloc[0] if 'Sample' in df_sample_2.columns else 'unknown_sample'
            output_file = os.path.join(output, f"df_match_3_{sample_value}.parquet")
            output_file = os.path.normpath(output_file)  # Normalize the path to remove any redundant slashes
            OzON_results_df.to_parquet(output_file, index=False)
            logging.info(f"Output file saved successfully to {output_file}")

        logging.info("Lipid matching process completed.")
        return OzON_results_df

# Example usage:
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python match_3_test.py <OzOFF_database.parquet> <OzON_database.parquet> <sample_file.parquet> <output_dir>")
        sys.exit(1)

    OzOFF_database_path = sys.argv[1]
    OzON_database_path = sys.argv[2]
    sample_path = sys.argv[3]
    output_dir = sys.argv[4]

    tolerance = 0.9
    retention_time_window = 1.5

    logging.info(f"Loading OzOFF database from {OzOFF_database_path}")
    OzOFF_database = pd.read_parquet(OzOFF_database_path)
    
    logging.info(f"Loading OzON database from {OzON_database_path}")
    OzON_database = pd.read_parquet(OzON_database_path)
    
    logging.info(f"Loading sample data from {sample_path}")
    df_sample_2 = pd.read_parquet(sample_path)

    logging.info(f"Starting matching process for sample: {sample_path}")
    OzON_results = MatchLipids.run_match_lipids(
        OzOFF_database, 
        OzON_database, 
        df_sample_2, 
        tolerance=tolerance, 
        retention_time_window=retention_time_window, 
        output=output_dir
    )

    logging.info(f"Filtered OzON Database Matches for {sample_path}:")
    logging.info(OzON_results.head().to_string())
