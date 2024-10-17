import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import logging

# Setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

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
        # Strip any leading/trailing spaces in df column names
        df.columns = df.columns.str.strip()
        
        # Log the columns to check if 'Isomer' exists
        logging.info(f"Columns in mrm_database: {self.mrm_database.columns}")
        logging.info(f"Columns in input DataFrame: {df.columns}")
        
        # Round the ion values in the input DataFrame for easier comparison
        df.loc[:, 'Parent_Ion'] = np.round(df['Parent_Ion'], 1)
        df.loc[:, 'Product_Ion'] = np.round(df['Product_Ion'], 1)
        df.loc[:, 'OzESI_Intensity'] = np.round(df['OzESI_Intensity'], 0)
        
        matched_lipids = []
        matched_classes = []
        matched_isomers = []  # Track matched isomers
        
        # Check if the Isomer column exists in mrm_database
        isomer_exists = 'Isomer' in self.mrm_database.columns

        for _, row in tqdm(df.iterrows(), total=len(df), desc="Matching Lipids"):
            parent_ion_matches = self.within_tolerance(self.mrm_database['Parent_Ion'].values, row['Parent_Ion'])
            product_ion_matches = self.within_tolerance(self.mrm_database['Product_Ion'].values, row['Product_Ion'])
            matches = parent_ion_matches & product_ion_matches
            
            if np.any(matches):
                matched_lipids.append(' | '.join(self.mrm_database.loc[matches, 'Lipid'].values))
                matched_classes.append(' | '.join(self.mrm_database.loc[matches, 'Class'].values))
                
                # If the Isomer column exists, capture Isomer values
                if isomer_exists:
                    matched_isomers.append(' | '.join(self.mrm_database.loc[matches, 'Isomer'].astype(str).values))
                else:
                    matched_isomers.append('')  # Leave empty if Isomer column doesn't exist
            else:
                matched_lipids.append('')
                matched_classes.append('')
                matched_isomers.append('')  # No match found for Isomer
            
        df['Lipid'] = matched_lipids
        df['Class'] = matched_classes
        df['Isomer'] = matched_isomers  # Add the matched Isomer column
        
        logging.info("Lipid matching parser completed.")
        return df


    @staticmethod
    def run_match_lipids(OzOFF_dir, OzON_database, df_sample_2, tolerance=0.3, retention_time_window=0.5, output=None):
        logging.info("Starting lipid matching process.")
        
        sample_value = df_sample_2['Sample'].iloc[0]
        logging.info(f"Sample value: {sample_value}")

        matching_OzOFF_file = MatchLipids.find_matching_OzOFF_file(OzOFF_dir, sample_value)
        
        if matching_OzOFF_file:
            logging.info(f"Found matching OzOFF file: {matching_OzOFF_file}")
            OzOFF_database = pd.read_parquet(matching_OzOFF_file)

            logging.info(f"Columns in OzOFF_database: {OzOFF_database.columns}")
            logging.info(f"Data in OzOFF_database:\n{OzOFF_database.head().to_string()}")

            OzOFF_database['STD_RT_Dif'] = OzOFF_database['STD_RT_OFF'] - df_sample_2['STD_RT_ON']
            df_sample_2['STD_RT_Dif'] = OzOFF_database['STD_RT_OFF'] - df_sample_2['STD_RT_ON']

            OzON_results_df = pd.DataFrame()
            species_list = OzOFF_database['Species'].tolist()
            temp_OzON_database = OzON_database[OzON_database['Species'].isin(species_list)].copy()
            
            matcher = MatchLipids(temp_OzON_database, tolerance=tolerance)

            for index, row in tqdm(OzOFF_database.iterrows(), total=len(OzOFF_database), desc="Processing OzOFF Database"):
                species = row['Species']
                retention_time_value = row['Retention_Time']
                sample_value = row['Sample']
                std_rt_dif = row['STD_RT_Dif']
                isomer_value = row['Isomer']  # Get the Isomer value

                adjusted_retention_time_value = retention_time_value + std_rt_dif

                filtered_df_sample_2 = df_sample_2.loc[
                    (df_sample_2['Sample'] == sample_value) &
                    (df_sample_2['Retention_Time'] >= (adjusted_retention_time_value - retention_time_window)) &
                    (df_sample_2['Retention_Time'] <= (adjusted_retention_time_value + retention_time_window))
                ].copy()

                if not filtered_df_sample_2.empty:
                    # Assign species, isomer, and Adjusted_RT to the filtered data
                    filtered_df_sample_2['Species'] = species
                    filtered_df_sample_2['Isomer'] = isomer_value  # Add the Isomer value here
                    filtered_df_sample_2['Adjusted_RT'] = adjusted_retention_time_value  # Add Adjusted_RT column

                    matched_temp_OzON_data = matcher.match_lipids_parser(filtered_df_sample_2)
                    matched_temp_OzON_data['STD_RT_OFF'] = row['STD_RT_OFF']
                    matched_temp_OzON_data['STD_RT_ON'] = filtered_df_sample_2['STD_RT_ON'].iloc[0]
                    matched_temp_OzON_data['STD_RT_Dif'] = std_rt_dif
                    matched_temp_OzON_data['Isomer'] = isomer_value  # Ensure Isomer is kept in the final results

                    # Concatenate the matched data to the final DataFrame
                    OzON_results_df = pd.concat([OzON_results_df, matched_temp_OzON_data], ignore_index=True)

            if output:
                output_file = os.path.join(output, f"df_match_3_{sample_value}.parquet")
                OzON_results_df.to_parquet(output_file, index=False)
                logging.info(f"Output file saved successfully to {output_file}")

            logging.info("Lipid matching process completed.")
            return OzON_results_df
        else:
            logging.error(f"No matching OzOFF file found for Sample: {sample_value}")
            sys.exit(1)



    @staticmethod
    def find_matching_OzOFF_file(OzOFF_dir, sample_value):
        logging.info(f"Listing files in directory: {OzOFF_dir}")
        files = os.listdir(OzOFF_dir)
        sample_key = sample_value
        for file in files:
            if sample_key in file and file.endswith('.parquet'):
                return os.path.join(OzOFF_dir, file)
        return None


# Example usage:
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python match_3_test_indi.py <OzOFF_directory> <OzON_database.parquet> <sample_file.parquet> <output_dir>")
        sys.exit(1)

    OzOFF_dir = sys.argv[1]
    OzON_database_path = sys.argv[2]
    sample_path = sys.argv[3]
    output_dir = sys.argv[4]

    tolerance = 0.3
    retention_time_window = 0.5

    logging.info(f"OzOFF directory: {OzOFF_dir}")

    logging.info(f"Loading sample data from {sample_path}")
    df_sample_2 = pd.read_parquet(sample_path)

    logging.info(f"Loading OzON database from {OzON_database_path}")
    OzON_database = pd.read_parquet(OzON_database_path)
    
    logging.info(f"Starting matching process for sample: {sample_path}")
    OzON_results = MatchLipids.run_match_lipids(
        OzOFF_dir, 
        OzON_database, 
        df_sample_2, 
        tolerance=tolerance, 
        retention_time_window=retention_time_window, 
        output=output_dir
    )

    logging.info(f"Filtered OzON Database Matches for {sample_path}:")
    logging.info(OzON_results.head().to_string())
