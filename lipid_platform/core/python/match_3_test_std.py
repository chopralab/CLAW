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
    def run_match_lipids(OzOFF_dir, OzON_database, df_sample_2, tolerance=0.3, retention_time_window=0.5, output=None):
        logging.info("Starting lipid matching process.")
        
        sample_value = df_sample_2['Sample'].iloc[0]  # Assuming Sample column exists
        logging.info(f"Sample value: {sample_value}")

        # Find the matching OzOFF file based on Sample value
        matching_OzOFF_file = MatchLipids.find_matching_OzOFF_file(OzOFF_dir, sample_value)
        
        if matching_OzOFF_file:
            logging.info(f"Found matching OzOFF file: {matching_OzOFF_file}")
            OzOFF_database = pd.read_parquet(matching_OzOFF_file)
            
            # Directly calculate the STD_RT_Dif by subtracting the STD_RT_ON from STD_RT_OFF
            OzOFF_database['STD_RT_Dif'] = OzOFF_database['STD_RT_OFF'] - df_sample_2['STD_RT_ON']
            
            # Log the calculated STD_RT_OFF, STD_RT_ON, and STD_RT_Dif for debugging
            for _, row in OzOFF_database.iterrows():
                logging.info(f"STD_RT_OFF: {row['STD_RT_OFF']}, STD_RT_ON: {df_sample_2['STD_RT_ON'].iloc[0]}, STD_RT_Dif: {row['STD_RT_Dif']}")

            # Add STD_RT_Dif to df_sample_2
            df_sample_2['STD_RT_Dif'] = OzOFF_database['STD_RT_OFF'] - df_sample_2['STD_RT_ON']
            
            # Proceed with the lipid matching process
            OzON_results_df = pd.DataFrame()

            species_list = OzOFF_database['Species'].tolist()
            logging.info(f'SPECIES LIST: {species_list}')
            temp_OzON_database = OzON_database[OzON_database['Species'].isin(species_list)].copy()
            logging.info(f'TEMP_OzON DATABASE: {temp_OzON_database}')
            
            logging.info("Filtered OzON database based on OzOFF species.")

            matcher = MatchLipids(temp_OzON_database, tolerance=tolerance)

            for index, row in tqdm(OzOFF_database.iterrows(), total=len(OzOFF_database), desc="Processing OzOFF Database"):
                species = row['Species']
                retention_time_value = row['Retention_Time']
                sample_value = row['Sample']
                std_rt_dif = row['STD_RT_Dif']
                std_rt_off = row['STD_RT_OFF']  # This is the OFF retention time

                # Calculate the adjusted retention time value
                adjusted_retention_time_value = retention_time_value + std_rt_dif

                logging.info(f"Processing row {index}: Species: {species}, Sample: {sample_value}")
                logging.info(f"OzOFF RT (STD_RT_OFF) being referenced: {std_rt_off}")
                logging.info(f"Retention Time Value: {retention_time_value}, STD_RT_Dif: {std_rt_dif}")
                logging.info(f"Adjusted Retention Time Value: {adjusted_retention_time_value}")

                # Adjust the retention time window with the STD_RT_Dif
                adjusted_retention_time_start = adjusted_retention_time_value - retention_time_window
                adjusted_retention_time_end = adjusted_retention_time_value + retention_time_window

                # Log the adjusted retention time window and the OzOFF reference retention time (STD_RT_OFF)
                logging.debug(f"Adjusted Retention Time Window for Species '{species}' and Sample '{sample_value}': Start = {adjusted_retention_time_start}, End = {adjusted_retention_time_end}")
                logging.info(f"OzOFF Reference Retention Time (STD_RT_OFF) for Species '{species}': {std_rt_off}")

                # Filtering df_sample_2 based on adjusted retention time window and sample value
                filtered_df_sample_2 = df_sample_2.loc[
                    (df_sample_2['Sample'] == sample_value) &
                    (df_sample_2['Retention_Time'] >= adjusted_retention_time_start) &
                    (df_sample_2['Retention_Time'] <= adjusted_retention_time_end)
                ].copy()

                if filtered_df_sample_2.empty:
                    logging.info(f"No matching entries found in df_sample_2 for Species: {species}, Sample: {sample_value} within the adjusted retention time window.")
                else:
                    logging.info(f"Found {len(filtered_df_sample_2)} matching entries in df_sample_2 for Species: {species}, Sample: {sample_value} within the adjusted retention time window.")

                    # Assign species and add Adjusted_RT to the filtered data
                    filtered_df_sample_2.loc[:, 'Species'] = species
                    filtered_df_sample_2['Adjusted_RT'] = adjusted_retention_time_value  # Add Adjusted_RT column with the adjusted retention time
                    temp_OzON_data = filtered_df_sample_2.copy()

                    logging.info(f"Assigned Species '{species}' and Adjusted_RT to the filtered data. Preparing for matching...")

                    # Ensure that the matched DataFrame includes the necessary columns
                    matched_temp_OzON_data = matcher.match_lipids_parser(temp_OzON_data)
                    matched_temp_OzON_data['STD_RT_OFF'] = row['STD_RT_OFF']
                    matched_temp_OzON_data['STD_RT_ON'] = filtered_df_sample_2['STD_RT_ON'].iloc[0]
                    matched_temp_OzON_data['STD_RT_Dif'] = std_rt_dif
                    matched_temp_OzON_data['Adjusted_RT'] = adjusted_retention_time_value  # Add Adjusted_RT to the matched data

                    # Log the matched retention times and difference
                    logging.info(f"Matched Lipid: {matched_temp_OzON_data['Lipid'].iloc[0]}")
                    logging.info(f"Matched STD_RT_OFF: {matched_temp_OzON_data['STD_RT_OFF'].iloc[0]}, STD_RT_ON: {matched_temp_OzON_data['STD_RT_ON']}, Calculated STD_RT_Dif: {std_rt_dif}")
                    logging.info(f"Adjusted Retention Time (Adjusted_RT) for this match: {adjusted_retention_time_value}")

                    # Log additional verification details
                    logging.info(f"Row Index: {index}, Species: {species}")
                    logging.info(f"Adjusted Retention Time Start: {adjusted_retention_time_start}")
                    logging.info(f"Adjusted Retention Time End: {adjusted_retention_time_end}")
                    logging.info(f"Filtered Rows: {len(filtered_df_sample_2)}")

                    # Add matched data to the results DataFrame
                    OzON_results_df = pd.concat([OzON_results_df, matched_temp_OzON_data], ignore_index=True)

            # After processing, log the DataFrame column names to ensure Adjusted_RT is present
            logging.info(f"Final DataFrame Columns: {OzON_results_df.columns}")


            if output:
                output_file = os.path.join(output, f"df_match_3_{sample_value}.parquet")
                output_file = os.path.normpath(output_file)  # Normalize the path to remove any redundant slashes
                OzON_results_df.to_parquet(output_file, index=False)
                logging.info(f"Output file saved successfully to {output_file}")

            logging.info("Lipid matching process completed.")
  
            # Return results df
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
