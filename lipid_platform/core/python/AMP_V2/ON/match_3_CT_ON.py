import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import logging
import argparse

# Setup logging (default level set to DEBUG; can be overridden by --log_level)
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

class MatchLipids:
    def __init__(self, mrm_database, tolerance=0.3, max_peaks=False, height=500, width=2, rel_height=0.5):
        self.mrm_database = mrm_database.copy()
        self.tolerance = tolerance
        self.max_peaks = max_peaks
        self.height = height
        self.width = width
        self.rel_height = rel_height
        
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
    def run_match_lipids(OzOFF_dir, OzON_database, df_sample_2, tolerance=0.3, retention_time_window=0.5, output=None, max_peaks=False, height=500, width=2, rel_height=0.5):
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

            matcher = MatchLipids(temp_OzON_database, tolerance=tolerance, max_peaks=max_peaks, height=height, width=width, rel_height=rel_height)

            for index, row in tqdm(OzOFF_database.iterrows(), total=len(OzOFF_database), desc="Processing OzOFF Database"):
                species = row['Species']
                retention_time_value = row['Retention_Time']
                sample_value = row['Sample']
                std_rt_dif = row['STD_RT_Dif']
                std_rt_off = row['STD_RT_OFF']  # This is the OFF retention time

                logging.info(f"Processing row {index}: Species: {species}, Sample: {sample_value}")
                logging.info(f"OzOFF RT (STD_RT_OFF) being referenced: {std_rt_off}")
                logging.info(f"Retention Time Value: {retention_time_value}, STD_RT_Dif: {std_rt_dif}")

                # Define the retention time window around the adjusted retention time value
                adjusted_retention_time_value = retention_time_value + std_rt_dif
                retention_time_window_start = adjusted_retention_time_value - retention_time_window
                retention_time_window_end = adjusted_retention_time_value + retention_time_window

                # Filter df_sample_2 based on Sample value
                filtered_df_sample_2 = df_sample_2.loc[
                    (df_sample_2['Sample'] == sample_value)
                ].copy()

                # Calculate Adjusted_RT for each row in filtered_df_sample_2
                filtered_df_sample_2['Adjusted_RT'] = filtered_df_sample_2['Retention_Time'] + std_rt_dif

                # Log Adjusted_RT for debugging
                logging.debug(f"Adjusted_RT values for filtered_df_sample_2: {filtered_df_sample_2['Adjusted_RT'].tolist()}")

                # Filter based on the retention time window using Adjusted_RT
                filtered_df_sample_2 = filtered_df_sample_2.loc[
                    (filtered_df_sample_2['Adjusted_RT'] >= retention_time_window_start) &
                    (filtered_df_sample_2['Adjusted_RT'] <= retention_time_window_end)
                ]

                if filtered_df_sample_2.empty:
                    logging.info(f"No matching entries found in df_sample_2 for Species: {species}, Sample: {sample_value} within the retention time window.")
                else:
                    logging.info(f"Found {len(filtered_df_sample_2)} matching entries in df_sample_2 for Species: {species}, Sample: {sample_value} within the retention time window.")

                    # Assign species to the filtered data
                    filtered_df_sample_2.loc[:, 'Species'] = species

                    temp_OzON_data = filtered_df_sample_2.copy()

                    logging.info(f"Assigned Species '{species}' and calculated Adjusted_RT for each row. Preparing for matching...")

                    # Ensure that the matched DataFrame includes the necessary columns
                    matched_temp_OzON_data = matcher.match_lipids_parser(temp_OzON_data)
                    matched_temp_OzON_data['STD_RT_OFF'] = row['STD_RT_OFF']
                    matched_temp_OzON_data['STD_RT_ON'] = filtered_df_sample_2['STD_RT_ON'].iloc[0]
                    matched_temp_OzON_data['STD_RT_Dif'] = std_rt_dif

                    # Log the matched retention times and difference
                    logging.info(f"Matched Lipid: {matched_temp_OzON_data['Lipid'].iloc[0]}")
                    logging.info(f"Matched STD_RT_OFF: {matched_temp_OzON_data['STD_RT_OFF'].iloc[0]}, STD_RT_ON: {matched_temp_OzON_data['STD_RT_ON']}, Calculated STD_RT_Dif: {std_rt_dif}")
                    logging.info(f"Adjusted Retention Time (Adjusted_RT) calculated individually for each matching row.")

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
    parser = argparse.ArgumentParser(description="Match Lipids between OzOFF and OzON databases.")
    
    parser.add_argument("--ozoff_dir", required=True, help="Directory containing OzOFF parquet files.")
    parser.add_argument("--ozon_database", required=True, help="Path to OzON database parquet file.")
    parser.add_argument("--sample_file", required=True, help="Path to sample parquet file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save output parquet files.")
    parser.add_argument("--tolerance", type=float, default=0.3, help="Tolerance value for matching (default: 0.3).")
    parser.add_argument("--retention_time_window", type=float, default=0.5, help="Retention time window (default: 0.5).")
    
    # Add the missing arguments
    parser.add_argument("--log_level", type=str, default="DEBUG", help="Logging level (e.g., DEBUG, INFO, WARNING).")
    parser.add_argument("--max_peaks", action="store_true", help="Flag to enable max peaks processing.")
    parser.add_argument("--height", type=float, default=500.0, help="Height parameter for processing.")
    parser.add_argument("--width", type=float, default=2.0, help="Width parameter for processing.")
    parser.add_argument("--rel_height", type=float, default=0.5, help="Relative height parameter for processing.")
    
    args = parser.parse_args()
    
    # Configure logging based on the log_level argument
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {args.log_level}')
    logging.getLogger().setLevel(numeric_level)

    # Log all configuration parameters
    logging.info("Configuration Parameters:")
    logging.info(f"OzOFF Directory: {args.ozoff_dir}")
    logging.info(f"OzON Database Path: {args.ozon_database}")
    logging.info(f"Sample File Path: {args.sample_file}")
    logging.info(f"Output Directory: {args.output_dir}")
    logging.info(f"Tolerance: {args.tolerance}")
    logging.info(f"Retention Time Window: {args.retention_time_window}")
    logging.info(f"Log Level: {args.log_level}")
    logging.info(f"Max Peaks: {args.max_peaks}")
    logging.info(f"Height: {args.height}")
    logging.info(f"Width: {args.width}")
    logging.info(f"Relative Height: {args.rel_height}")
    
    # Load data
    logging.info(f"Loading sample data from {args.sample_file}")
    df_sample_2 = pd.read_parquet(args.sample_file)

    logging.info(f"Loading OzON database from {args.ozon_database}")
    OzON_database = pd.read_parquet(args.ozon_database)
    
    # Run matching process
    logging.info(f"Starting matching process for sample: {args.sample_file}")
    OzON_results = MatchLipids.run_match_lipids(
        args.ozoff_dir, 
        OzON_database, 
        df_sample_2, 
        tolerance=args.tolerance, 
        retention_time_window=args.retention_time_window, 
        output=args.output_dir,
        max_peaks=args.max_peaks,
        height=args.height,
        width=args.width,
        rel_height=args.rel_height
    )

    # Log results
    logging.info(f"Filtered OzON Database Matches for {args.sample_file}:")
    logging.info(OzON_results.head().to_string())
