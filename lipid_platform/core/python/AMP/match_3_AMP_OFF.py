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
        """
        Initialize the MatchLipids class.

        :param mrm_database: DataFrame containing the MRM database with ion information.
        :param tolerance: Tolerance value for matching ion values.
        """
        self.mrm_database = mrm_database
        self.tolerance = tolerance
        
        # Round the ion values in the MRM database for easier comparison
        self.mrm_database['Parent_Ion'] = np.round(self.mrm_database['Parent_Ion'], 1)
        self.mrm_database['Product_Ion'] = np.round(self.mrm_database['Product_Ion'], 1)

    def within_tolerance(self, values1, values2):
        """
        Checks if two arrays of values are within a specified tolerance.

        :param values1: First array of values to compare.
        :param values2: Second array of values to compare.
        :return: Array of booleans indicating if the values are within tolerance.
        """
        return np.abs(values1 - values2) <= self.tolerance

    def match_lipids_parser(self, df):
        """
        Performs lipid matching by vectorized comparison of ion values.

        :param df: DataFrame containing ion information to be matched.
        :return: DataFrame with matched lipid and class information if matches were found.
        """
        logging.info("Starting lipid matching parser.")
        
        # Round the ion values in the input DataFrame for easier comparison
        df['Parent_Ion'] = np.round(df['Parent_Ion'], 1)
        df['Product_Ion'] = np.round(df['Product_Ion'], 1)
        df['OzESI_Intensity'] = np.round(df['OzESI_Intensity'], 0)
        
        matched_lipids = []  # List to store matched lipid names
        matched_classes = []  # List to store matched lipid classes

        # Iterate through each row of the input DataFrame with progress bar
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Matching Lipids"):
            # Check for matches within the specified tolerance
            parent_ion_matches = self.within_tolerance(self.mrm_database['Parent_Ion'].values, row['Parent_Ion'])
            product_ion_matches = self.within_tolerance(self.mrm_database['Product_Ion'].values, row['Product_Ion'])
            matches = parent_ion_matches & product_ion_matches
            
            # If matches are found, append the matched lipid and class information
            if np.any(matches):
                matched_lipids.append(' | '.join(self.mrm_database.loc[matches, 'Lipid'].values))
                matched_classes.append(' | '.join(self.mrm_database.loc[matches, 'Class'].values))
            else:
                matched_lipids.append('')
                matched_classes.append('')
        
        # Add the matched lipid and class information to the DataFrame
        df['Lipid'] = matched_lipids
        df['Class'] = matched_classes
        
        logging.info("Lipid matching parser completed.")
        return df

    @staticmethod
    def run_match_lipids(matcher, df_sample_2, output_dir):
        """
        Run the lipid matching process and save the results.

        Parameters:
        matcher (MatchLipids): Instance of MatchLipids initialized with the MRM database.
        df_sample_2 (pd.DataFrame): DataFrame containing sample data to match.
        output_dir (str): Directory to save the matched results.

        Returns:
        pd.DataFrame: DataFrame containing the matched and filtered results.
        """
        logging.info("Starting lipid matching process.")
        
        # Match lipids for the entire sample DataFrame
        matched_results_df = matcher.match_lipids_parser(df_sample_2)
        
        # Save the matched results to a Parquet file in the output directory
        sample_value = df_sample_2['Sample'].iloc[0] if 'Sample' in df_sample_2.columns else 'unknown_sample'
        output_file = os.path.join(output_dir, f"df_match_3_{sample_value}.parquet")
        matched_results_df.to_parquet(output_file, index=False)
        logging.info(f"Matched lipid results saved to {output_file}")
        
        logging.info("Lipid matching process completed.")
        return matched_results_df

# Example usage:
if __name__ == "__main__":
    # Ensure the correct number of arguments are provided
    if len(sys.argv) < 4:
        print("Usage: python match_3_OFF.py <OzOFF_database.parquet> <sample_file.parquet> <output_dir>")
        sys.exit(1)

    # Parse command-line arguments
    OzOFF_database_path = sys.argv[1]
    sample_path = sys.argv[2]
    output_dir = sys.argv[3]

    # Load the input DataFrames
    logging.info(f"Loading OzOFF database from {OzOFF_database_path}")
    OzOFF_database = pd.read_parquet(OzOFF_database_path)
    
    logging.info(f"Loading sample data from {sample_path}")
    df_sample_2 = pd.read_parquet(sample_path)

    # Initialize the MatchLipids class with the OzOFF database
    matcher = MatchLipids(OzOFF_database)

    # Run the matching process and save the results
    logging.info(f"Starting matching process for sample: {sample_path}")
    matched_results = MatchLipids.run_match_lipids(matcher, df_sample_2, output_dir)

    # Print results for the current sample
    logging.info(f"Matched Lipid Results for {sample_path}:")
    logging.info(matched_results.head().to_string())
