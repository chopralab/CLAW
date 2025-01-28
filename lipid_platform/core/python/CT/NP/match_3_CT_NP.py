import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import logging
import argparse

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
    def run_match_lipids(matcher, df_sample, output_dir, sample_name):
        """
        Run the lipid matching process and save the results.

        Parameters:
        matcher (MatchLipids): Instance of MatchLipids initialized with the MRM database.
        df_sample (pd.DataFrame): DataFrame containing sample data to match.
        output_dir (str): Directory to save the matched results.
        sample_name (str): Name of the sample for file naming.

        Returns:
        pd.DataFrame: DataFrame containing the matched and filtered results.
        """
        logging.info("Starting lipid matching process.")
        
        # Match lipids for the entire sample DataFrame
        matched_results_df = matcher.match_lipids_parser(df_sample)
        
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Define the output file path
        output_file = os.path.join(output_dir, f"df_match_3_{sample_name}.parquet")
        matched_results_df.to_parquet(output_file, index=False)
        logging.info(f"Matched lipid results saved to {output_file}")
        
        logging.info("Lipid matching process completed.")
        return matched_results_df

def parse_arguments():
    """
    Parse command-line arguments.

    :return: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Match lipids based on ion values.")
    parser.add_argument(
        '--mrm_database',
        type=str,
        required=True,
        help="Path to the MRM database Parquet file."
    )
    parser.add_argument(
        '--sample_file',
        type=str,
        required=True,
        help="Path to the sample Parquet file to be matched."
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help="Directory to save the matched results."
    )
    parser.add_argument(
        '--tolerance',
        type=float,
        default=0.3,
        help="Tolerance value for matching ion values (default: 0.3)."
    )
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Load the MRM database
    logging.info(f"Loading MRM database from {args.mrm_database}")
    mrm_database = pd.read_parquet(args.mrm_database)

    # Load the sample data
    logging.info(f"Loading sample data from {args.sample_file}")
    df_sample = pd.read_parquet(args.sample_file)

    # Extract sample name for output file naming
    sample_name = os.path.splitext(os.path.basename(args.sample_file))[0]

    # Initialize the MatchLipids class with the MRM database and specified tolerance
    matcher = MatchLipids(mrm_database, tolerance=args.tolerance)

    # Run the matching process and save the results
    matched_results = MatchLipids.run_match_lipids(
        matcher=matcher,
        df_sample=df_sample,
        output_dir=args.output_dir,
        sample_name=sample_name
    )

    # Print the first few rows of the matched results
    logging.info(f"Matched Lipid Results for {args.sample_file}:")
    logging.info(matched_results.head().to_string())

if __name__ == "__main__":
    main()
