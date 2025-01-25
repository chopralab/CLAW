#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm

# Setup logging
logging.basicConfig(
    level=logging.INFO,  # Set to DEBUG for more detailed logs
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

class MatchLipids:
    def __init__(self, mrm_database: pd.DataFrame, tolerance: float = 0.3):
        self.mrm_database = mrm_database.copy()  # Work with a copy to prevent altering original data
        self.tolerance = tolerance

        # Round the ion values in the MRM database for easier comparison
        self.mrm_database['Parent_Ion'] = np.round(self.mrm_database['Parent_Ion'], 1)
        self.mrm_database['Product_Ion'] = np.round(self.mrm_database['Product_Ion'], 1)

    def within_tolerance(self, values1: np.ndarray, value2: float) -> np.ndarray:
        """Check if elements in values1 are within tolerance of value2."""
        return np.abs(values1 - value2) <= self.tolerance

    def match_lipids_parser(self, df: pd.DataFrame) -> pd.DataFrame:
        """Match lipids based on Parent_Ion and Product_Ion within the specified tolerance."""
        logging.info("Starting lipid matching parser.")

        # Round the ion values in the input DataFrame for easier comparison
        df = df.copy()
        df['Parent_Ion'] = np.round(df['Parent_Ion'], 1)
        df['Product_Ion'] = np.round(df['Product_Ion'], 1)
        df['OzESI_Intensity'] = np.round(df['OzESI_Intensity'], 0)

        matched_lipids = []
        matched_classes = []

        for _, row in tqdm(df.iterrows(), total=len(df), desc="Matching Lipids"):
            parent_matches = self.within_tolerance(self.mrm_database['Parent_Ion'].values, row['Parent_Ion'])
            product_matches = self.within_tolerance(self.mrm_database['Product_Ion'].values, row['Product_Ion'])
            matches = parent_matches & product_matches

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
    def find_matching_OzOFF_file(OzOFF_dir: str, sample_value: str) -> str:
        """Find the OzOFF file that matches the sample value."""
        logging.info(f"Searching for matching OzOFF file in directory: {OzOFF_dir}")
        try:
            files = os.listdir(OzOFF_dir)
        except FileNotFoundError:
            logging.error(f"OzOFF directory not found: {OzOFF_dir}")
            sys.exit(1)

        for file in files:
            if sample_value in file and file.endswith('.parquet'):
                matching_file = os.path.join(OzOFF_dir, file)
                logging.info(f"Found matching OzOFF file: {matching_file}")
                return matching_file

        logging.error(f"No matching OzOFF file found for Sample: {sample_value}")
        sys.exit(1)

    @staticmethod
    def run_match_lipids(
        OzOFF_dir: str,
        OzON_database_path: str,
        sample_path: str,
        tolerance: float,
        retention_time_window: float,
        output_dir: str
    ) -> pd.DataFrame:
        """Run the lipid matching process."""
        logging.info("Starting lipid matching process.")

        # Load sample data
        logging.info(f"Loading sample data from: {sample_path}")
        try:
            df_sample = pd.read_parquet(sample_path)
        except Exception as e:
            logging.error(f"Failed to load sample file: {e}")
            sys.exit(1)

        if 'Sample' not in df_sample.columns:
            logging.error("Sample file must contain a 'Sample' column.")
            sys.exit(1)

        sample_value = df_sample['Sample'].iloc[0]
        logging.info(f"Sample value: {sample_value}")

        # Find matching OzOFF file
        matching_OzOFF_file = MatchLipids.find_matching_OzOFF_file(OzOFF_dir, sample_value)

        # Load OzOFF database
        logging.info(f"Loading OzOFF database from: {matching_OzOFF_file}")
        try:
            OzOFF_database = pd.read_parquet(matching_OzOFF_file)
        except Exception as e:
            logging.error(f"Failed to load OzOFF database: {e}")
            sys.exit(1)

        # Load OzON database
        logging.info(f"Loading OzON database from: {OzON_database_path}")
        try:
            OzON_database = pd.read_parquet(OzON_database_path)
        except Exception as e:
            logging.error(f"Failed to load OzON database: {e}")
            sys.exit(1)

        # Calculate STD_RT_Dif
        OzOFF_database['STD_RT_Dif'] = OzOFF_database['STD_RT_OFF'] - df_sample['STD_RT_ON'].iloc[0]
        logging.debug("Calculated STD_RT_Dif for OzOFF database.")

        # Log STD_RT values for debugging
        for _, row in OzOFF_database.iterrows():
            logging.debug(
                f"STD_RT_OFF: {row['STD_RT_OFF']}, "
                f"STD_RT_ON: {df_sample['STD_RT_ON'].iloc[0]}, "
                f"STD_RT_Dif: {row['STD_RT_Dif']}"
            )

        # Add STD_RT_Dif to sample DataFrame
        df_sample['STD_RT_Dif'] = OzOFF_database['STD_RT_Dif'].values[0]  # Assuming single STD_RT_ON

        # Filter OzON database based on Species
        species_list = OzOFF_database['Species'].unique().tolist()
        logging.info(f"Filtering OzON database for species: {species_list}")
        temp_OzON_database = OzON_database[OzON_database['Species'].isin(species_list)].copy()
        logging.debug(f"Filtered OzON database contains {len(temp_OzON_database)} entries.")

        # Initialize matcher
        matcher = MatchLipids(temp_OzON_database, tolerance=tolerance)

        OzON_results_df = pd.DataFrame()

        # Iterate through OzOFF database entries
        for index, row in tqdm(OzOFF_database.iterrows(), total=len(OzOFF_database), desc="Processing OzOFF Entries"):
            species = row['Species']
            retention_time = row['Retention_Time']
            std_rt_dif = row['STD_RT_Dif']
            std_rt_off = row['STD_RT_OFF']

            logging.debug(
                f"Processing Species: {species}, "
                f"Retention_Time: {retention_time}, "
                f"STD_RT_Dif: {std_rt_dif}, "
                f"STD_RT_OFF: {std_rt_off}"
            )

            # Define retention time window
            adjusted_rt = retention_time + std_rt_dif
            rt_start = adjusted_rt - retention_time_window
            rt_end = adjusted_rt + retention_time_window

            # Filter sample data based on Sample and Retention Time window
            filtered_sample = df_sample[
                (df_sample['Sample'] == sample_value) &
                (df_sample['Retention_Time'] + df_sample['STD_RT_Dif'] >= rt_start) &
                (df_sample['Retention_Time'] + df_sample['STD_RT_Dif'] <= rt_end)
            ].copy()

            if filtered_sample.empty:
                logging.info(
                    f"No matching entries in sample for Species: {species} within RT window ({rt_start} - {rt_end})."
                )
                continue

            logging.info(
                f"Found {len(filtered_sample)} matching entries for Species: {species} within RT window."
            )

            # Assign Species to the filtered data
            filtered_sample['Species'] = species

            # Perform lipid matching
            matched_data = matcher.match_lipids_parser(filtered_sample)
            matched_data['STD_RT_OFF'] = std_rt_off
            matched_data['STD_RT_ON'] = df_sample['STD_RT_ON'].iloc[0]
            matched_data['STD_RT_Dif'] = std_rt_dif

            OzON_results_df = pd.concat([OzON_results_df, matched_data], ignore_index=True)

        # Log final DataFrame columns
        logging.debug(f"Final results DataFrame columns: {OzON_results_df.columns.tolist()}")

        # Save results if output directory is provided
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"df_match_3_{sample_value}.parquet")
            try:
                OzON_results_df.to_parquet(output_file, index=False)
                logging.info(f"Results saved to: {output_file}")
            except Exception as e:
                logging.error(f"Failed to save results: {e}")
                sys.exit(1)

        logging.info("Lipid matching process completed successfully.")
        return OzON_results_df

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Match lipids between OzOFF and OzON databases based on ion values and retention times."
    )
    parser.add_argument(
        "--ozoff_dir",
        required=True,
        help="Path to the OzOFF directory containing parquet files."
    )
    parser.add_argument(
        "--ozon_database",
        required=True,
        help="Path to the OzON database parquet file."
    )
    parser.add_argument(
        "--sample_file",
        required=True,
        help="Path to the sample parquet file."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save the output parquet file."
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.3,
        help="Tolerance for ion matching (default: 0.3)."
    )
    parser.add_argument(
        "--retention_time_window",
        type=float,
        default=0.5,
        help="Retention time window for matching (default: 0.5)."
    )
    parser.add_argument(
        "--log_level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO)."
    )
    return parser.parse_args()

def main():
    """Main function to execute lipid matching."""
    args = parse_arguments()

    # Update logging level based on user input
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        logging.error(f"Invalid log level: {args.log_level}")
        sys.exit(1)
    logging.getLogger().setLevel(numeric_level)

    # Run lipid matching
    results_df = MatchLipids.run_match_lipids(
        OzOFF_dir=args.ozoff_dir,
        OzON_database_path=args.ozon_database,
        sample_path=args.sample_file,
        tolerance=args.tolerance,
        retention_time_window=args.retention_time_window,
        output_dir=args.output_dir
    )

    # Optionally, perform further analysis or export additional results here
    logging.info("Process completed successfully.")

if __name__ == "__main__":
    main()
