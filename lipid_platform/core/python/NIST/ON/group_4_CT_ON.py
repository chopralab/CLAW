import pandas as pd
import os
import re
from tqdm import tqdm
import sys
import logging
import argparse

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class LipidGrouper:
    def __init__(self, new_columns=None):
        """
        Initialize the LipidGrouper class.

        :param new_columns: Dictionary defining new columns to be created from sample names.
        """
        self.new_columns = new_columns
        logging.debug(f"Initialized LipidGrouper with columns: {new_columns}")

    def create_columns_from_sample(self, df):
        """
        Create new columns in the DataFrame based on the sample names. If the 'STD' column contains 'FAME',
        Biology, Genotype, Cage, and Mouse columns are not added.

        :param df: Input DataFrame containing a 'Sample' and 'STD' column.
        :return: DataFrame with new columns extracted from sample names if 'STD' is not 'FAME'.
        """
        if not self.new_columns:
            logging.debug("No new columns to create from sample names.")
            return df

        logging.debug("Creating new columns from sample names...")
        df_copy = df.copy()
        logging.debug(f"Initial DataFrame columns: {df_copy.columns}")

        # Vectorized check if the 'STD' column contains 'FAME'
        if 'STD' in df_copy.columns and (df_copy['STD'] == 'FAME').any():
            logging.debug("STD is 'FAME'. Skipping Biology, Genotype, Cage, and Mouse columns for relevant rows.")
        else:
            # Only apply extraction for samples where STD is not 'FAME'
            extracted_df = df_copy['Sample'].apply(self.extract_values_from_sample)
            extracted_df = pd.DataFrame(extracted_df.tolist(), index=df_copy.index)

            for col in self.new_columns.keys():
                if col in df_copy.columns:
                    df_copy.drop(columns=[col], inplace=True)

            df_copy = pd.concat([df_copy, extracted_df], axis=1)
            logging.debug(f"Columns after adding extracted values: {df_copy.columns}")

        return df_copy

    def extract_values_from_sample(self, sample):
        """
        Extract specific values from a sample name based on predefined columns.

        :param sample: The sample name from which to extract values.
        :return: Dictionary of extracted values.
        """
        extracted_values = {}
        for col, values in self.new_columns.items():
            extracted_values[col] = next((value for value in values if value in sample), '')
        logging.debug(f"Extracted values from sample '{sample}': {extracted_values}")
        return extracted_values

    def extract_matching_lipids(self, OzON_results_df):
        """
        Extract matching lipids from the 'Possible_Lipids' column based on the 'Species' column.

        :param OzON_results_df: Input DataFrame containing 'Possible_Lipids', 'Class', and 'Species' columns.
        :return: DataFrame with extracted matching lipids and their classes.
        """
        logging.debug("Extracting matching lipids...")
        df_copy = OzON_results_df.copy()
        df_copy.rename(columns={'Lipid': 'Possible_Lipids'}, inplace=True)

        def find_matches(lipid_string, class_string, species):
            if pd.isna(lipid_string):
                return [(None, None)]
            species_pattern = rf'FA\({species}\)_[^|]+'
            lipid_matches = re.findall(species_pattern, lipid_string)

            if lipid_matches:
                lipid_list = [lipid.strip() for lipid in lipid_string.split(' | ')]
                class_list = [cls.strip() for cls in class_string.split(' | ')]
                matched_lipids_classes = [(lipid, class_list[lipid_list.index(lipid.strip())]) for lipid in lipid_matches]
                return matched_lipids_classes
            return [(None, None)]

        df_copy['Matches'] = df_copy.apply(
            lambda row: find_matches(row['Possible_Lipids'], row['Class'], row['Species']),
            axis=1
        )

        exploded_df = df_copy.explode('Matches')
        valid_matches = exploded_df['Matches'].dropna()
        if valid_matches.empty:
            logging.error("No valid matches found in the 'Matches' column.")
            raise ValueError("No valid matches found in the 'Matches' column.")

        exploded_df[['Lipid', 'Class']] = pd.DataFrame(exploded_df['Matches'].tolist(), index=exploded_df.index)
        matching_results_df = exploded_df[exploded_df['Lipid'].notnull()].copy()
        matching_results_df.drop(columns=['Matches'], inplace=True)

        first_column = matching_results_df.pop('Lipid')
        matching_results_df.insert(0, 'Lipid', first_column)
        possible_lipids_column = matching_results_df.pop('Possible_Lipids')
        matching_results_df['Possible_Lipids'] = possible_lipids_column

        logging.debug("Matching lipids extracted successfully.")
        return matching_results_df

    def group_by_ion(self, df):
        """
        Group DataFrame rows by ion information and create a new column for group IDs.

        :param df: Input DataFrame containing 'Parent_Ion', 'Product_Ion', and 'Sample_ID' columns.
        :return: DataFrame with a new 'group_by_ion' column.
        """
        logging.debug("Grouping by ion...")
        df['group_by_ion'] = df.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()
        logging.debug("DataFrame after grouping by ion:\n%s", df.head())
        return df

    def group_by_lipid(self, df, group_columns=None):
        """
        Group DataFrame rows by lipid information and create a new column for group IDs.

        :param df: Input DataFrame.
        :param group_columns: List of columns to group by. Defaults to ['Lipid', 'Biology', 'Genotype', 'Mouse', 'Cage'].
        :return: DataFrame with a new 'group_by_lipid' column.
        """
        # Check if any row in the 'STD' column is 'FAME'
        if 'STD' in df.columns and (df['STD'] == 'FAME').any():
            logging.debug("STD is 'FAME'. Skipping Biology, Genotype, Cage, and Mouse columns for grouping.")
            group_columns = ['Lipid']  # Only group by 'Lipid' if STD is 'FAME'
        else:
            if group_columns is None:
                group_columns = ['Lipid', 'Biology', 'Genotype', 'Mouse', 'Cage']

        logging.debug(f"Grouping by columns: {group_columns}")
        df['group_by_lipid'] = df.groupby(group_columns).ngroup()
        logging.debug("DataFrame after grouping by lipid:\n%s", df.head())
        return df

    def group_by_func(self, df, group_columns=None):
        """
        Perform a series of grouping operations on the DataFrame and sort by retention time.

        :param df: Input DataFrame.
        :param group_columns: List of columns to group by for the lipid grouping step.
        :return: DataFrame with new group ID columns and sorted by retention time.
        """
        logging.debug("Starting grouping operations...")
        df = self.create_columns_from_sample(df)
        df = self.group_by_ion(df)
        df = self.group_by_lipid(df, group_columns)

        # Sort the DataFrame by Retention_Time after grouping
        logging.debug("Sorting by Retention_Time...")
        df = df.sort_values(by=['group_by_lipid', 'Retention_Time'])
        logging.debug("DataFrame after sorting:\n%s", df.head())
        return df

    def save_grouped_results(self, df, file_path):
        """
        Save the grouped DataFrame to a Parquet file.

        :param df: DataFrame to be saved.
        :param file_path: Path where the file will be saved.
        """
        logging.debug(f"Saving grouped DataFrame to {file_path}...")
        df.to_parquet(file_path, index=False)
        logging.debug(f"File successfully saved to: {file_path}")

    def create_folder(self, folder_path):
        """
        Create a folder if it does not already exist.

        :param folder_path: Path of the folder to create.
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            logging.debug(f"Folder created: {folder_path}")
        else:
            logging.debug(f"Folder already exists: {folder_path}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Group Lipid Data by Various Criteria.")
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the input Parquet file."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory where the output Parquet file will be saved."
    )
    parser.add_argument(
        "--python_file",
        type=str,
        default="group_4_CT_ON.py",
        help="Name of the Python script (default: group_4_CT_ON.py)."
    )
    parser.add_argument(
        "--log_level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO)."
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    # Update logging level based on arguments
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        logging.error(f'Invalid log level: {args.log_level}')
        sys.exit(1)
    logging.getLogger().setLevel(numeric_level)

    # Initialize LipidGrouper with new columns definition
    grouper = LipidGrouper(new_columns={
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
    })

    # Load the input DataFrame
    logging.info(f"Loading input DataFrame from: {args.input_file}")
    try:
        OzON_results = pd.read_parquet(args.input_file)
        logging.debug(f"Input DataFrame loaded. Columns: {OzON_results.columns}")
    except Exception as e:
        logging.error(f"Failed to load input file: {e}")
        sys.exit(1)

    # Extract matching lipids
    try:
        matching_results_df = grouper.extract_matching_lipids(OzON_results)
    except Exception as e:
        logging.error(f"Error during lipid extraction: {e}")
        sys.exit(1)

    # Group the lipids
    try:
        df_grouped = grouper.group_by_func(matching_results_df)
    except Exception as e:
        logging.error(f"Error during grouping: {e}")
        sys.exit(1)

    # Get the sample value from the DataFrame for naming the output file
    sample_value = df_grouped['Sample'].iloc[0] if 'Sample' in df_grouped.columns else 'unknown_sample'

    # Define output file path
    output_file_path = os.path.join(args.output_dir, f"df_group_4_{sample_value}.parquet")
    
    # Create the output directory if it doesn't exist
    grouper.create_folder(args.output_dir)

    # Save the grouped results
    try:
        grouper.save_grouped_results(df_grouped, output_file_path)
    except Exception as e:
        logging.error(f"Failed to save grouped results: {e}")
        sys.exit(1)

    # View the final grouped DataFrame
    logging.info(f"Final grouped DataFrame saved to {output_file_path}")
    logging.debug(f"Final grouped DataFrame preview:\n{df_grouped.head()}")
