import pandas as pd
import os
import re
from tqdm import tqdm
import sys

class LipidGrouper:
    def __init__(self, new_columns=None):
        """
        Initialize the LipidGrouper class.

        :param new_columns: Dictionary defining new columns to be created from sample names.
        """
        self.new_columns = new_columns

    def extract_values_from_sample(self, sample):
        """
        Extract specific values from a sample name based on predefined columns.

        :param sample: The sample name from which to extract values.
        :return: Dictionary of extracted values.
        """
        extracted_values = {}
        for col, values in self.new_columns.items():
            extracted_values[col] = next((value for value in values if value in sample), '')
        return extracted_values

    def create_columns_from_sample(self, df):
        """
        Create new columns in the DataFrame based on the sample names.

        :param df: Input DataFrame containing a 'Sample' column.
        :return: DataFrame with new columns extracted from sample names.
        """
        if not self.new_columns:
            return df

        df_copy = df.copy()
        extracted_df = df_copy['Sample'].apply(self.extract_values_from_sample)
        extracted_df = pd.DataFrame(extracted_df.tolist(), index=df_copy.index)

        for col in self.new_columns.keys():
            if col in df_copy.columns:
                df_copy.drop(columns=[col], inplace=True)

        df_copy = pd.concat([df_copy, extracted_df], axis=1)
        return df_copy

    def group_by_ion(self, df):
        """
        Group DataFrame rows by ion information and create a new column for group IDs.

        :param df: Input DataFrame containing 'Parent_Ion', 'Product_Ion', and 'Sample_ID' columns.
        :return: DataFrame with a new 'group_by_ion' column.
        """
        df['group_by_ion'] = df.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()
        return df

    def group_by_lipid(self, df, group_columns=None):
        """
        Group DataFrame rows by lipid information and create a new column for group IDs.

        :param df: Input DataFrame.
        :param group_columns: List of columns to group by. Defaults to ['Lipid', 'Biology', 'Genotype', 'Mouse', 'Cage'].
        :return: DataFrame with a new 'group_by_lipid' column.
        """
        if group_columns is None:
            group_columns = ['Lipid', 'Biology', 'Genotype', 'Mouse', 'Cage']

        df['group_by_lipid'] = df.groupby(group_columns).ngroup()
        return df

    def group_by_func(self, df, group_columns=None):
        """
        Perform a series of grouping operations on the DataFrame and sort by retention time.

        :param df: Input DataFrame.
        :param group_columns: List of columns to group by for the lipid grouping step.
        :return: DataFrame with new group ID columns and sorted by retention time.
        """
        df = self.create_columns_from_sample(df)
        df = self.group_by_ion(df)
        df = self.group_by_lipid(df, group_columns)
        
        # Sort the DataFrame by Retention_Time after grouping
        df = df.sort_values(by=['group_by_lipid', 'Retention_Time'])
        return df

    def save_grouped_results(self, df, file_path):
        """
        Save the grouped DataFrame to a Parquet file.

        :param df: DataFrame to be saved.
        :param file_path: Path where the file will be saved.
        """
        df.to_parquet(file_path, index=False)
        print(f"File successfully saved to: {file_path}")

    def create_folder(self, folder_path):
        """
        Create a folder if it does not already exist.

        :param folder_path: Path of the folder to create.
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Folder created: {folder_path}")
        else:
            print(f"Folder already exists: {folder_path}")

    def extract_matching_lipids(self, OzON_results_df):
        """
        Extract matching lipids from the 'Possible_Lipids' column based on the 'Species' column.

        :param OzON_results_df: Input DataFrame containing 'Possible_Lipids', 'Class', and 'Species' columns.
        :return: DataFrame with extracted matching lipids and their classes.
        """
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
            raise ValueError("No valid matches found in the 'Matches' column.")

        exploded_df[['Lipid', 'Class']] = pd.DataFrame(exploded_df['Matches'].tolist(), index=exploded_df.index)
        matching_results_df = exploded_df[exploded_df['Lipid'].notnull()].copy()
        matching_results_df.drop(columns=['Matches'], inplace=True)

        first_column = matching_results_df.pop('Lipid')
        matching_results_df.insert(0, 'Lipid', first_column)
        possible_lipids_column = matching_results_df.pop('Possible_Lipids')
        matching_results_df['Possible_Lipids'] = possible_lipids_column

        return matching_results_df

if __name__ == "__main__":
    # Ensure the correct number of arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python group_4.py <input_file_path>")
        sys.exit(1)

    # Parse command-line arguments
    input_file_path = sys.argv[1]

    # Initialize LipidGrouper with new columns definition
    grouper = LipidGrouper(new_columns={
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
    })

    # Load the input DataFrame
    OzON_results = pd.read_parquet(input_file_path)

    # Extract matching lipids and perform grouping
    matching_results_df = grouper.extract_matching_lipids(OzON_results)
    matching_results_df['Lipid'] = matching_results_df['Lipid'].fillna('Unknown')
    df_grouped = grouper.group_by_func(matching_results_df)

    # Get the sample value from the DataFrame for naming the output file
    sample_value = df_grouped['Sample'].iloc[0] if 'Sample' in df_grouped.columns else 'unknown_sample'

    # Save the grouped results
    output_file_path = f"Projects/AMP/group/df_group_4_{sample_value}.parquet"
    grouper.create_folder("Projects/AMP/group")
    grouper.save_grouped_results(df_grouped, output_file_path)

    # View the final grouped DataFrame
    print(df_grouped)
