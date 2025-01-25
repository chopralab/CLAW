import pandas as pd
import os
import re
from tqdm import tqdm
import sys

class LipidGrouper:
    def __init__(self, new_columns=None):
        """
        Initialize the LipidGrouper class.
        """
        print("Initializing LipidGrouper...", file=sys.stderr)
        self.new_columns = new_columns
        self.unknown_count = 1  # Initialize counter for unknown species

    def extract_species(self, lipid):
        """
        Extract the species from the Lipid column using the FA(##:#) pattern.
        Special handling for cases like FA(d2-16:0).
        """
        if pd.isna(lipid) or lipid.strip() == "":
            return f"Unknown_{self.unknown_count}"  # Assign a unique Unknown species
        if lipid.startswith('FA(d2-'):
            lipid = lipid.replace('d2-', '')  # Remove 'd2-' prefix
        match = re.search(r'FA\((\d+:\d+)\)', lipid)  # Looking for FA(##:#) pattern
        if match:
            return match.group(1)
        else:
            print(f"Warning: Could not extract species from lipid: {lipid}", file=sys.stderr)
            return f"Unknown_{self.unknown_count}"

    def species_create(self, df):
        """
        Create a Species column from the Lipid column.
        """
        print(f"Checking 'Lipid' column for non-empty values. Sample values from the 'Lipid' column:", file=sys.stderr)
        print(df['Lipid'].unique(), file=sys.stderr)  # Print unique values for debugging

        # Convert Lipid column to string explicitly
        df['Lipid'] = df['Lipid'].astype(str)
        
        # Apply the species extraction function
        df['Species'] = df['Lipid'].apply(lambda lipid: self.extract_species(lipid))

        # Update the unknown counter for each row that has 'Unknown' species
        df['Species'] = df.apply(lambda row: self.increment_unknown_count(row), axis=1)

        # Print some rows of the Lipid and Species columns for further debugging
        print("Species column created. Sample output:", file=sys.stderr)
        print(df[['Lipid', 'Species']].head(), file=sys.stderr)
        
        return df  # Ensure the modified DataFrame is returned

    def increment_unknown_count(self, row):
        if 'Unknown' in row['Species']:
            self.unknown_count += 1
        return row['Species']

    def extract_values_from_sample(self, sample):
        """
        Extract specific values from a sample name based on predefined columns.
        """
        print(f"Extracting values from sample name: {sample}", file=sys.stderr)
        extracted_values = {}
        for col, values in self.new_columns.items():
            extracted_values[col] = next((value for value in values if value in sample), '')
        print(f"Extracted values: {extracted_values}", file=sys.stderr)
        return extracted_values

    def create_columns_from_sample(self, df):
        """
        Create new columns in the DataFrame based on the sample names.
        """
        print(f"Creating new columns from 'Sample' column for DataFrame with {len(df)} rows.", file=sys.stderr)
        df_copy = df.copy()

        extracted_df = df_copy['Sample'].apply(self.extract_values_from_sample)
        extracted_df = pd.DataFrame(extracted_df.tolist(), index=df_copy.index)

        for col in self.new_columns.keys():
            if col in df_copy.columns:
                df_copy.drop(columns=[col], inplace=True)

        df_copy = pd.concat([df_copy, extracted_df], axis=1)
        print("New columns created from Sample. DataFrame now has the following columns:", file=sys.stderr)
        print(df_copy.columns, file=sys.stderr)
        return df_copy

    def group_by_ion(self, df):
        """
        Group DataFrame rows by ion information and create a new column for group IDs.
        """
        print(f"Grouping by ion information for DataFrame with {len(df)} rows.", file=sys.stderr)
        df['group_by_ion'] = df.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()
        print("Grouping by ion complete. Sample output:", file=sys.stderr)
        print(df[['Parent_Ion', 'Product_Ion', 'Sample_ID', 'group_by_ion']].head(), file=sys.stderr)
        return df

    def group_by_lipid(self, df, group_columns=None):
        """
        Group DataFrame rows by lipid information and create a new column for group IDs.
        """
        print(f"Grouping by lipid information using columns: {group_columns}", file=sys.stderr)
        df['group_by_lipid'] = df.groupby(group_columns).ngroup()
        print("Grouping by lipid complete. Sample output:", file=sys.stderr)
        print(df[['Lipid', 'group_by_lipid']].head(), file=sys.stderr)
        return df

    def group_by_func(self, df, group_columns=None, STD_Only=None):
        """
        Perform a series of grouping operations on the DataFrame and sort by retention time.
        If STD_Only is set to 'STD', it will skip grouping by certain columns.
        """
        print(f"Starting grouping process for DataFrame with {len(df)} rows.", file=sys.stderr)

        # Check if 'Sample' or 'Std' column contains 'FAME'
        if 'Sample' in df.columns and df['Sample'].str.contains('FAME').any():
            print("Sample column contains 'FAME'. Using simplified grouping.", file=sys.stderr)
            group_columns = ['Lipid']
        elif 'Std' in df.columns and df['Std'].str.contains('FAME').any():
            print("Std column contains 'FAME'. Using simplified grouping.", file=sys.stderr)
            group_columns = ['Lipid']
        else:
            if STD_Only == 'STD':
                print(f"STD_Only is set to '{STD_Only}', skipping Biology, Genotype, Mouse, and Cage columns.", file=sys.stderr)
                group_columns = ['Lipid']  # Only group by Lipid if STD_Only is 'STD'
            else:
                if group_columns is None:
                    group_columns = ['Lipid', 'Biology', 'Genotype', 'Mouse', 'Cage']
                group_columns = [col for col in group_columns if col in df.columns]
                print(f"Using group columns: {group_columns}", file=sys.stderr)

        df = self.create_columns_from_sample(df)
        df = self.group_by_ion(df)
        df = self.group_by_lipid(df, group_columns)

        df = df.sort_values(by=['group_by_lipid', 'Retention_Time'])

        if STD_Only == 'STD':
            # Drop the unnecessary columns if STD_Only is 'STD'
            columns_to_drop = ['Biology', 'Genotype', 'Cage', 'Mouse']
            df = df.drop(columns=[col for col in columns_to_drop if col in df.columns])
            print(f"Columns {columns_to_drop} dropped as STD_Only is set to 'STD'.", file=sys.stderr)

        print("DataFrame sorted by 'group_by_lipid' and 'Retention_Time'.", file=sys.stderr)
        return df



    def save_grouped_results(self, df, file_path):
        """
        Save the grouped DataFrame to a Parquet file.
        """
        print(f"Saving DataFrame with {len(df)} rows to file: {file_path}", file=sys.stderr)
        df.to_parquet(file_path, index=False)
        print(f"File successfully saved to: {file_path}", file=sys.stderr)

    def create_folder(self, folder_path):
        """
        Create a folder if it does not already exist.
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Folder created: {folder_path}", file=sys.stderr)
        else:
            print(f"Folder already exists: {folder_path}", file=sys.stderr)


if __name__ == "__main__":
    # Ensure the correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python group_4.py <input_file_path> <STD_Only>", file=sys.stderr)
        sys.exit(1)

    # Parse command-line arguments
    input_file_path = sys.argv[1]
    STD_Only = sys.argv[2]

    # Initialize LipidGrouper with new columns definition
    grouper = LipidGrouper(new_columns={
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
    })

    # Load the input DataFrame
    print(f"Loading input DataFrame from file: {input_file_path}", file=sys.stderr)
    OzON_results = pd.read_parquet(input_file_path)

    # Convert relevant columns to strings
    OzON_results['Lipid'] = OzON_results['Lipid'].astype(str)
    OzON_results['Sample'] = OzON_results['Sample'].astype(str)
    print(f"DataFrame loaded. Head of the DataFrame:", file=sys.stderr)
    print(OzON_results.head(), file=sys.stderr)

    # Create Species column from the Lipid column
    OzON_results = grouper.species_create(OzON_results)

    # Perform grouping
    df_grouped = grouper.group_by_func(OzON_results, STD_Only=STD_Only)

    # Handle empty DataFrame before accessing Sample
    if df_grouped.empty:
        print("Warning: df_grouped is empty. Assigning 'unknown_sample'.", file=sys.stderr)
        sample_value = 'unknown_sample'
    else:
        sample_value = df_grouped['Sample'].iloc[0] if 'Sample' in df_grouped.columns else 'unknown_sample'

    # Save the grouped results
    output_file_path = f"Projects/AMP/group/OFF/df_group_4_{sample_value}_OFF.parquet"
    grouper.create_folder("Projects/AMP/group/OFF")
    grouper.save_grouped_results(df_grouped, output_file_path)

    # View the final grouped DataFrame
    print(f"Final grouped DataFrame:", file=sys.stderr)
    print(df_grouped.head(), file=sys.stderr)
