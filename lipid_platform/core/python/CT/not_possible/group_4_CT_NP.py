import pandas as pd
import os
import re
import sys
import argparse
from tqdm import tqdm

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
        Ignores any additional characters following the main FA(##:#) structure.
        """
        if pd.isna(lipid) or lipid.strip() == "":
            species = f"Unknown_{self.unknown_count}"
            self.unknown_count += 1
            return species  # Assign a unique Unknown species
        if lipid.startswith('FA(d2-'):
            lipid = lipid.replace('d2-', '')  # Remove 'd2-' prefix

        # Match pattern for FA(##:#), ignoring anything after it
        match = re.search(r'FA\((\d+:\d+)\)', lipid)
        if match:
            return match.group(1)  # Returns just the core FA(##:#) part
        else:
            print(f"Warning: Could not extract species from lipid: {lipid}", file=sys.stderr)
            species = f"Unknown_{self.unknown_count}"
            self.unknown_count += 1
            return species

    def species_create(self, df):
        """
        Create a Species column from the Lipid column, extracting entries by '|' and keeping all species in a single row.
        """
        print("Creating Species column from 'Lipid' column...", file=sys.stderr)
        df['Lipid_list'] = df['Lipid'].astype(str).str.split('|')

        # Extract species for each lipid part and combine them in a single entry
        df['Species'] = df['Lipid_list'].apply(
            lambda lipids: '|'.join([self.extract_species(lipid) for lipid in lipids])
        )

        # Drop the helper column
        df.drop(columns=['Lipid_list'], inplace=True)

        print("Species column created successfully.", file=sys.stderr)
        return df

    def extract_values_from_sample(self, sample):
        """
        Extract specific values from a sample name based on predefined columns.
        """
        extracted_values = {}
        for col, values in self.new_columns.items():
            extracted_values[col] = next((value for value in values if value in sample), 'Unknown')
        return extracted_values

    def create_columns_from_sample(self, df):
        """
        Create new columns in the DataFrame based on the sample names.
        """
        print("Creating new columns from 'Sample' column...", file=sys.stderr)
        extracted_df = df['Sample'].apply(self.extract_values_from_sample).apply(pd.Series)

        # Replace missing values with 'Unknown' if necessary
        extracted_df.fillna('Unknown', inplace=True)

        df = pd.concat([df, extracted_df], axis=1)
        return df

    def group_by_ion(self, df):
        """
        Group DataFrame rows by ion information and create a new column for group IDs.
        """
        print("Grouping by ion information...", file=sys.stderr)
        df['group_by_ion'] = df.groupby(['Parent_Ion', 'Product_Ion', 'Sample_ID']).ngroup()
        return df

    def group_by_lipid(self, df, group_columns):
        """
        Group DataFrame rows by lipid information and create a new column for group IDs.
        """
        print(f"Grouping by lipid information using columns: {group_columns}", file=sys.stderr)
        df['group_by_lipid'] = df.groupby(group_columns).ngroup()
        return df

    def group_by_func(self, df, group_columns=None, STD_Only=None):
        """
        Perform a series of grouping operations on the DataFrame and sort by retention time.
        If STD_Only is set to 'STD', it will skip grouping by certain columns.
        """
        print("Starting grouping process...", file=sys.stderr)

        # Determine grouping columns based on STD_Only flag and content
        if ('Sample' in df.columns and df['Sample'].str.contains('FAME').any()) or \
           ('Std' in df.columns and df['Std'].str.contains('FAME').any()):
            print("Detected 'FAME' in 'Sample' or 'Std' column. Using simplified grouping.", file=sys.stderr)
            group_columns = ['Lipid']
        else:
            if STD_Only == 'STD':
                print("STD_Only is set. Using only 'Lipid' for grouping.", file=sys.stderr)
                group_columns = ['Lipid']
            else:
                if group_columns is None:
                    group_columns = ['Lipid', 'Biology', 'Genotype', 'Mouse', 'Cage']
                group_columns = [col for col in group_columns if col in df.columns]
                print(f"Using group columns: {group_columns}", file=sys.stderr)

        df = self.create_columns_from_sample(df)
        df = self.group_by_ion(df)
        df = self.group_by_lipid(df, group_columns)

        df.sort_values(by=['group_by_lipid', 'Retention_Time'], inplace=True)

        if STD_Only == 'STD':
            columns_to_drop = ['Biology', 'Genotype', 'Cage', 'Mouse']
            df.drop(columns=[col for col in columns_to_drop if col in df.columns], inplace=True)
            print(f"Dropped columns {columns_to_drop} as STD_Only is set.", file=sys.stderr)

        print("Grouping and sorting completed.", file=sys.stderr)
        return df

    def save_grouped_results(self, df, output_dir):
        """
        Save the grouped DataFrame to a Parquet file in the specified output directory.
        """
        sample_value = df['Sample'].iloc[0] if not df.empty and 'Sample' in df.columns else 'unknown_sample'
        output_file_path = os.path.join(output_dir, f"df_grouped_{sample_value}_OFF.parquet")
        
        print(f"Saving DataFrame to {output_file_path}...", file=sys.stderr)
        df.to_parquet(output_file_path, index=False)
        print(f"File saved to {output_file_path}.", file=sys.stderr)

    def create_folder(self, folder_path):
        """
        Create a folder if it does not already exist.
        """
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Created folder: {folder_path}", file=sys.stderr)
        else:
            print(f"Folder already exists: {folder_path}", file=sys.stderr)


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Group lipid data based on specified parameters.")
    parser.add_argument(
        "--input_file",
        required=True,
        help="Path to the input Parquet file."
    )
    parser.add_argument(
        "--std_only",
        required=True,
        choices=['STD', 'Sample'],
        help="Flag to indicate STD_ONLY mode ('STD' or 'Sample')."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where the output Parquet file will be saved."
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Initialize LipidGrouper with new columns definition
    grouper = LipidGrouper(new_columns={
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
    })

    # Load the input DataFrame
    print(f"Loading input DataFrame from file: {args.input_file}", file=sys.stderr)
    try:
        OzON_results = pd.read_parquet(args.input_file)
    except Exception as e:
        print(f"Error loading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Convert relevant columns to strings
    OzON_results['Lipid'] = OzON_results['Lipid'].astype(str)
    OzON_results['Sample'] = OzON_results['Sample'].astype(str)
    print("DataFrame loaded successfully. Here's a preview:", file=sys.stderr)
    print(OzON_results.head(), file=sys.stderr)

    # Create Species column from the Lipid column
    OzON_results = grouper.species_create(OzON_results)

    # Perform grouping
    df_grouped = grouper.group_by_func(OzON_results, STD_Only=args.std_only)

    # Create output directory if it doesn't exist
    grouper.create_folder(args.output_dir)

    # Save the grouped results
    grouper.save_grouped_results(df_grouped, args.output_dir)

    print("Final grouped DataFrame preview:", file=sys.stderr)
    print(df_grouped.head(), file=sys.stderr)


if __name__ == "__main__":
    main()
