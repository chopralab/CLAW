import pandas as pd
from tqdm import tqdm
import os

class SampleIDExtract:
    def __init__(self, new_columns=None):
        """
        Initialize the SampleIDExtract class with predefined column mappings.

        Input:
        - new_columns: A dictionary containing lists of possible values for each category.

        Output:
        - An instance of the SampleIDExtract class with the new_columns attribute set.
        """
        if new_columns is None:
            new_columns = {
                'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
                'Genotype': ['5xFAD', 'WT'],
                'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
                'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
                'Other': ['Blank', 'blank']
            }
        self.new_columns = new_columns

    def extract_sample_parts(self, sample_id):
        """
        Extract parts of the sample ID that match predefined categories.

        Input:
        - sample_id: A string representing the sample ID.

        Output:
        - A string with matched parts joined by underscores.
        """
        parts = sample_id.replace('-', '_').split('_')
        matched_parts = []
        added_parts = set()

        for part in parts:
            for key in self.new_columns:
                for value in self.new_columns[key]:
                    if value in part and part not in added_parts:
                        matched_parts.append(value)
                        added_parts.add(part)
                        break

        return '_'.join(matched_parts)

    def apply_extraction(self, df):
        """
        Apply the sample ID extraction to a DataFrame.

        Input:
        - df: A pandas DataFrame with a 'Sample_ID' column.

        Output:
        - The input DataFrame with an additional 'Sample' column.
        """
        tqdm.pandas(desc="Extracting Sample Parts")
        df['Sample'] = df['Sample_ID'].progress_apply(self.extract_sample_parts)
        return df

def main():
    """
    Main function to load a Parquet file, apply sample ID extraction, and save separate Parquet files for each sample type.

    Input:
    - None

    Output:
    - Separate Parquet files for each sample type, named based on the 'Sample' column, saved in 'Projects/AMP/samples/'.
    """
    new_columns = {
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
        'Other': ['Blank', 'blank']
    }

    output_dir = 'Projects/AMP/samples/'
    os.makedirs(output_dir, exist_ok=True)

    print("Loading Parquet file...")
    OzESI_df = pd.read_parquet("df_mzml_parser_1.parquet")
    print("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(new_columns)

    print("Applying extraction...")
    OzON_Data = sample_extractor.apply_extraction(OzESI_df)
    print("Extraction applied successfully.")

    unique_samples = OzON_Data['Sample'].unique()
    for sample in tqdm(unique_samples, desc="Processing Samples"):
        print(f"Processing sample: {sample}")
        sample_df = OzON_Data[OzON_Data['Sample'] == sample]
        filename = os.path.join(output_dir, f"df_sample_2_{sample}.parquet")
        print(f"Saving {filename}...")
        sample_df.to_parquet(filename, index=False, compression="brotli")
        print(f"File {filename} saved successfully.")

if __name__ == "__main__":
    main()
