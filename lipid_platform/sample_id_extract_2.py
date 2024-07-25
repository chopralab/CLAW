import pandas as pd
from tqdm import tqdm

class SampleIDExtract:
    def __init__(self, new_columns=None):
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
        tqdm.pandas(desc="Extracting Sample Parts")
        df['Sample'] = df['Sample_ID'].progress_apply(self.extract_sample_parts)
        return df

def main():
    new_columns = {
        'Biology': ['cortex', 'dienc', 'hippo', 'cereb'],
        'Genotype': ['5xFAD', 'WT'],
        'Cage': ['FAD231', 'FAD259', 'FAD257', 'FAD263', 'FAD249', 'FAD246', 'FAD245'],
        'Mouse': ['m1', 'm2', 'm3', 'm4', 'm5'],
        'Other': ['Blank', 'blank']
    }

    print("Loading Parquet file...")
    OzESI_df = pd.read_parquet("OzON_Data_2.parquet")
    print("Parquet file loaded successfully.")

    sample_extractor = SampleIDExtract(new_columns)

    print("Applying extraction...")
    OzON_Data = sample_extractor.apply_extraction(OzESI_df)
    print("Extraction applied successfully.")

    print("Saving the result to OzON_Data_SIE_3.parquet...")
    OzON_Data.to_parquet("OzON_Data_SIE_3.parquet", index=False, compression="brotli")
    print("File saved successfully.")

if __name__ == "__main__":
    main()
