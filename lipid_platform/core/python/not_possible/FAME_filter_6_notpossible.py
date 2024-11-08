# core/python/fame_filter_AMP_7.py

import pandas as pd
import os
import re
from tqdm import tqdm

def filter_OzOFF_by_fame_std(input_dir, fame_std, output_dir, fame_rt_window=0.5):
    """
    Filters parquet files in the input directory based on fame_std retention times and species,
    using Adjusted_RT_FAME from the input files to compare against fame_std Retention_Time.
    Updates Species to only contain the matched species and saves the filtered result as parquet files
    in the output directory. Logs and removed lipids CSVs are saved in a log directory inside the output directory.
    """

    def filter_fa_species(df, species_pattern):
        """
        Filters DataFrame rows based on the provided fatty acid species pattern (e.g., 'FA(16:1)').
        Keeps only entries that match the pattern in the Lipid column.
        """
        pattern = rf'\b{species_pattern}(?:_[^|]*)?'
        filtered_df = df[df['Lipid'].str.contains(pattern, regex=True)]
        filtered_df['Lipid'] = filtered_df['Lipid'].apply(lambda x: '|'.join(re.findall(pattern, x)))
        return filtered_df

    print("Starting filter_OzOFF_by_fame_std function")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Retention time window: {fame_rt_window}")
    print(f"Fame standard DataFrame columns: {fame_std.columns.tolist()}")
    print(f"Total fame_std entries: {len(fame_std)}")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Create a log subdirectory inside the output directory
    log_dir = os.path.join(output_dir, 'log')
    os.makedirs(log_dir, exist_ok=True)

    # Check if input directory has parquet files
    input_files = [f for f in os.listdir(input_dir) if f.endswith(".parquet")]
    if not input_files:
        print("No parquet files found in input directory. Exiting function.")
        return

    # Iterate over all parquet files in the input directory with progress tracking
    for parquet_file in tqdm(input_files, desc="Processing parquet files"):
        file_path = os.path.join(input_dir, parquet_file)
        filtered_ozON = pd.read_parquet(file_path)
        filtered_ozON['Lipid_Possible'] = filtered_ozON['Lipid']

        # Detailed logging information for input data
        print(f"\nProcessing file: {parquet_file}")
        print(f"Number of entries in filtered_ozON: {len(filtered_ozON)}")
        print(f"Columns in filtered_ozON: {filtered_ozON.columns.tolist()}")

        # Filter based on all FA species patterns from fame_std
        all_filtered_dataframes = []
        for _, fame_entry in fame_std.iterrows():
            species_pattern = fame_entry['Species']
            print(f"Filtering by species pattern: {species_pattern}")
            filtered_species_df = filter_fa_species(filtered_ozON, species_pattern)
            all_filtered_dataframes.append(filtered_species_df)
            print(f"Number of entries after filtering by {species_pattern}: {len(filtered_species_df)}")

        # Concatenate all filtered DataFrames
        filtered_ozON = pd.concat(all_filtered_dataframes).drop_duplicates().reset_index(drop=True)
        print(f"Number of entries after concatenating filtered species dataframes: {len(filtered_ozON)}")

        # Prepare lists for logging
        drop_indices = []
        removed_lipids = []

        # Get sample name from Sample column
        sample_name = filtered_ozON['Sample'].unique()[0]
        log_filename = f"{sample_name}_fame_filter_log_debug.txt"
        log_file_path = os.path.join(log_dir, log_filename)
        csv_filename = f"{sample_name}_fame_filter_removed_lipids_debug.csv"
        csv_file_path = os.path.join(log_dir, csv_filename)

        # Start logging detailed entries
        print(f"Logging details to: {log_file_path}")
        print(f"Logging removed lipids to: {csv_file_path}")

        # Iterate through each entry in fame_std with progress tracking
        for _, fame_entry in tqdm(fame_std.iterrows(), total=len(fame_std), desc=f"Filtering entries for {sample_name}"):
            fame_std_rt = fame_entry['Retention_Time']
            fame_std_species = fame_entry['Species']
            rt_lower_bound = fame_std_rt - fame_rt_window
            rt_upper_bound = fame_std_rt + fame_rt_window
            print(f"\nFame entry - Species: {fame_std_species}, RT: {fame_std_rt}, Window: ({rt_lower_bound}, {rt_upper_bound})")

            # Filter entries within the defined Adjusted_RT_FAME window and matching Species
            for index, row in filtered_ozON.iterrows():
                adjusted_rt = row['Adjusted_RT_FAME']
                species_list = row['Species'].split('|')
                print(f"Row index: {index}, Adjusted_RT_FAME: {adjusted_rt}, Species list: {species_list}")

                if fame_std_species in species_list and rt_lower_bound <= adjusted_rt <= rt_upper_bound:
                    print(f"Match found for Species: {fame_std_species} within RT window.")
                    filtered_ozON.at[index, 'Species'] = fame_std_species
                else:
                    rt_diff = min(abs(adjusted_rt - rt_lower_bound), abs(adjusted_rt - rt_upper_bound))
                    drop_indices.append(index)
                    removed_lipids.append({
                        'Lipid': row['Lipid'],
                        'Adjusted_RT_FAME': adjusted_rt,
                        'Ground_Truth_RT': fame_std_rt,
                        'Species': row['Species'],
                        'Ground_Truth_Species': fame_std_species,
                        'RT_Difference': rt_diff
                    })
                    print(f"Non-matching entry added to removal list - RT difference: {rt_diff}")

        # Drop all collected indices
        filtered_ozON = filtered_ozON.drop(drop_indices).reset_index(drop=True)
        print(f"Number of entries after dropping non-matching rows: {len(filtered_ozON)}")

        # Save filtered DataFrame and removed lipids log to CSV
        filtered_ozON.to_parquet(os.path.join(output_dir, f"{sample_name}_filtered.parquet"))
        pd.DataFrame(removed_lipids).to_csv(csv_file_path, index=False)

        # Write log entries to a text file
        with open(log_file_path, 'w') as log_file:
            for entry in removed_lipids:
                log_file.write(str(entry) + '\n')

    print("filter_OzOFF_by_fame_std function completed.")
