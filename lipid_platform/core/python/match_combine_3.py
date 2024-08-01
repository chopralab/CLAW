import os
import pandas as pd
import glob
import gc
import logging
from tqdm import tqdm
from datetime import datetime

def setup_logger():
    """
    Sets up a logger to log messages with timestamps to a file named 'concatenate_parquet.log'.
    
    Returns:
    logger (Logger): Configured logger instance.
    """
    logger = logging.getLogger('concatenate_parquet')
    logger.setLevel(logging.DEBUG)
    
    # Create file handler which logs even debug messages
    fh = logging.FileHandler('concatenate_parquet.log')
    fh.setLevel(logging.DEBUG)
    
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    
    # Add the handlers to the logger
    logger.addHandler(fh)
    
    return logger

def concatenate_parquet_chunks_in_batches(output_directory, temp_directory, output_filename, batch_size=5):
    """
    Concatenates Parquet chunk files in batches and saves intermediate results to avoid memory overflow.
    
    Parameters:
    output_directory (str): Directory where the Parquet chunk files are located.
    temp_directory (str): Directory to save intermediate Parquet files.
    output_filename (str): Filename for the output combined Parquet file.
    batch_size (int): Number of files to process in each batch.
    """
    logger = setup_logger()
    logger.info('Starting the concatenation process.')

    # Get a list of all Parquet files in the directory with 'match' and 'chunk' in the file name
    parquet_files = glob.glob(os.path.join(output_directory, '*match*chunk*.parquet'))

    # Ensure temp directory exists
    os.makedirs(temp_directory, exist_ok=True)

    intermediate_files = []
    batch_count = 0

    # Process files in batches
    for i in tqdm(range(0, len(parquet_files), batch_size), desc="Processing batches"):
        batch_files = parquet_files[i:i + batch_size]
        data_frames = []

        for parquet_file in batch_files:
            logger.info(f"Loading {parquet_file}")
            
            # Load each Parquet file into a DataFrame
            df = pd.read_parquet(parquet_file)
            data_frames.append(df)

        # Concatenate all DataFrames in the current batch
        batch_combined_df = pd.concat(data_frames, ignore_index=True)

        # Save the intermediate batch result to a Parquet file
        intermediate_filename = os.path.join(temp_directory, f'intermediate_{batch_count}.parquet')
        batch_combined_df.to_parquet(intermediate_filename, index=False)
        logger.info(f"Saved intermediate batch {batch_count} to {intermediate_filename}")

        intermediate_files.append(intermediate_filename)
        batch_count += 1

        # Explicitly free memory
        del batch_combined_df, data_frames
        gc.collect()

    # Load all intermediate files and concatenate them into a single DataFrame
    final_data_frames = []
    for file in tqdm(intermediate_files, desc="Loading intermediate files"):
        final_data_frames.append(pd.read_parquet(file))

    final_combined_df = pd.concat(final_data_frames, ignore_index=True)

    # Save the final combined DataFrame to a Parquet file
    final_combined_df.to_parquet(output_filename, index=False)
    logger.info(f"Combined DataFrame saved to {output_filename}")
    logger.info('Finished the concatenation process.')

if __name__ == "__main__":
    # Define the directory containing the Parquet chunk files
    output_directory = './'

    # Define the directory to save intermediate Parquet files
    temp_directory = './temp/'

    # Define the output filename for the combined Parquet file
    output_filename = 'df_match_3_combined.parquet'

    # Run the concatenation function
    concatenate_parquet_chunks_in_batches(output_directory, temp_directory, output_filename)
