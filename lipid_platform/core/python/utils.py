import pandas as pd
import os
import glob

def list_files_in_dirs(dir1, dir2, extension='*.parquet'):
    """
    Lists file names from two directories and creates keys for matching.

    Parameters:
    dir1 (str): Path to the first directory.
    dir2 (str): Path to the second directory.
    extension (str): File extension to match (default: '*.parquet').

    Returns:
    DataFrame: DataFrame containing file names from both directories and their keys.
    """
    # Get list of files in both directories
    files1 = glob.glob(os.path.join(dir1, extension))
    files2 = glob.glob(os.path.join(dir2, extension))

    # Create separate DataFrames for each directory
    df1 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files1],
        'Key': ["_".join(os.path.basename(f).split('_')[2:7]).replace('.parquet', '') for f in files1]
    })
    df2 = pd.DataFrame({
        'File': [os.path.basename(f) for f in files2],
        'Key': ["_".join(os.path.basename(f).split('_')[2:7]).replace('.parquet', '') for f in files2]
    })

    # Remove '5_' prefix from keys if present
    df1['Key'] = df1['Key'].str.replace('^5_', '', regex=True)
    df2['Key'] = df2['Key'].str.replace('^5_', '', regex=True)

    # Print all keys before merging
    print("\nOzOFF Keys:")
    for idx, row in df1.iterrows():
        print(f"{idx + 1}. File: {row['File']}")
        print(f"   Key: {row['Key']}\n")

    print("\nOzON Keys:")
    for idx, row in df2.iterrows():
        print(f"{idx + 1}. File: {row['File']}")
        print(f"   Key: {row['Key']}\n")

    # Merge the DataFrames on the Key column
    merged_df = pd.merge(
        df1, 
        df2, 
        on='Key', 
        how='outer',
        suffixes=('_OzOFF', '_OzON')
    )

    return merged_df
