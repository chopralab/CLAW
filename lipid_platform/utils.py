import pandas as pd

def save_df_csv(df: pd.DataFrame, path: str) -> None:
    """
    Saves a DataFrame to a CSV file at the specified path, appending '.csv' to the file name if necessary.

    Args:
        df (pd.DataFrame): The DataFrame to save.
        path (str): The file path and name where the DataFrame should be saved, without the '.csv' extension.

    Returns:
        None
    """
    # Check if the path ends with '.csv', append it if not
    if not path.endswith('.csv'):
        path += '.csv'
    
    # Save the DataFrame to a CSV file
    df.to_csv(path, index=False)
    print(f"DataFrame saved to {path}")



import pandas as pd
import matplotlib.pyplot as plt

def plot_matched_scatter_and_create_df(df1, df2, match_columns, ozesi_intensity_threshold=50):
    """
    Plots a scatter plot for each match between two dataframes based on specified columns and creates a new DataFrame
    with one matching row from each DataFrame for each group. Removes matches where the difference in Correct_RT
    is greater than 3 and OzESI_Intensity is below 300.

    Args:
    df1 (pandas.DataFrame): The first DataFrame (OzON).
    df2 (pandas.DataFrame): The second DataFrame (OzOFF).
    match_columns (list): List of column names to match between df1 and df2.
    ozesi_intensity_threshold (float): The threshold for OzESI_Intensity. Default is 300.

    Returns:
    pandas.DataFrame: A DataFrame with one matching row from each DataFrame for each group.
    """
    # Apply OzESI_Intensity threshold
    df1 = df1[df1['OzESI_Intensity'] >= ozesi_intensity_threshold]
    df2 = df2[df2['OzESI_Intensity'] >= ozesi_intensity_threshold]

    # Select the first row for each group in both dataframes
    df1_grouped = df1.groupby(match_columns).first().reset_index()
    df2_grouped = df2.groupby(match_columns).first().reset_index()

    # Find common rows based on the match_columns
    common_rows = pd.merge(df1_grouped, df2_grouped, on=match_columns, suffixes=('_OzON', '_OzOFF'))

    # Calculate the difference in Correct_RT and filter out rows with difference > 0.5
    common_rows['RT_Diff'] = (common_rows['Correct_RT_OzON'] - common_rows['Correct_RT_OzOFF']).abs()
    common_rows = common_rows[common_rows['RT_Diff'] <= 0.8]

    # Assign a unique group number for each combination of match_columns
    common_rows['Group_Sample'] = common_rows.groupby(match_columns).ngroup()

    # Set the figure size
    plt.figure(figsize=(24, 6))

    # Plotting
    plt.scatter(common_rows['Group_Sample'], common_rows['Correct_RT_OzON'], label='OzON', alpha=0.7)
    plt.scatter(common_rows['Group_Sample'], common_rows['Correct_RT_OzOFF'], label='OzOFF', alpha=0.7, color='r')
    plt.xlabel('Samples', fontsize=16)
    plt.ylabel('Retention Time', fontsize=16)
    plt.title('Brain5xFAD OzOFF vs OzON Sample Comparison', fontsize=20)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)
    plt.show()

    # Create a new DataFrame with Group_Sample and the match_columns
    df_OzOFFvsOzON_Matching = common_rows[['Group_Sample'] + match_columns + ['Correct_RT_OzON', 'Correct_RT_OzOFF']]

    return df_OzOFFvsOzON_Matching


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  # Ensure numpy is imported for isinstance check

def plot_lipid_groups(d_OzON2):
    """
    Plots OzESI Intensity for each lipid group in a DataFrame over the Retention Time.

    Args:
    d_OzON2 (pandas.DataFrame): DataFrame containing 'Lipid', 'Biology', 'Cage', 'Genotype', 
                                'Retention_Time', and 'OzESI_Intensity' columns.

    """
    # Group by Lipid
    lipid_groups = d_OzON2.groupby("Lipid")

    # Determine the number of unique lipids to create subplots
    num_lipids = len(lipid_groups)
    fig, axes = plt.subplots(num_lipids, 1, figsize=(10, 6 * num_lipids), squeeze=False)

    # Plot each lipid group in a separate subplot
    for (lipid, group), ax in zip(lipid_groups, axes.flatten()):
        for _, sub_group in group.groupby(["Biology", "Cage", "Genotype"]):
            ax.plot(sub_group["Retention_Time"], sub_group["OzESI_Intensity"], marker='o', linestyle='', ms=8)
        ax.set_title(f"Lipid: {lipid}")
        ax.set_xlabel("Retention Time")
        ax.set_ylabel("OzESI Intensity")

    plt.tight_layout()
    plt.show()
