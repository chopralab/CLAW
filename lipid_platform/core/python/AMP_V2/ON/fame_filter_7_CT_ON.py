import pandas as pd
import os


def FAME_RT_Adjustment(ozon_df, fame_df):
    """
    Adjusts retention times in the ozON DataFrame based on the difference between
    fame STD_RT_OFF (from the row with Sample_ID='FAME_off') and ozON STD_RT_ON values.
    Creates two new columns in the ozON DataFrame:
    - STD_RT_ONvsFAME_Dif: The difference between fame STD_RT_OFF and ozON STD_RT_ON.
    - STD_RT_ONvsFAME_Adjusted: The adjusted retention time for each row in the ozON DataFrame.

    Parameters:
    ozon_df (pd.DataFrame): The ozON dataset containing retention time and STD_RT_ON values.
    fame_df (pd.DataFrame): The fame dataset containing STD_RT_OFF values.

    Returns:
    pd.DataFrame: The updated ozON DataFrame with new columns for differences and adjusted retention times.
    """

    # Ensure Sample_ID columns are cleaned
    fame_df['Sample_ID'] = fame_df['Sample_ID'].astype(str).str.strip()
    ozon_df['Sample_ID'] = ozon_df['Sample_ID'].astype(str).str.strip()

    # Get the fame STD_RT_OFF value where Sample_ID == 'FAME_off'
    fame_off_row = fame_df[fame_df['Sample_ID'] == 'FAME_off']
    if fame_off_row.empty:
        raise ValueError("No row with Sample_ID == 'FAME_off' found in fame DataFrame.")

    # Extract the STD_RT_OFF value
    fame_std_rt_off = fame_off_row['STD_RT_OFF'].values[0]

    # Create a mapping of Sample_ID in ozON to STD_RT_ONvsFAME_Dif
    std_diff_map = {}
    for sample_id in ozon_df['Sample_ID'].unique():
        subset = ozon_df[ozon_df['Sample_ID'] == sample_id]
        
        # Ensure there is exactly one unique STD_RT_ON value for this Sample_ID
        std_rt_on_values = subset['STD_RT_ON'].dropna().unique()
        if len(std_rt_on_values) != 1:
            raise ValueError(f"Sample_ID '{sample_id}' must have exactly one unique STD_RT_ON value.")

        # Calculate the difference between fame_std_rt_off and STD_RT_ON
        std_rt_on = std_rt_on_values[0]
        std_diff_map[sample_id] = fame_std_rt_off - std_rt_on

    # Add the STD_RT_ONvsFAME_Dif column to the ozON DataFrame
    ozon_df['STD_RT_ONvsFAME_Dif'] = ozon_df['Sample_ID'].map(std_diff_map)

    # Calculate the adjusted retention time for each row
    ozon_df['STD_RT_ONvsFAME_Adjusted'] = ozon_df['Retention_Time'] + ozon_df['STD_RT_ONvsFAME_Dif']

    return ozon_df





def filter_ozon_by_fame_std(input_dir, fame_std, output_dir, fame_rt_window=0.5):
    """
    Filters parquet files in the input directory based on fame_std retention times
    and saves the filtered result as parquet files in the output directory.
    Logs and removed lipids CSVs are saved in a log directory inside the output directory.

    Parameters:
    input_dir (str): Directory where input parquet files are located.
    fame_std (pd.DataFrame): DataFrame containing the ground truth species and retention times.
    output_dir (str): Directory where filtered parquet files will be saved.
    fame_rt_window (float): Tolerance for filtering retention times (default: 0.5).

    Returns:
    None: Saves filtered parquet files and logs in the specified output directories.
    """

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Create a log subdirectory inside the output directory
    log_dir = os.path.join(output_dir, 'log')
    os.makedirs(log_dir, exist_ok=True)

    # Iterate over all parquet files in the input directory
    for parquet_file in os.listdir(input_dir):
        if parquet_file.endswith(".parquet"):
            # Read the parquet file
            file_path = os.path.join(input_dir, parquet_file)
            ozon_df = pd.read_parquet(file_path)

            # Adjust retention times using FAME_RT_Adjustment
            ozon_df = FAME_RT_Adjustment(ozon_df, fame_std)

            # Create a new DataFrame to store the filtered ozON data
            filtered_ozON_result = ozon_df.copy()

            # List to store the indices of rows to be dropped and the details of removed lipids
            drop_indices = []
            removed_lipids = []

            # Prepare log entries for writing to the log file
            log_entries = []

            # Get the initial number of lipids (rows) in the filtered_ozON DataFrame
            initial_count = len(filtered_ozON_result)

            # Get the sample name from the Sample column
            sample_name = filtered_ozON_result['Sample'].unique()[0]
            log_filename = f"{sample_name}_fame_filter_log_7.txt"
            log_file_path = os.path.join(log_dir, log_filename)

            # CSV filename for removed lipids
            csv_filename = f"{sample_name}_fame_filter_removed_lipids_7.csv"
            csv_file_path = os.path.join(log_dir, csv_filename)

            # Iterate through each unique Species in fame_std
            for species in fame_std['Species'].unique():
                log_entries.append(f"Processing Species: {species}")
                
                # Get the fame_std entries for this Species
                fame_std_species = fame_std[fame_std['Species'] == species]
                
                # Get the filtered_ozON entries for this Species
                ozON_species = filtered_ozON_result[filtered_ozON_result['Species'] == species]
                
                # Iterate through each isomer in fame_std for this Species
                for isomer in fame_std_species['Isomer'].unique():
                    log_entries.append(f"  Checking Isomer: {isomer}")
                    
                    # Get fame_std Retention_Time for this isomer
                    fame_std_isomer = fame_std_species[fame_std_species['Isomer'] == isomer]
                    fame_std_rt = fame_std_isomer['Retention_Time'].values[0]
                    
                    # Define the retention time window for this isomer
                    lower_bound = fame_std_rt - fame_rt_window
                    upper_bound = fame_std_rt + fame_rt_window
                    log_entries.append(f"    Retention Time Window: {lower_bound:.4f} to {upper_bound:.4f}")
                    
                    # Filter ozON_species based on STD_RT_ONvsFAME_Adjusted within the window
                    for index, row in ozON_species.iterrows():
                        adjusted_rt = row['STD_RT_ONvsFAME_Adjusted']
                        if lower_bound <= adjusted_rt <= upper_bound:
                            log_entries.append(f"    Match found for {row['Lipid']} with STD_RT_ONvsFAME_Adjusted {adjusted_rt:.4f} within the window.")
                        else:
                            difference = min(abs(adjusted_rt - lower_bound), abs(adjusted_rt - upper_bound))
                            log_entries.append(f"    {row['Lipid']} with STD_RT_ONvsFAME_Adjusted {adjusted_rt:.4f} is outside the window by {difference:.4f}. Marking for removal.")
                            # Add the index and the lipid details to the list of rows to drop later and removed lipids
                            drop_indices.append(index)
                            removed_lipids.append({
                                'Lipid': row['Lipid'],
                                'Adjusted_RT': adjusted_rt,
                                'Ground_Truth_RT': fame_std_rt,
                                'Difference': difference
                            })

            # Drop all collected indices at once
            filtered_ozON_result = filtered_ozON_result.drop(drop_indices).reset_index(drop=True)

            # Get the final number of lipids (rows) in the filtered_ozON_result DataFrame
            final_count = len(filtered_ozON_result)
            
            # Add summary to log entries
            log_entries.append(f"\nStarted with {initial_count} Lipids")
            log_entries.append(f"Ended with {final_count} Lipids\n")
            
            # Add removed lipids details to log entries
            if removed_lipids:
                log_entries.append("List of removed lipids:")
                for lipid in removed_lipids:
                    log_entries.append(f"{lipid['Lipid']} with STD_RT_ONvsFAME_Adjusted {lipid['Adjusted_RT']:.4f} is outside the window by {lipid['Difference']:.4f}. Ground Truth Retention_Time: {lipid['Ground_Truth_RT']:.4f}")
            else:
                log_entries.append("No lipids were removed.")

            # Save log entries to text file
            with open(log_file_path, 'w') as log_file:
                for entry in log_entries:
                    log_file.write(entry + '\n')

            # If there are removed lipids, save them to a CSV file
            if removed_lipids:
                removed_lipids_df = pd.DataFrame(removed_lipids)
                removed_lipids_df.columns = ['Lipid', 'Adjusted_RT', 'Ground_Truth_Retention_Time', 'Outside_the_window_by']
                removed_lipids_df.to_csv(csv_file_path, index=False)

            print(f"Log saved to {log_file_path}")
            if removed_lipids:
                print(f"Removed lipids saved to {csv_file_path}")

            # Save the final filtered DataFrame as a parquet file
            parquet_output_path = os.path.join(output_dir, f"{sample_name}_fame_filter_7.parquet")
            filtered_ozON_result.to_parquet(parquet_output_path, index=False)
            print(f"Filtered data saved to {parquet_output_path}")
