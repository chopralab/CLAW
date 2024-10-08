import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys

class MatchLipids:
    def __init__(self, mrm_database, tolerance=0.3):
        """
        Initialize the Match_Lipids class.

        :param mrm_database: DataFrame containing the MRM database with ion information.
        :param tolerance: Tolerance value for matching ion values.
        """
        self.mrm_database = mrm_database
        self.tolerance = tolerance
        
        # Round the ion values in the MRM database for easier comparison
        self.mrm_database['Parent_Ion'] = np.round(self.mrm_database['Parent_Ion'], 1)
        self.mrm_database['Product_Ion'] = np.round(self.mrm_database['Product_Ion'], 1)

    def within_tolerance(self, values1, values2):
        """
        Checks if two arrays of values are within a specified tolerance.

        :param values1: First array of values to compare.
        :param values2: Second array of values to compare.
        :return: Array of booleans indicating if the values are within tolerance.
        """
        return np.abs(values1 - values2) <= self.tolerance

    def match_lipids_parser(self, df):
        """
        Performs lipid matching by vectorized comparison of ion values.

        :param df: DataFrame containing ion information to be matched.
        :return: DataFrame with matched lipid and class information if matches were found.
        """
        # Round the ion values in the input DataFrame for easier comparison
        df['Parent_Ion'] = np.round(df['Parent_Ion'], 1)
        df['Product_Ion'] = np.round(df['Product_Ion'], 1)
        # Round OzESI_Intensity
        df['OzESI_Intensity'] = np.round(df['OzESI_Intensity'], 0)
        
        matched_lipids = []  # List to store matched lipid names
        matched_classes = []  # List to store matched lipid classes

        # Iterate through each row of the input DataFrame with progress bar
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Matching Lipids"):
            # Check for matches within the specified tolerance
            parent_ion_matches = self.within_tolerance(self.mrm_database['Parent_Ion'].values, row['Parent_Ion'])
            product_ion_matches = self.within_tolerance(self.mrm_database['Product_Ion'].values, row['Product_Ion'])
            matches = parent_ion_matches & product_ion_matches
            
            # If matches are found, append the matched lipid and class information
            if np.any(matches):
                matched_lipids.append(' | '.join(self.mrm_database.loc[matches, 'Lipid'].values))
                matched_classes.append(' | '.join(self.mrm_database.loc[matches, 'Class'].values))
            else:
                matched_lipids.append('')
                matched_classes.append('')
        
        # Add the matched lipid and class information to the DataFrame
        df['Lipid'] = matched_lipids
        df['Class'] = matched_classes
        
        return df

    @staticmethod
    def run_match_lipids(OzOFF_database, OzON_database, df_sample_2, tolerance=0.3, retention_time_window=0.5, output=None, match_on_median=False, start_index=0, end_index=None):
        """
        Function to match lipids between OzOFF and OzON databases and filter df_sample_2 data based on retention time.
        
        Parameters:
        OzOFF_database (pd.DataFrame): DataFrame containing OzOFF database information.
        OzON_database (pd.DataFrame): DataFrame containing OzON database information.
        df_sample_2 (pd.DataFrame): DataFrame containing df_sample_2 data.
        tolerance (float): Tolerance value for matching lipids.
        retention_time_window (float): Retention time window for filtering df_sample_2 data.
        output (str): Path to save the matched results Parquet file.
        match_on_median (bool): If True, match on Retention_Time_Median; otherwise, match on Retention_Time.
        start_index (int): Start index for chunk processing.
        end_index (int): End index for chunk processing.
        
        Returns:
        pd.DataFrame: DataFrame containing the matched and filtered results.
        """
        # Initialize an empty DataFrame to store the master results
        OzON_results_df = pd.DataFrame()

        # Use isin() to filter OzON_database where Lipid matches any Species in OzOFF_database
        species_list = OzOFF_database['Species'].tolist()
        temp_OzON_database = OzON_database[OzON_database['Species'].isin(species_list)]

        # Initialize the Match_Lipids class
        matcher = MatchLipids(temp_OzON_database, tolerance=tolerance)

        # Select the appropriate retention time column
        retention_time_column = 'Retention_Time_Median' if match_on_median else 'Retention_Time'

        # Adjust the end index if it is None or out of bounds
        end_index = end_index or len(OzOFF_database)
        end_index = min(end_index, len(OzOFF_database))

        # Iterate through rows in OzOFF_database with tqdm progress bar
        for index, row in tqdm(OzOFF_database.iloc[start_index:end_index].iterrows(), total=end_index-start_index, desc="Processing OzOFF Database"):
            species = row['Species']
            retention_time_value = row[retention_time_column]
            
            # Further filter the df_sample_2 by Retention_Time within ±retention_time_window
            filtered_df_sample_2 = df_sample_2[
                (df_sample_2['Retention_Time'] >= (retention_time_value - retention_time_window)) &
                (df_sample_2['Retention_Time'] <= (retention_time_value + retention_time_window))
            ].copy()  # Use .copy() to avoid SettingWithCopyWarning
            
            # Add a new column 'Species' to store the respective information
            filtered_df_sample_2.loc[:, 'Species'] = species
            
            # Initialize a temporary DataFrame for the current iteration
            temp_OzON_data = filtered_df_sample_2.copy()

            # Match lipids for the current filtered data
            matched_temp_OzON_data = matcher.match_lipids_parser(temp_OzON_data)
            
            # Append the results to the master DataFrame
            OzON_results_df = pd.concat([OzON_results_df, matched_temp_OzON_data], ignore_index=True)

        # Save the results to a Parquet file if an output path is provided
        if output:
            OzON_results_df.to_parquet(output, index=False)

        # Return the results DataFrame
        return OzON_results_df

# Example usage:
if __name__ == "__main__":
    # Ensure the correct number of arguments are provided
    if len(sys.argv) < 6:
        print("Usage: python match_3.py <OzOFF_database.parquet> <OzON_database.parquet> <df_sample_2.parquet> <start_index> <end_index>")
        sys.exit(1)

    # Parse command-line arguments
    OzOFF_database_path = sys.argv[1]
    OzON_database_path = sys.argv[2]
    df_sample_2_path = sys.argv[3]
    start_index = int(sys.argv[4])
    end_index = int(sys.argv[5])

    # User-defined variables
    tolerance = 0.3  # Set your tolerance value here
    retention_time_window = 0.5  # Set your retention time window here
    matched_data = f"df_match_3_chunk_{start_index}_{end_index}.parquet"  # Save MatchedData output Parquet file
    match_on_median = True  # Set to True if you want to match on Retention_Time_Median

    # Load the input DataFrames
    OzOFF_database = pd.read_parquet(OzOFF_database_path)
    OzON_database = pd.read_parquet(OzON_database_path)
    df_sample_2 = pd.read_parquet(df_sample_2_path)

    # Run the function with user-defined variables, including match_on_median
    OzON_results = MatchLipids.run_match_lipids(
        OzOFF_database, 
        OzON_database, 
        df_sample_2, 
        tolerance=tolerance, 
        retention_time_window=retention_time_window, 
        output=matched_data, 
        match_on_median=match_on_median,
        start_index=start_index,
        end_index=end_index
    )

    # Print results
    print("Filtered OzON Database Matches:")
    print(OzON_results.head())
