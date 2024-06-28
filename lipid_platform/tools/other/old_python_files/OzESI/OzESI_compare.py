# import pandas as pd
# import ast

# class OzoneCompare:
#     def __init__(self, directory, file_name_on, file_name_off):
#         self.directory = directory
#         self.file_name_on = file_name_on
#         self.file_name_off = file_name_off
#         self.load_data()

#     def load_data(self):
#         # Load data from Excel files
#         self.OzON_df = pd.read_excel(self.directory + self.file_name_on)
#         self.OzOFF_df = pd.read_excel(self.directory + self.file_name_off)
#         # Drop 'Unnamed: 0' columns if they exist
#         self.OzON_df.drop(columns='Unnamed: 0', errors='ignore', inplace=True)
#         self.OzOFF_df.drop(columns='Unnamed: 0', errors='ignore', inplace=True)

#     def match_dataframes(self):
#         matched = pd.DataFrame(columns=self.OzON_df.columns.tolist() + ['FAC_OFF', 'Retention_Time_OFF'])
#         for _, on_row in self.OzON_df.iterrows():
#             on_fac_list = ast.literal_eval(on_row['FAC'])
#             for _, off_row in self.OzOFF_df.iterrows():
#                 off_fac_list = ast.literal_eval(off_row['FAC'])
#                 for on_fac in on_fac_list:
#                     if on_fac in off_fac_list and abs(on_row['Retention_Time'] - off_row['Retention_Time']) <= 0.5:
#                         new_row = on_row.copy()
#                         new_row['FAC_OFF'] = on_fac
#                         new_row['Retention_Time_OFF'] = off_row['Retention_Time']
#                         matched = matched.append(new_row, ignore_index=True)
#         return matched

#     def print_fac_and_off_match(self, df, match_group):
#         group_data = df[df['Match_Group'] == match_group]
#         if group_data.empty:
#             print(f"No entries found for Match Group {match_group}.")
#         else:
#             print(f"Entries for Match Group {match_group}:")
#             for index, row in group_data.iterrows():
#                 print(f"Index {index}: FAC - {row['FAC']}, FAC_OFF - {row['FAC_OFF']}")

#     def filter_contains_colon_zero(self, df):
#         return df[df['FAC_OFF'].apply(lambda x: ':0' not in x)]
###########################
#########################3

############ USE CSV FILE NOW TO IMPORT


# OzESI_compare.py
# OzESI_compare.py
# OzESI_compare.py
# OzESI_compare.py
# OzESI_compare.py

import pandas as pd
import ast
import os

class OzoneCompare:
    def __init__(self, csv_data_folder, file_name_on, file_name_off, output_file_name):
        self.csv_data_folder = csv_data_folder
        self.file_name_on = file_name_on
        self.file_name_off = file_name_off
        self.output_file_name = output_file_name
        self.load_data()

    def load_data(self):
        # Ensure directory has a trailing slash
        if not self.csv_data_folder.endswith('/'):
            self.csv_data_folder += '/'
        
        # Load data from CSV files
        self.OzON_df = pd.read_csv(os.path.join(self.csv_data_folder, self.file_name_on))
        print(f"Loaded ON data: {self.OzON_df.head()}")  # Debugging print statement
        self.OzOFF_df = pd.read_csv(os.path.join(self.csv_data_folder, self.file_name_off))
        print(f"Loaded OFF data: {self.OzOFF_df.head()}")  # Debugging print statement
        
        # Drop 'Unnamed: 0' columns if they exist
        self.OzON_df.drop(columns='Unnamed: 0', errors='ignore', inplace=True)
        self.OzOFF_df.drop(columns='Unnamed: 0', errors='ignore', inplace=True)

    def match_dataframes(self):
        matched = pd.DataFrame(columns=self.OzON_df.columns.tolist() + ['FAC_OFF', 'Retention_Time_OFF'])
        for _, on_row in self.OzON_df.iterrows():
            try:
                on_fac_list = ast.literal_eval(on_row['FAC'])
            except (ValueError, SyntaxError):
                continue  # Skip rows where 'FAC' cannot be parsed
                
            for _, off_row in self.OzOFF_df.iterrows():
                try:
                    off_fac_list = ast.literal_eval(off_row['FAC'])
                except (ValueError, SyntaxError):
                    continue  # Skip rows where 'FAC' cannot be parsed
                
                for on_fac in on_fac_list:
                    if on_fac in off_fac_list and abs(on_row['Retention_Time'] - off_row['Retention_Time']) <= 0.5:
                        new_row = on_row.copy()
                        new_row['FAC_OFF'] = on_fac
                        new_row['Retention_Time_OFF'] = off_row['Retention_Time']
                        matched = matched.append(new_row, ignore_index=True)
        
        print(f"Matched data: {matched.head()}")  # Debugging print statement
        return matched

    def filter_contains_colon_zero(self, df):
        filtered_df = df[df['FAC_OFF'].apply(lambda x: ':0' not in x)]
        print(f"Filtered data: {filtered_df.head()}")  # Debugging print statement
        return filtered_df
    
    def save_matched_data(self, matched_df):
        output_path = os.path.join(self.csv_data_folder, self.output_file_name)
        matched_df.to_csv(output_path, index=False)
        print(f"Matched data saved to {output_path}")  # Debugging print statement
