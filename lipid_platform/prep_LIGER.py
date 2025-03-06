import re
import pandas as pd
import os

# Function to split 'lipid' column by '|' and create new rows for each split
def split_lipid_column(df):
    df_exploded = df.assign(lipid=df['lipid'].str.split('|')).explode('lipid')
    df_exploded = df.assign(lipid=df['lipid'].str.split('|')).explode('lipid')
    return df_exploded

# Function to process Ceramide (Cer) and Ceramide-1-phosphate (CerP) lipids
def process_cer_lipids_v4(df):
    def process_lipid_entry(lipid):
        """
        Process lipid names to:
        - Change 'CerP' to 'Cer1P'
        - Remove 'd' within parentheses
        - Remove content of any nested parentheses
        - Sum numbers in format (XX:Y/ZZ:W) → (XX+ZZ:Y+W)
        """
        # Change 'CerP' to 'Cer1P'
        lipid = lipid.replace('CerP', 'Cer1P')

        # Find the indexes of the first '(' and the last ')'
        start_index = lipid.find('(')
        end_index = lipid.rfind(')')

        if start_index != -1 and end_index != -1:
            content = lipid[start_index+1:end_index]  # Extract content within parentheses

            # Remove 'd' within the main parentheses
            content = content.replace('d', '')

            # Remove content of any nested parentheses
            content = re.sub(r'\([^()]*\)', '', content)

            # Sum numbers in format (XX:Y/ZZ:W) → (XX+ZZ:Y+W)
            if '/' in content:
                parts = content.split('/')
                try:
                    nums = [list(map(int, re.findall(r'\d+', part))) for part in parts]  # <-- FIXED LINE
                    summed_nums = [str(sum(pair)) for pair in zip(*nums)]
                    content = ':'.join(summed_nums)
                except ValueError:
                    return lipid  # Return original lipid if conversion fails

            return lipid[:start_index+1] + content + lipid[end_index:]

        return lipid

    df['lipid'] = df.apply(lambda row: process_lipid_entry(row['lipid']) if row['type'] == 'Cer' else row['lipid'], axis=1)
    return df

# Function to process SM lipids
def process_sm_lipids(df):
    def process_lipid_entry(lipid):
        start_index = lipid.find('(')
        end_index = lipid.rfind(')')

        if start_index != -1 and end_index != -1:
            content = lipid[start_index+1:end_index]
            content = content.replace('d', '')  # Remove 'd'
            content = re.sub(r'\([^()]*\)', '', content)  # Remove nested parentheses

            # Sum numbers in format (XX:Y/ZZ:W) → (XX+ZZ:Y+W)
            if '/' in content:
                try:
                    parts = content.split('/')
                    nums = [list(map(int, re.findall(r'\d+', part))) for part in parts]  # <-- FIXED LINE
                    summed_nums = [str(sum(pair)) for pair in zip(*nums)]
                    content = ':'.join(summed_nums)
                except ValueError:
                    return lipid  # Return original lipid if conversion fails

            return lipid[:start_index+1] + content + lipid[end_index:]

        return lipid

    df['lipid'] = df.apply(lambda row: process_lipid_entry(row['lipid']) if row['type'] == 'SM' else row['lipid'], axis=1)
    return df

# Functions for DAG and TAG lipids
def process_lipid_DAG(lipid):
    lipid = re.sub(r"DG\(O-", "O-DG(", lipid)
    lipid = re.sub(r"_C.*$", "", lipid)
    return lipid

def clean_lipids_DAG(df):
    df['lipid'] = df.apply(lambda row: process_lipid_DAG(row['lipid']) if row['type'] == 'DAG' else row['lipid'], axis=1)
    return df

def process_lipid_TAG(lipid):
    lipid = re.sub(r"TG\(O-", "O-TG(", lipid)
    lipid = re.sub(r"_FA.*$", "", lipid)
    lipid = lipid.replace("[", "").replace("]", "")
    return lipid

def clean_lipids_TAG(df):
    df['lipid'] = df.apply(lambda row: process_lipid_TAG(row['lipid']) if row['type'] == 'TAG' else row['lipid'], axis=1)
    return df

# Functions for PC, PE, PG, PI, PS lipids
def process_lipid_PC(lipid):
    lipid = re.sub(r"PC\(P-", "P-PC(", lipid)
    lipid = re.sub(r"PC\(O-", "O-PC(", lipid)
    return lipid

def process_lipid_PE(lipid):
    lipid = re.sub(r"PE\(P-", "P-PE(", lipid)
    lipid = re.sub(r"PE\(O-", "O-PE(", lipid)
    return lipid

def process_lipid_PI(lipid):
    lipid = re.sub(r"PI\(P-", "P-PI(", lipid)
    lipid = re.sub(r"PI\(O-", "O-PI(", lipid)
    return lipid

def process_lipid_PS(lipid):
    lipid = re.sub(r"PS\(P-", "P-PS(", lipid)
    lipid = re.sub(r"PS\(O-", "O-PS(", lipid)
    return lipid

def process_lipid_PG(lipid):
    lipid = re.sub(r"PG\(P-", "P-PG(", lipid)
    lipid = re.sub(r"PG\(O-", "O-PG(", lipid)
    return lipid

def clean_lipids_PC_PE_PG_PI_PS(df):
    df['lipid'] = df.apply(lambda row: process_lipid_PC(row['lipid']) if row['type'] == 'PC' else row['lipid'], axis=1)
    df['lipid'] = df.apply(lambda row: process_lipid_PE(row['lipid']) if row['type'] == 'PE' else row['lipid'], axis=1)
    df['lipid'] = df.apply(lambda row: process_lipid_PS(row['lipid']) if row['type'] == 'PS' else row['lipid'], axis=1)
    df['lipid'] = df.apply(lambda row: process_lipid_PI(row['lipid']) if row['type'] == 'PI' else row['lipid'], axis=1)
    df['lipid'] = df.apply(lambda row: process_lipid_PG(row['lipid']) if row['type'] == 'PG' else row['lipid'], axis=1)
    return df



file_path = "/scratch/gilbreth/cbeveri/CLAW/lipid_platform/Variable_Storage/folder_path.txt"

# Read the file and store the value in path_variable
with open(file_path, "r") as file:
    path_variable = file.read().strip()

# Processing all CSV files in a folder
# input_folder = path_variable+"results"
# output_folder = path_variable+"csv_clean_normalized"
    
import os

input_folder = os.path.join(path_variable, "results")
output_folder = os.path.join(path_variable, "csv_clean_normalized")

os.makedirs(input_folder, exist_ok=True)
os.makedirs(output_folder, exist_ok=True)


for filename in os.listdir(input_folder):
    if filename.endswith("full.csv"):
        df = pd.read_csv(os.path.join(input_folder, filename))
        # df = df.iloc[:, :-1]
        df = split_lipid_column(df)
        df = process_cer_lipids_v4(df)
        df = process_sm_lipids(df)
        df = clean_lipids_DAG(df)
        df = clean_lipids_TAG(df)
        df = clean_lipids_PC_PE_PG_PI_PS(df)

        # Normalize all numerical columns (ignoring first 11 columns)
        # for col in df.columns[11:]:
        #     column_sum = df[col].sum()
        #     if column_sum != 0:
        #         df[col] = df[col] / column_sum

        # Save cleaned data
        df.to_csv(os.path.join(output_folder, filename), index=False)
