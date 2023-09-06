import pandas as pd
import os

# def process_csv_file(filepath):
    # Load the data from the CSV file into a DataFrame
import os
import pandas as pd

def process_csv(file_path):
    # Load the data
    df = pd.read_csv(file_path)
    print(file_path)

    # Fetch titles and lengths
    title1 = df['Title1'][0].replace(":", "_").replace("|", "_").replace(" ", "_").replace("__","_")
    title2 = df['Title2'][0].replace(":", "_").replace("|", "_").replace(" ", "_")
    length1 = int(df['length1'][0])
    length2 = int(df['length2'][0])

    # Drop columns
    drop_cols = ['Class', 'Title1', 'Title2', 'length1', 'length2',"Title","Blank_name"]
    drop_cols.append(df['Blank_name'][0])  # Drop the Blank value column
    df = df.drop(columns=drop_cols)

    print("")

    # Prepare new column names using title1 and title2
    # value_columns = [col for col in df.columns if col != 'Lipid']
    new_col_names_title1 = [title1 for _ in range(length1)]
    new_col_names_title2 = [title2 for _ in range(length2)]
    new_col_names = new_col_names_title1 + new_col_names_title2
    
    print(list(df))
    # exit()
    # Rename the columns and save
    df.columns = ['Lipid'] + new_col_names
    output_name = f"Biopan_{title1}_vs_{title2}.csv"
    df.to_csv("biopan/"+output_name, index=False)
    print(f"Processed and saved as {output_name}")

# Iterate through CSV files in the "Pre_EdgeR" folder
directory_path = "/home/cbeveri/lipid_parser2/Lipidomics/lipid_platform/Projects/NEW_5xFAD_LIVER_BRAIN_LD/Pre_EdgeR"
for file_name in os.listdir(directory_path):
    if file_name.endswith(".csv"):
        full_file_path = os.path.join(directory_path, file_name)
        process_csv(full_file_path)
#