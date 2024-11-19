import re
import os
import pandas as pd

def parse_chromatogram_data(input_dir, output_dir):
    """
    Parse chromatogram data from text files in the specified input directory
    and save the results as a CSV in the specified output directory.

    Args:
        input_dir (str): Directory containing the input text files.
        output_dir (str): Directory to save the output CSV file.

    Returns:
        pd.DataFrame: DataFrame containing the parsed chromatogram data.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'parsed_chromatogram_data.csv')

    # Initialize lists to store parsed data
    filenames = []
    q1_values = []
    q3_values = []
    summed_intensities = []
    lipids = []
    dates = []          # List for Date
    sample_names = []   # List for Sample_Name
    samples = []        # List for Sample

    # Iterate through all .txt files in the specified directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.txt'):
            file_path = os.path.join(input_dir, file_name)

            # Extract Date and Sample_Name from the filename
            base_name = os.path.splitext(file_name)[0]  # Removes the .txt extension
            parts = base_name.split('_', 1)             # Split only on the first underscore
            if len(parts) == 2:
                date_str = parts[0]                     # e.g., '20241115'
                sample_name = parts[1]                  # e.g., 'Plasma_Acyl-Carnitines'
            else:
                # Handle unexpected filename formats
                date_str = ''
                sample_name = ''

            # Determine the Sample value based on Sample_Name
            if 'Blank' in sample_name:
                sample = 'Blank'
            else:
                sample = 'Sample'

            # Open and read the file
            with open(file_path, 'r') as file:
                lines = file.readlines()

                current_filename = ""
                current_q1 = None
                current_q3 = None
                current_intensity_sum = 0
                current_lipid = ""
                parsing_intensity = False

                for line in lines:
                    # Extract the filename (assumes filename appears earlier in the file)
                    if 'sourceFile:' in line or 'name:' in line:
                        match = re.search(r'name:\s+([\w.]+)', line)
                        if match:
                            current_filename = match.group(1)

                    # Find Q1, Q3 values, and Lipid
                    if 'id: SRM SIC Q1=' in line:
                        match = re.search(r'Q1=(\d+\.\d+).*Q3=(\d+\.\d+).*name=([^\s]+)', line)
                        if match:
                            current_q1 = float(match.group(1))
                            current_q3 = float(match.group(2))
                            current_lipid = match.group(3)

                    # Check if we are parsing intensity array data
                    if 'cvParam: intensity array' in line:
                        parsing_intensity = True
                        current_intensity_sum = 0
                    elif parsing_intensity and 'binary: [' in line:
                        # Extract and sum intensity values
                        match = re.search(r'binary:\s+\[\d+\]\s+([\d\s]+)', line)
                        if match:
                            intensities = map(int, match.group(1).split())
                            current_intensity_sum = sum(intensities)
                        parsing_intensity = False

                        # Append the extracted data to the lists
                        if current_filename and current_q1 is not None and current_q3 is not None:
                            filenames.append(current_filename)
                            q1_values.append(current_q1)
                            q3_values.append(current_q3)
                            summed_intensities.append(current_intensity_sum)
                            lipids.append(current_lipid)
                            dates.append(date_str)              # Append Date
                            sample_names.append(sample_name)    # Append Sample_Name
                            samples.append(sample)              # Append Sample

    # Create a DataFrame
    chromatogram_df = pd.DataFrame({
        'Date': dates,                        # Date column
        'Sample_Name': sample_names,          # Sample_Name column
        'Sample': samples,                    # Sample column
        'Lipid': lipids,
        'Q1': q1_values,
        'Q3': q3_values,
        'Intensity': summed_intensities,
        'Filename': filenames,
    })

    # Create Base_Sample_Name by removing 'Blank_' prefix if present
    chromatogram_df['Base_Sample_Name'] = chromatogram_df['Sample_Name'].str.replace('Blank_', '', regex=False)

    # Assign group numbers based on Base_Sample_Name
    chromatogram_df['Blank_Group'] = pd.factorize(chromatogram_df['Base_Sample_Name'])[0] + 1  # Start groups at 1

    # Optionally, drop the Base_Sample_Name column if not needed
    chromatogram_df.drop(columns=['Base_Sample_Name'], inplace=True)

    # Optionally, convert Date to datetime format
    # chromatogram_df['Date'] = pd.to_datetime(chromatogram_df['Date'], format='%Y%m%d')

    # Reorder columns for better readability
    columns_order = ['Date', 'Sample_Name', 'Sample', 'Blank_Group', 'Lipid', 'Q1', 'Q3', 'Intensity', 'Filename']
    chromatogram_df = chromatogram_df[columns_order]

    # Save the DataFrame to a CSV file
    chromatogram_df.to_csv(output_file, index=False)

    return chromatogram_df
