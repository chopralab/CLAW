import re
import os
import pandas as pd

def TIC_RSD(intensities):
    """
    Calculate the Relative Standard Deviation (RSD) for a given intensity array.

    Args:
        intensities (list of int): The intensity values.

    Returns:
        float: The RSD percentage.
    """
    if not intensities:
        return 0.0
    mean = sum(intensities) / len(intensities)
    if mean == 0:
        return 0.0
    variance = sum((x - mean) ** 2 for x in intensities) / len(intensities)
    std_dev = variance ** 0.5
    rsd = (std_dev / mean) * 100
    return rsd

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
    lipids = []
    dates = []          # List for Date
    sample_names = []   # List for Sample_Name
    samples = []        # List for Sample
    summed_intensities = []  # List for Summed_Intensity
    tic_rsd_values = [] # List for TIC_RSD

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
            sample = 'Blank' if 'Blank' in sample_name else 'Sample'

            # Open and read all lines of the file
            with open(file_path, 'r') as file:
                lines = file.readlines()

            # -------------------
            # Extract TIC Data
            # -------------------
            TIC_intensities = []
            for i, line in enumerate(lines):
                if 'id: TIC' in line:
                    # Look for the intensity array within TIC section
                    for j in range(i, len(lines)):
                        if 'binaryDataArray:' in lines[j] and 'intensity array' in lines[j]:
                            # The intensity values are in the next line
                            if j + 1 < len(lines):
                                match = re.search(r'binary:\s+\[\d+\]\s+([\d\s]+)', lines[j + 1])
                                if match:
                                    TIC_intensities = list(map(int, match.group(1).split()))
                                    break
                    break  # Assuming only one TIC section per file

            # Calculate TIC_RSD using the TIC_RSD function
            TIC_rsd = TIC_RSD(TIC_intensities)

            # -------------------
            # Parse Lipid Data
            # -------------------
            current_filename = ""
            current_q1 = None
            current_q3 = None
            current_lipid = ""
            parsing_intensity = False
            intensities = []

            for line in lines:
                # Extract the filename (assumes filename appears earlier in the file)
                if 'sourceFile:' in line or 'name:' in line:
                    match = re.search(r'name:\s+([\w.]+)', line)
                    if match:
                        current_filename = match.group(1)

                # Extract Q1, Q3 values, and Lipid name
                if 'id: SRM SIC Q1=' in line:
                    match = re.search(r'Q1=(\d+\.\d+).*Q3=(\d+\.\d+).*name=([^\s]+)', line)
                    if match:
                        current_q1 = float(match.group(1))
                        current_q3 = float(match.group(2))
                        current_lipid = match.group(3)

                # Check if we are parsing intensity array data for lipids
                if 'cvParam: intensity array' in line:
                    parsing_intensity = True
                    intensities = []
                elif parsing_intensity and 'binary: [' in line:
                    # Extract intensity values
                    match = re.search(r'binary:\s+\[\d+\]\s+([\d\s]+)', line)
                    if match:
                        intensities = list(map(int, match.group(1).split()))
                        current_intensity_sum = sum(intensities)

                        # Append the extracted and calculated data to the lists
                        if current_filename and current_q1 is not None and current_q3 is not None:
                            filenames.append(current_filename)
                            q1_values.append(current_q1)
                            q3_values.append(current_q3)
                            lipids.append(current_lipid)
                            dates.append(date_str)              # Append Date
                            sample_names.append(sample_name)    # Append Sample_Name
                            samples.append(sample)              # Append Sample
                            summed_intensities.append(current_intensity_sum)  # Append Summed_Intensity
                            tic_rsd_values.append(TIC_rsd)      # Append TIC_RSD

                    parsing_intensity = False

    # Create a DataFrame
    chromatogram_df = pd.DataFrame({
        'Date': dates,                        # Date column
        'Sample_Name': sample_names,          # Sample_Name column
        'Sample': samples,                    # Sample column
        'Lipid': lipids,
        'Q1': q1_values,
        'Q3': q3_values,
        'Summed_Intensity': summed_intensities,  # Summed_Intensity column
        'TIC_RSD': tic_rsd_values,            # TIC_RSD column
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
    columns_order = [
        'Date', 
        'Sample_Name', 
        'Sample', 
        'Blank_Group', 
        'Lipid', 
        'Q1', 
        'Q3', 
        'Summed_Intensity',
        'TIC_RSD',                           # Include TIC_RSD in the order
        'Filename'
    ]
    chromatogram_df = chromatogram_df[columns_order]

    # Save the DataFrame to a CSV file
    chromatogram_df.to_csv(output_file, index=False)

    return chromatogram_df
