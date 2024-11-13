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

    # Iterate through all .txt files in the specified directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.txt'):
            file_path = os.path.join(input_dir, file_name)

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

    # Create a DataFrame
    chromatogram_df = pd.DataFrame({
        'Lipid': lipids,
        'Q1': q1_values,
        'Q3': q3_values,
        'Intensity': summed_intensities,
        'Filename': filenames,
    })

    # Save the DataFrame to a CSV file
    chromatogram_df.to_csv(output_file, index=False)

    return chromatogram_df
