import os
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

def plot_lipid_intensities(input_csv, output_dir, sorting='Lipid'):
    """
    Create bar plots for lipid intensities for each unique file in the input CSV.
    Each file gets its own directory, and multiple plots are created if there are more than 10 lipids.

    Args:
        input_csv (str): Path to the CSV file containing parsed chromatogram data.
        output_dir (str): Directory to save the plots.
        sorting (str): Sorting criterion, either 'Lipid' or 'Intensity'. Default is 'Lipid'.

    Returns:
        None
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load the CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # Group by Filename to create separate plots for each file
    grouped = df.groupby('Filename')

    for filename, group in tqdm(grouped, desc="Processing files"):
        # Create a directory for the current filename
        file_dir = os.path.join(output_dir, filename)
        os.makedirs(file_dir, exist_ok=True)

        # Sort the group based on the specified criterion
        if sorting == 'Intensity':
            group = group.sort_values(by='Intensity', ascending=False)
        elif sorting == 'Lipid':
            group = group.sort_values(by='Lipid', ascending=True)
        else:
            raise ValueError("Invalid sorting criterion. Use 'Lipid' or 'Intensity'.")

        # Split the lipids into chunks of 10
        lipid_chunks = [group[i:i + 10] for i in range(0, len(group), 10)]

        for i, chunk in enumerate(lipid_chunks):
            plt.figure(figsize=(10, 6))
            plt.bar(chunk['Lipid'], chunk['Intensity'], alpha=0.7, color='blue')
            plt.title(f"Lipid Intensities - {filename} (Chunk {i + 1})")
            plt.xlabel('Lipid')
            plt.ylabel('Intensity')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()

            # Save the plot in the specific directory for the current filename
            plot_file = os.path.join(file_dir, f"{filename}_chunk_{i + 1}_lipid_intensities.png")
            plt.savefig(plot_file)
            plt.close()

            print(f"Plot saved: {plot_file}")

