import pandas as pd
import re
import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks, peak_widths

def create_all_folders(*folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"Created new folder: {folder}")
        else:
            print(f"Folder already exists: {folder}")



# Function to create a folder if it doesn't exist
def create_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory created at {directory}")
    else:
        print(f"Directory already exists at {directory}")



# Create the directory if it doesn't exist
def create_base_directory(base_plot_directory):
    if not os.path.exists(base_plot_directory):
        os.makedirs(base_plot_directory)
        print(f"Directory created at {base_plot_directory}")
    else:
        print(f"Directory already exists at {base_plot_directory}")

# Define a function to generate filenames based on lipid names
def generate_filename(base_plot_directory, lipid_name):
    safe_lipid_name = lipid_name.replace("/", "-").replace(" ", "_").replace(":", "-")  # Ensure filename is safe for filesystems
    return f"{base_plot_directory}{safe_lipid_name}_OzON.png"

#### for ozone compare

# Function to save a DataFrame to an Excel file
def save_for_ozone_compare(peaks_df, project_results_directory, save_df_name):
    # Ensure the directory exists
    if not os.path.exists(project_results_directory):
        os.makedirs(project_results_directory)
        print(f"Directory created at {project_results_directory}")
    else:
        print(f"Directory already exists at {project_results_directory}")
    
    # Save the DataFrame to an Excel file
    file_path = os.path.join(project_results_directory, save_df_name)
    peaks_df.csv(file_path, index=False)
    print(f'peaks_df saved to excel in results folder {file_path}')

#Import all the necessary libraries

import os








import shutil


def filter_o3mzml_files(project_name):
    """
    Filters .mzML files from the given source subdirectory (located within the project folder)
    and moves them into appropriate subfolders ('o3on' or 'o2only') based on file names.
    """

    # Ensure source subdirectory is within Projects folder and mzml subfolder
    src_directory = os.path.join("Projects", project_name, "mzml")

    # Check if the source directory exists
    if not os.path.exists(src_directory):
        print(f"The directory '{src_directory}' does not exist. Please try again with a valid directory.")
        return

    # Set the destination directories for o3on and o2only (at the same level as the 'mzml' folder)
    dst_directory_o3on = os.path.join("Projects", project_name, "o3on")
    dst_directory_o2only = os.path.join("Projects", project_name, "o2only")

    # If destination directories don't exist, warn the user and exit the function
    if not os.path.exists(dst_directory_o3on) or not os.path.exists(dst_directory_o2only):
        print(f"One or both of the destination directories '{dst_directory_o3on}' and '{dst_directory_o2only}' do not exist.")
        return

    # Loop through the files in the source directory
    for filename in os.listdir(src_directory):
        # Check if the file ends with ".mzML"
        if filename.endswith(".mzML"):
            # Determine which folder to move the file to based on its name
            if "_o3on.mzML" in filename:
                shutil.move(os.path.join(src_directory, filename), os.path.join(dst_directory_o3on, filename))
            elif "_o2only.mzML" in filename:
                shutil.move(os.path.join(src_directory, filename), os.path.join(dst_directory_o2only, filename))






def create_project_folder():
    """
    Creates project folder with subdirectories 'mzml', 'plots', and 'results'.
    """

    while True:
        # Ask user for the project name
        project_name = input("Enter the project name (e.g. 'Mouse_data'): ")

        # Set the main project directory
        project_directory = os.path.join("Projects", project_name)

        # Check if the project directory exists
        if os.path.exists(project_directory):
            user_input = input(f"The directory '{project_directory}' already exists. type exit to exit | or try another name ").strip().lower()
            if user_input != 'exit':
                print("Please enter a different directory.")
                continue
        else:
            print(f"The directory '{project_directory}' does not exist. Creating it now.")

        # Set the subdirectories for mzml, plots, and results
        mzml_directory = os.path.join(project_directory, "mzml")
        plots_directory = os.path.join(project_directory, "plots")
        results_directory = os.path.join(project_directory, "results")

        # Create the subdirectories
        os.makedirs(mzml_directory, exist_ok=True)
        os.makedirs(plots_directory, exist_ok=True)
        os.makedirs(results_directory, exist_ok=True)

        break



