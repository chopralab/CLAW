#Import all the necessary libraries
import pymzml
import csv
import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
import re
import plotly.express as px
from collections import defaultdict
import plotly.io as pio
import json
import plotly.graph_objs as go
import matplotlib.colors as mcolors
import json
import ipywidgets as widgets
import warnings
import time
import shutil



def create_project_folder():
    """
    Filters .mzML files from the given source subdirectory (located within Projects folder)
    and moves them into appropriate subfolders ('o3on' or 'o2only') based on file names.
    """

    while True:
        # Ask user for the source subdirectory (e.g. "04_29_23/mzml")
        src_subdirectory = input("Enter the the project_name (e.g. 'Mouse_data/mzml'): ")

        # Ensure source subdirectory is within Projects folder
        src_directory = os.path.join("Projects", src_subdirectory)

        # Check if the source directory exists
        if os.path.exists(src_directory):
            user_input = input(f"The directory '{src_directory}' already exists. type exit to exit | or try another name ").strip().lower()
            if user_input != 'exit':
                print("Please enter a different directory.")
                continue
        else:
            print(f"The directory '{src_directory}' does not exist. I am creating it. Meow.")
            os.makedirs(src_directory)

        # Set the destination directories for o3on and o2only
        dst_directory_o3on = os.path.join(src_directory, "o3on")
        dst_directory_o2only = os.path.join(src_directory, "o2only")

        # Create the destination directories if they don't exist
        os.makedirs(dst_directory_o3on, exist_ok=True)
        os.makedirs(dst_directory_o2only, exist_ok=True)

        # Loop through the files in the source directory
        for filename in os.listdir(src_directory):
            # Check if the file ends with ".mzML"
            if filename.endswith(".mzML"):
                # Determine which folder to move the file to based on its name
                if "_o3on.mzML" in filename:
                    shutil.move(os.path.join(src_directory, filename), os.path.join(dst_directory_o3on, filename))
                elif "_o2only.mzML" in filename:
                    shutil.move(os.path.join(src_directory, filename), os.path.join(dst_directory_o2only, filename))
        break
