# from msilib.schema import Component
import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA

from matplotlib.colors import ListedColormap, LinearSegmentedColormap

threshold = 'FDR'
# threshold = 'PValue'

import re


# Original names
names222 = [
    "DOD100-M2-5xFAD-cerebellum10x_01202023", "DOD100-M2-5xFAD-cortex10x_01202023",
    "DOD100-M2-5xFAD-diencephalon10x_01202023", "DOD100-M2-5xFAD-hippo10x_01202023",
    "DOD100-M3-5xFAD-cerebellum10x_01202023", "DOD100-M3-5xFAD-diencephalon10x_01202023",
    "DOD100-M3-5xFAD-hippo10x_01202023", "DOD99-M1-5xFAD-cerebellum10x_01202023",
    "DOD99-M1-5xFAD-cortex10x_01202023", "DOD99-M1-5xFAD-diencephalon10x_01202023",
    "DOD99-M1-5xFAD-hippo10x_01202023", "M1FAD173-5xFAD-Cerebellum10x_01182023",
    "M1FAD173-5xFAD-Cortex10x_01182023", "M1FAD173-5xFAD-Diencephalon10x_01182023",
    "M1FAD173-5xFAD-Hippo10x_01182023", "M2FAD173-5xFAD-Cerebellum10x_01192023",
    "M2FAD173-5xFAD-Cortex10x_01192023", "M2FAD173-5xFAD-Diencephalon10x_01192023",
    "M2FAD173-5xFAD-Hippo10x_01192023", "DOD100-M1-WT-cerebellum10x_01202023",
    "DOD100-M1-WT-cortex10x_01202023", "DOD100-M1-WT-diencephalon10x_01202023",
    "DOD100-M1-WT-hippo10x_01202023", "DOD99-M2-WT-cerebellum10x_01202023",
    "DOD99-M2-WT-cortex10x_01202023", "DOD99-M2-WT-diencephalon10x_01202023",
    "DOD99-M2-WT-hippo10x_01202023", "M5FAD173-WT-Cerebellum10x_01192023",
    "M5FAD173-WT-Cortex10x_01192023", "M5FAD173-WT-Diencephalon10x_01192023",
    "M5FAD173-WT-Hippo10x_01192023", "10xBlank_01202023"
]

# New hardcoded names
new_names = [
    "5xFAD_Male_2_Cerebellum", "5xFAD_Male_2_Cortex",
    "5xFAD_Male_2_Diencephalon", "5xFAD_Male_2_Hippocampus",
    "5xFAD_Male_3_Cerebellum", "5xFAD_Male_3_Diencephalon",
    "5xFAD_Male_3_Hippocampus", "5xFAD_Male_1_Cerebellum",
    "5xFAD_Male_1_Cortex", "5xFAD_Male_1_Diencephalon",
    "5xFAD_Male_1_Hippocampus", "5xFAD_Male_1_Cerebellum",
    "5xFAD_Male_1_Cortex", "5xFAD_Male_1_Diencephalon",
    "5xFAD_Male_1_Hippocampus", "5xFAD_Male_2_Cerebellum",
    "5xFAD_Male_2_Cortex", "5xFAD_Male_2_Diencephalon",
    "5xFAD_Male_2_Hippocampus", "WT_Male_1_Cerebellum",
    "WT_Male_1_Cortex", "WT_Male_1_Diencephalon",
    "WT_Male_1_Hippocampus", "WT_Male_2_Cerebellum",
    "WT_Male_2_Cortex", "WT_Male_2_Diencephalon",
    "WT_Male_2_Hippocampus", "WT_Male_5_Cerebellum",
    "WT_Male_5_Cortex", "WT_Male_5_Diencephalon",
    "WT_Male_5_Hippocampus", "10xBlank_01202023"
]


# Set up the file path
file_name = 'heatmap_test.csv'
full_file_path = 'heatmap_test/' + file_name

# Load data from CSV or Excel file
try:
    main_dataframe = pd.read_csv(full_file_path)
except:
    main_dataframe = pd.read_excel(full_file_path)

# print(file_name[:-5])

# Display unique types
unique_types = main_dataframe['type'].unique()
# print(unique_types)

# Loop through each unique value in the 'type' column
for type_value in unique_types:
    # Filter the dataframe based on the current type
    filtered_dataframe = main_dataframe[main_dataframe['type'] == type_value]
    
    # Getting column names and their count
    column_names = list(filtered_dataframe)
    num_columns = len(column_names)
    
    # print(column_names)
    # print(num_columns)
    
    last_column_name = column_names[-1]
    # print(last_column_name)

    for index in range(3, num_columns-1):
        print(column_names[index])

    columns_to_process = column_names[3:-1]

    # print(last_column_name)
    # print(column_names[3])
    
    PCA_title = f"PCA_{4}_componets_DIvided_by_blank{file_name[:-4]}"
    heatmap_title = f"heat_maps_by_class/{file_name[:-4]} {type_value}.svg"
    heatmap_title2 = f"heat_maps_by_class/{file_name[:-4]}_Log10_heatmap_division.png"

    blank_values = np.array(filtered_dataframe[last_column_name])

    # print(blank_values)

    intensity_values = []
    normalized_intensity_values = []
    non_normalized_values = []
    sample_names = []
    line_values = []

for column_name in columns_to_process:
    print(type(filtered_dataframe))
    print(column_name)
    
    delta_intensity = np.array(filtered_dataframe[column_name]) - blank_values
    normalized_intensity = np.array(filtered_dataframe[column_name]) / blank_values
    raw_intensity = filtered_dataframe[column_name].tolist()

    sample_names.append(column_name)
    sorted_raw_intensity = sorted(raw_intensity, reverse=True)

    line_values.append(raw_intensity)
    intensity_values.append(delta_intensity)
    normalized_intensity_values.append(normalized_intensity)
    non_normalized_values.append(sorted_raw_intensity)


    # print(len(intensity_values))
    # print(len(sample_names))

    intensity_values = np.array(intensity_values)
    intensity_values[intensity_values < 0] = 0

    print(sample_names)
    sorted_blank_values = sorted(blank_values, reverse=True)
    print(non_normalized_values[0])

    x_values = list(range(len(non_normalized_values[0])))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def save_and_clear_plot(filename, dpi=500):
    plt.savefig(filename, dpi=dpi)
    plt.clf()
    plt.close()

def generate_heatmap(data, colormap, vmin, vmax, y_labels, title, filename_prefix):
    heatmap = plt.pcolor(data, cmap=colormap, vmin=vmin, vmax=vmax)
    plt.yticks(ticks=np.arange(len(y_labels)), labels=y_labels)
    plt.colorbar(heatmap)
    plt.title(title)
    plt.tight_layout()
    plt.show()
    save_and_clear_plot(f"{filename_prefix}.png")
    save_and_clear_plot(f"{filename_prefix}.svg")

def generate_line_plot(x, y_data, y_label, title, filename):
    plt.figure()
    plt.plot(x, y_data, color='red', label=y_label)
    plt.plot(x, blank_values.tolist(), color='blue', label="blank")
    plt.title(title)
    plt.xlabel("Unsorted of Lipid Intensities")
    plt.ylabel("Intensity")
    plt.legend()
    plt.tight_layout()
    save_and_clear_plot(filename)

# Convert blank values to list and sort in descending order
sorted_blank_values = sorted(blank_values, reverse=True)

x_values = np.arange(len(non_normalized_values[0])).tolist()

# Colormap configurations
colormap1 = LinearSegmentedColormap.from_list("mycmap", [[0,"blue"], [0.5,"yellow"], [1,"red"]])
colormap2 = LinearSegmentedColormap.from_list("mycmap", [(0,"black"), (.149999,"black"), (.15,"purple"), 
                                                        (.2,"purple"), (.2001,"blue"), (.3,"blue"), 
                                                        (0.3000000001,"red"), (.5,"yellow"), (1,"green")])

# Categorize lipid types
lipid_categories = ["TAGs" if "TAG" in lipid_type else lipid_type for lipid_type in filtered_dataframe["type"]]

# Create a list with lipid types, replacing consecutive duplicates with spaces
deduplicated_lipid_categories = [lipid_categories[0]] + [
    lipid_categories[i] if lipid_categories[i] != lipid_categories[i-1] else " " 
    for i in range(1, len(lipid_categories))
]

# Generate the first heatmap
generate_heatmap(np.log10(intensity_values), colormap1, 0, 6, sample_names, type_value, heatmap_title)

# Generate the second heatmap
generate_heatmap(normalized_intensity_values, colormap2, 0, 10, sample_names, f"Intensity/Blank {type_value}", heatmap_title2)

# Generate line plots
for idx in range(len(line_values)):
    generate_line_plot(x_values, line_values[idx], sample_names[idx], 
                       f"{sample_names[idx]}{type_value} UnSorted List of Intensities", 
                       f"Ion_Distribution_by_class/{sample_names[idx]}{type_value}_Non_sorted_.png")

print(np.amin(intensity_values))
