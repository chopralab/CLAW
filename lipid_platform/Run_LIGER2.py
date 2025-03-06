import os
import glob
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, norm
import math
from collections import defaultdict



file_path = "/scratch/gilbreth/cbeveri/CLAW/lipid_platform/Variable_Storage/folder_path.txt"

# Read the file and store the value in path_variable
with open(file_path, "r") as file:
    path_variable = file.read().strip()

# ---------------------------
# Configuration
# ---------------------------
# Input directory containing all normalized CSV files
input_dir = path_variable+"csv_clean_normalized"

# Output directory for results (will be created if it doesn't exist)
output_dir = path_variable+"FA_output"
os.makedirs(output_dir, exist_ok=True)

# Path to gene database (remains constant)
gene_database_path = 'tools/gene_database.csv'  # Replace with actual path
gene_database = pd.read_csv(gene_database_path)

# List all CSV files in the input folder
csv_files = glob.glob(os.path.join(input_dir, "*.csv"))

# ---------------------------
# Helper Function for Lipid Parsing
# ---------------------------
def parse_lipid_v2(lipid):
    """Extracts the lipid headgroup and fatty acid composition.
       Returns a tuple (headgroup, chain_length, saturation) if possible."""
    if '(' in lipid and ')' in lipid:
        parts = lipid.split("(")[1].split(")")[0].split(":")
        if parts[0].isdigit():
            return (lipid.split("(")[0], int(parts[0]), int(parts[1]))
        else:
            return lipid.split("(")[1].split(")")[0]
    else:
        return lipid

# ---------------------------
# Analysis for Each Input CSV File
# ---------------------------
for csv_file in csv_files:
    print("Processing file:", csv_file)
    
    # Load the current BioPAN data file
    bio_pan_data = pd.read_csv(csv_file)
    
    # Drop last column if not needed
    bio_pan_data = bio_pan_data.drop(bio_pan_data.columns[-1], axis=1)
    
    # Get the number of columns for group1 intensities (assumed stored in column "Length1")
    length1 = bio_pan_data["Length1"].iloc[0]
    
    # ---------------------------
    # 1. Individual FA Conversion Analysis
    # ---------------------------
    results_revised = []
    for index, row in gene_database.iterrows():
        # Parse lipid info from gene_database
        lipid1_info = parse_lipid_v2(row['Lipid 1'])
        lipid2_info = parse_lipid_v2(row['Lipid 2'])
        
        # Identify matching BioPAN lipids by comparing parsed info
        matching_lipids1 = bio_pan_data['lipid'].apply(lambda x: parse_lipid_v2(x) == lipid1_info)
        matching_lipids2 = bio_pan_data['lipid'].apply(lambda x: parse_lipid_v2(x) == lipid2_info)
        
        if any(matching_lipids1) and any(matching_lipids2):
            lipid1_data = bio_pan_data[matching_lipids1]
            lipid2_data = bio_pan_data[matching_lipids2]
            
            # Define columns for the two groups (intensity data starts at column 11)
            start_col = 11
            split_index = start_col + int(length1)
            
            group1_columns = lipid1_data.columns[start_col:split_index]
            group2_columns = lipid1_data.columns[split_index:]
            
            # Calculate ratios for each group (Lipid2 / Lipid1)
            group1_ratios = lipid2_data[group1_columns].values / lipid1_data[group1_columns].values
            group2_ratios = lipid2_data[group2_columns].values / lipid1_data[group2_columns].values
            
            # Flatten arrays and replace NaNs with zeros (or use nan_policy in t-test)
            ratios_group1 = np.nan_to_num(group1_ratios.flatten())
            ratios_group2 = np.nan_to_num(group2_ratios.flatten())
            
            # Perform Welch's t-test and compute z-score (two-tailed test)
            t_stat, p_value = ttest_ind(ratios_group1, ratios_group2, equal_var=False, nan_policy='omit')
            z_score = norm.ppf(1 - p_value/2)
            
            results_revised.append({
                'Reaction': row['Reaction'],
                'Gene Code': row['Gene Code'] if 'Gene Code' in row else None,
                'Gene': row['Gene'],
                'Lipid 1': row['Lipid 1'],
                'Lipid 2': row['Lipid 2'],
                'T-Statistic': t_stat,
                'P-Value': p_value,
                'Z-Score': z_score
            })
    
    results_df_revised = pd.DataFrame(results_revised)
    
    # Create an output file name based on the input file's base name
    base_name = os.path.splitext(os.path.basename(csv_file))[0]
    individual_output_filename = os.path.join(output_dir, f"FA_individual_{base_name}.csv")
    results_df_revised.to_csv(individual_output_filename, index=False)
    print("Individual analysis completed for:", csv_file)
    
    # ---------------------------
    # 2. Combo Reactions Analysis (from Individual Results)
    # ---------------------------
    def find_combo_reactions(graph):
        """
        Perform DFS on the reaction graph to chain reactions together.
        The combined Z-score is computed as the sum of Z-scores divided by sqrt(n),
        and gene and reaction info is concatenated along the path.
        """
        combos = []
        def dfs(node, path, z_sum, count, visited, gene_path, reaction_path):
            for edge in graph.get(node, []):
                next_node, z, gene_edge, reaction_edge = edge
                if next_node in visited:
                    continue
                new_path = path + [next_node]
                new_z_sum = z_sum + z
                new_count = count + 1
                new_gene_path = gene_path + [gene_edge]
                new_reaction_path = reaction_path + [reaction_edge]
                combos.append({
                    'Start': path[0],
                    'End': next_node,
                    'Path': " -> ".join(new_path),
                    'Edge_Genes': " ; ".join(new_gene_path),
                    'Edge_Reactions': " ; ".join(new_reaction_path),
                    'Num_Reactions': new_count,
                    'Combo_Zscore': new_z_sum / math.sqrt(new_count)
                })
                dfs(next_node, new_path, new_z_sum, new_count, visited | {next_node}, new_gene_path, new_reaction_path)
        for start in graph.keys():
            dfs(start, [start], 0, 0, {start}, [], [])
        return combos
    
    indiv_graph = defaultdict(list)
    for idx, row in results_df_revised.iterrows():
        start = row['Lipid 1']
        end = row['Lipid 2']
        z = row['Z-Score']
        gene = str(row['Gene'])
        reaction = str(row['Reaction'])
        if z is not None:
            indiv_graph[start].append((end, z, gene, reaction))
    
    combo_indiv = find_combo_reactions(indiv_graph)
    combo_indiv_df = pd.DataFrame(combo_indiv)
    
    combo_output_filename = os.path.join(output_dir, f"FA_combo_{base_name}.csv")
    combo_indiv_df.to_csv(combo_output_filename, index=False)
    print("Combo reactions analysis completed for:", csv_file)

