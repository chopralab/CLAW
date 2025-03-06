import os
import glob
import re
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, norm
from math import sqrt
from collections import defaultdict
import re
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, norm
from math import sqrt
from collections import defaultdict


file_path = "/scratch/gilbreth/cbeveri/CLAW/lipid_platform/Variable_Storage/folder_path.txt"

# Read the file and store the value in path_variable
with open(file_path, "r") as file:
    path_variable = file.read().strip()



output_dirs = {
    "matched_pairs": os.path.join(path_variable, "matched_pairs"),
    "matched_ttest": os.path.join(path_variable, "matched_ttest"),
    "matched_ttest_group": os.path.join(path_variable, "matched_ttest_group"),
    "sum": os.path.join(path_variable, "sum"),
    "sum_combo": os.path.join(path_variable, "sum_combo")
}


for dir_path in output_dirs.values():
    os.makedirs(dir_path, exist_ok=True)

matched_pairs_string = os.path.join(path_variable, "matched_pairs/")
matched_ttest_string = os.path.join(path_variable, "matched_ttest/")
matched_ttest_group_string = os.path.join(path_variable, "matched_ttest_group/")
sum_dir_string = os.path.join(path_variable, "sum/")
sum_combo_string = os.path.join(path_variable, "sum_combo/")

csv_folder = os.path.join(path_variable, 'csv_clean_normalized')
csv_files = glob.glob(os.path.join(csv_folder, '*.csv'))
print("CSV files found:", csv_files)


# ---------------------------
# Load Gene Database (unchanged)
# ---------------------------
gene_database_path = 'tools/gene_database.csv'      # Update with your actual path
gene_database = pd.read_csv(gene_database_path)

# ---------------------------
# Loop over all CSV files in csv_clean_normalized (except gene_database.csv)
# ---------------------------
for file_path in glob.glob(os.path.join(csv_folder, '*.csv')):

    if os.path.basename(file_path) == "gene_database.csv":
        continue

    filename = os.path.splitext(os.path.basename(file_path))[0]
    print(filename)


    







    
    bio_pan_data = pd.read_csv(file_path)

    # Drop the last column if not needed
    bio_pan_data = bio_pan_data.drop(bio_pan_data.columns[-1], axis=1)
    length1 = bio_pan_data["Length1"].iloc[0]  # number of intensity columns for group 1

    print(bio_pan_data["lipid"].iloc[0])
    print(type(bio_pan_data["lipid"].iloc[0]))

    # ---------------------------
    # 2. Define FA Exclusion, FA_ratios, and Helper Functions
    # ---------------------------
    fa_pattern = re.compile(r"^FA\(\d+:\d+\)$")  # Exclude lipids that match exactly FA(XX:Y)




    FA_ratios = ['24:6', '24:5', '16:0', '24:1', '12:0', '12:1', '14:0', '14:1', '15:0', '15:1', '16:1', '17:0', '17:1', '18:0',
                    '18:1', '18:2', '18:3', '18:4', '19:0', '20:0', '20:1', '20:3', '20:4', '20:5', '22:0', '22:1', '22:4', 
                    '22:5', '22:6', '24:0', '26:0', '26:1', '28:0', '30:0', '32:0', '34:0']




    def find_matches_for_lipid(sample_data, lipid):
        lipid_prefix = lipid.split("(")[0]
        return [s for s in sample_data if s.startswith(lipid_prefix)]

    def parse_lipid_value(lipid_str):
        parts = lipid_str.split(":")
        return int(parts[0]), int(parts[1])

    def is_difference_in_FA_ratios(content1, content2):
        len1_val, sat1 = parse_lipid_value(content1)
        len2_val, sat2 = parse_lipid_value(content2)
        diff = f"{abs(len1_val - len2_val)}:{abs(sat1 - sat2)}"
        return diff in FA_ratios

    def headgroup(lipid):
        """Return the headgroup (prefix) of a lipid string (e.g., 'PC' from 'PC(32:0)')."""
        if '(' in lipid:
            return lipid.split('(')[0]
        return lipid

    # ---------------------------
    # 3. Generate Matched Pairs (from v1 and v2 methods)
    # ---------------------------
    matched_pairs_set = set()
    for index, row in gene_database.iterrows():
        lipid1 = row['Lipid 1']
        lipid2 = row['Lipid 2']
        number = row['Number']
        gene = row['Gene']
        reaction = row['Reaction']
        if fa_pattern.match(lipid1) or fa_pattern.match(lipid2):
            continue
        match1 = find_matches_for_lipid(bio_pan_data['lipid'], lipid1)
        match2 = find_matches_for_lipid(bio_pan_data['lipid'], lipid2)
        # Method A: Equality matching (v1 style)
        for m1 in match1:
            for m2 in match2:
                if fa_pattern.match(m1) or fa_pattern.match(m2):
                    continue
                if '(' in m1 and '(' in m2:
                    content_m1 = m1.split("(")[1].split(")")[0]
                    content_m2 = m2.split("(")[1].split(")")[0]
                    if content_m1 == content_m2:
                        matched_pairs_set.add((m1, m2, lipid1, lipid2, number, gene, reaction))
        # Method B: Difference matching (v2 style)
        for m1 in match1:
            for m2 in match2:
                if fa_pattern.match(m1) or fa_pattern.match(m2):
                    continue
                if '(' in m1 and '(' in m2:
                    content_m1 = m1.split("(")[1].split(")")[0]
                    content_m2 = m2.split("(")[1].split(")")[0]
                    if is_difference_in_FA_ratios(content_m1, content_m2):
                        matched_pairs_set.add((m1, m2, lipid1, lipid2, number, gene, reaction))
    matched_pairs = list(matched_pairs_set)
    columns_pairs = ['Match 1', 'Match 2', 'Lipid 1', 'Lipid 2', 'Number', 'Gene', 'Reaction']
    matched_pairs_df = pd.DataFrame(matched_pairs, columns=columns_pairs)
    
    save_path_1 = matched_pairs_string+filename+"matched_pairs.csv"
    matched_pairs_df.to_csv(save_path_1, index=False)

    # ---------------------------
    # 4. Calculate T-test Results for Individual Reactions (v1)
    # ---------------------------
    df = bio_pan_data
    def calculate_ratios(row_match1, row_match2):
        ratios = []
        start_col = 11  # assume first 11 columns are metadata
        for col in range(start_col, len(row_match1.columns)):
            val = row_match1.iloc[0, col]
            if val != 0:
                ratios.append(row_match2.iloc[0, col] / val)
            else:
                ratios.append(None)
        return ratios

    def calculate_z_score(p_value):
        try:
            return norm.ppf(1 - p_value / 2)
        except Exception:
            return np.nan

    all_results = []
    for m1, m2, l1, l2, n, g, r in matched_pairs:
        row_match1 = df[df['lipid'] == m1]
        row_match2 = df[df['lipid'] == m2]
        if not row_match1.empty and not row_match2.empty:
            pair_ratios = calculate_ratios(row_match1, row_match2)
            group1 = [v for v in pair_ratios[:int(length1)] if v is not None]
            group2 = [v for v in pair_ratios[int(length1):] if v is not None]
            if group1 and group2:
                t_stat, p_val = ttest_ind(group1, group2, equal_var=False, nan_policy='omit')
                z_val = calculate_z_score(p_val)
            else:
                t_stat, p_val, z_val = np.nan, np.nan, np.nan
            all_results.append({
                'Match 1': m1,
                'Match 2': m2,
                'Lipid 1': l1,
                'Lipid 2': l2,
                'Number': n,
                'Gene': g,
                'Reaction': r,
                'T-Statistic': t_stat,
                'P-Value': p_val,
                'Z-Score': z_val
            })
    results_df = pd.DataFrame(all_results)
    save_path_2 = matched_ttest_string+filename+"matched_pairs_ttest.csv"

    results_df.to_csv(save_path_2, index=False)

    # --- Aggregate duplicate reaction edges
    group_cols = ['Reaction', 'Lipid 1', 'Lipid 2', 'T-Statistic', 'P-Value', 'Z-Score']
    results_grouped = results_df.groupby(group_cols, as_index=False).agg({
        'Gene': lambda x: ','.join(sorted(set(x))),
        'Match 1': 'first',
        'Match 2': 'first',
        'Number': 'first'
    })
    save_path_3 = matched_ttest_group_string+filename+"matched_pairs_ttest_grouped.csv"

    results_grouped.to_csv(save_path_3, index=False)
    aggregated_df = results_grouped  # Use this aggregated data for the next steps

    # ---------------------------
    # 5. Sum Composition Results (v2)
    # ---------------------------
    start_col = 11
    sum_comp_results = []
    for key, group in aggregated_df.groupby(['Gene', 'Reaction', 'Lipid 1', 'Lipid 2']):
        substrate_sum = None
        product_sum = None
        for idx, row in group.iterrows():
            m1 = row['Match 1']
            m2 = row['Match 2']
            row1 = df[df['lipid'] == m1]
            row2 = df[df['lipid'] == m2]
            if row1.empty or row2.empty:
                continue
            intensities1 = row1.iloc[0, start_col:].astype(float)
            intensities2 = row2.iloc[0, start_col:].astype(float)
            if substrate_sum is None:
                substrate_sum = intensities1.copy()
                product_sum = intensities2.copy()
            else:
                substrate_sum += intensities1
                product_sum += intensities2
        if substrate_sum is None or product_sum is None:
            continue
        ratios = []
        for sub, prod in zip(substrate_sum, product_sum):
            if sub != 0:
                ratios.append(prod / sub)
            else:
                ratios.append(None)
        group1 = [v for v in ratios[:int(length1)] if v is not None]
        group2 = [v for v in ratios[int(length1):] if v is not None]
        if group1 and group2:
            t_stat, p_val = ttest_ind(group1, group2, equal_var=False, nan_policy='omit')
            z_val = calculate_z_score(p_val)
        else:
            t_stat, p_val, z_val = np.nan, np.nan, np.nan
        sum_comp_results.append({
            'Gene': key[0],
            'Reaction': key[1],
            'Lipid 1': key[2],
            'Lipid 2': key[3],
            'T-Statistic': t_stat,
            'P-Value': p_val,
            'Z-Score': z_val
        })
    sum_comp_df = pd.DataFrame(sum_comp_results)
    
    save_path_4 = sum_dir_string+filename+"sum.csv"
    
    sum_comp_df.to_csv(save_path_4, index=False)

    # ---------------------------
    # 6. Build Combo Reactions from the Sum Composition Results
    # ---------------------------
    # New approach: read the previously saved combined_sum_composition_results.csv and use it to build combo reactions.
    sum_comp_df = pd.read_csv(save_path_4)

    def build_graph_from_sum_comp(df_results, use_headgroup=False):
        """
        Build a reaction graph from the sum composition results.
        Each unique reaction edge is represented once.
        If use_headgroup is True, convert lipid names to headgroups.
        Expected columns: 'Gene', 'Reaction', 'Lipid 1', 'Lipid 2', 'T-Statistic', 'P-Value', 'Z-Score'
        """
        graph = defaultdict(list)
        for idx, row in df_results.iterrows():
            if use_headgroup:
                start = headgroup(row['Lipid 1'])
                end = headgroup(row['Lipid 2'])
            else:
                start = row['Lipid 1']
                end = row['Lipid 2']
            z = row['Z-Score']
            gene_val = row['Gene']
            reaction_val = row['Reaction']
            if not pd.isna(z):
                graph[start].append((end, z, gene_val, reaction_val))
        return graph

    # Build the reaction graph from the sum composition results.
    reaction_graph = build_graph_from_sum_comp(sum_comp_df, use_headgroup=True)

    def combine_z_scores(z_list, method="stouffer"):
        """Combine a list of Z-scores using the specified method."""
        if not z_list:
            return np.nan
        if method == "stouffer":
            return sum(z_list) / sqrt(len(z_list))
        elif method == "mean":
            return sum(z_list) / len(z_list)
        else:
            raise ValueError("Unknown combination method: choose 'stouffer' or 'mean'")

    def find_combo_reactions(graph, combine_method="stouffer"):
        """
        Run DFS from all nodes in the reaction graph.
        Each unique path (chain of reactions) is recorded with its combined Z-score.
        The final output includes columns: Start, End, Path, Edge_Genes, Edge_Reactions, Num_Reactions, Combo_Zscore.
        """
        all_paths = []
        def dfs(node, path, z_list, visited, gene_path, reaction_path):
            for edge in graph.get(node, []):
                next_node, z, gene_edge, reaction_edge = edge
                if next_node in visited:
                    continue
                new_path = path + [next_node]
                new_z_list = z_list + [z]
                new_gene_path = gene_path + [gene_edge]
                new_reaction_path = reaction_path + [reaction_edge]
                combo_z = combine_z_scores(new_z_list, method=combine_method)
                all_paths.append({
                    'Start': path[0],
                    'End': next_node,
                    'Path': " -> ".join(new_path),
                    'Edge_Genes': " ; ".join(str(x) for x in new_gene_path),
                    'Edge_Reactions': " ; ".join(str(x) for x in new_reaction_path),
                    'Num_Reactions': len(new_z_list),
                    'Combo_Zscore': combo_z
                })
                dfs(next_node, new_path, new_z_list, visited | {next_node}, new_gene_path, new_reaction_path)
        for start in graph.keys():
            dfs(start, [start], [], {start}, [], [])
        # Deduplicate paths by their "Path" string (keeping the one with the highest Combo_Zscore).
        unique = {}
        for entry in all_paths:
            key = entry['Path']
            if key not in unique or entry['Combo_Zscore'] > unique[key]['Combo_Zscore']:
                unique[key] = entry
        return list(unique.values())

    combo_results = find_combo_reactions(reaction_graph, combine_method="stouffer")
    combo_df = pd.DataFrame(combo_results)
    
    save_path_5 = sum_combo_string+filename+"_combosum.csv"
    
    combo_df.to_csv(save_path_5, index=False)

    print("All aggregated output files have been saved.")








