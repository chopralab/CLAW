# import os
# import pandas as pd
# import matplotlib.pyplot as plt

# class Plot:
#     def __init__(self, dataframe, plot_directory):
#         self.dataframe = dataframe
#         self.base_plot_directory = plot_directory
#         self.ensure_directory_exists()

#     def ensure_directory_exists(self):
#         """Ensure the base directory for plots exists."""
#         if not os.path.exists(self.base_plot_directory):
#             os.makedirs(self.base_plot_directory)
#             print(f"Directory created at {self.base_plot_directory}")
#         else:
#             print(f"Directory already exists at {self.base_plot_directory}")

#     def generate_filename(self, lipid_name):
#         """Generate a safe filename for saving plots based on the lipid name."""
#         safe_lipid_name = lipid_name.replace("/", "-").replace(" ", "_").replace(":", "-")
#         return f"{self.base_plot_directory}/{safe_lipid_name}_OzON.png"

#     def scatter(self):
#         """Create a scatter plot of Peak Area by Lipid."""
#         grouped = self.dataframe.groupby('Lipid')
#         fig, ax = plt.subplots(figsize=(10, 6))
#         for (key, group), color in zip(grouped, plt.cm.tab20.colors):
#             ax.scatter(group['Lipid'], group['Peak_Area'], label=str(key), color=color, s=100)

#         ax.set_xlabel('Lipid')
#         ax.set_ylabel('Peak Area')
#         ax.set_title('Peak Area by Lipid')
#         ax.set_xticks([i for i, _ in enumerate(grouped.groups.keys())])
#         ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
#         ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
#         plt.show()

#     def plot_bar(self):
#         """Create a bar plot of Peak Area by Lipid."""
#         grouped = self.dataframe.groupby('Lipid')
#         fig, ax = plt.subplots(figsize=(10, 6))
#         bar_width = 0.35
#         num_groups = len(grouped.groups.keys())
#         indices = list(range(num_groups))

#         for i, (key, group) in enumerate(grouped):
#             positions = [x + bar_width * i for x in indices]
#             normalized_peak_areas = [group[group['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
#             ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

#         ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
#         ax.set_xticks(ticks_positions)
#         ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
#         ax.set_xlabel('Lipid')
#         ax.set_ylabel('Peak Area')
#         ax.set_title('Peak Area by Lipid')
#         ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
#         plt.show()

#     def plot_bar_by_group_sample(self):
#         """Create bar plots of Peak Area by Lipid for each combination of Biology, Genotype, Mouse, and Cage."""
#         group_cols = ['Biology', 'Genotype', 'Mouse', 'Cage']
#         unique_groups = self.dataframe[group_cols].drop_duplicates()

#         for _, group in unique_groups.iterrows():
#             # Filter the dataframe for the current group combination
#             group_filter = (self.dataframe['Biology'] == group['Biology']) & \
#                            (self.dataframe['Genotype'] == group['Genotype']) & \
#                            (self.dataframe['Mouse'] == group['Mouse']) & \
#                            (self.dataframe['Cage'] == group['Cage'])
#             group_df = self.dataframe[group_filter]

#             # Group by 'Lipid'
#             grouped = group_df.groupby('Lipid')
#             fig, ax = plt.subplots(figsize=(10, 6))
#             bar_width = 0.35
#             num_groups = len(grouped.groups.keys())
#             indices = list(range(num_groups))

#             for i, (key, group_data) in enumerate(grouped):
#                 positions = [x + bar_width * i for x in indices]
#                 # Sum the peak areas for each lipid
#                 normalized_peak_areas = [group_data[group_data['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
#                 ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

#             ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
#             ax.set_xticks(ticks_positions)
#             ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
#             ax.set_xlabel('Lipid')
#             ax.set_ylabel('Peak Area')
#             title = f'Peak Area by Lipid for Biology {group["Biology"]}, Genotype {group["Genotype"]}, Mouse {group["Mouse"]}, Cage {group["Cage"]}'
#             ax.set_title(title)
#             ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
#             plt.show()



### above is old working code take input dataframe

#### new code below for input csv 

### csv input

# OzESI_plot_v2.py

# import os
# import pandas as pd
# import matplotlib.pyplot as plt

# class Plot:
#     def __init__(self, dataframe, plot_directory):
#         self.dataframe = dataframe
#         self.base_plot_directory = plot_directory
#         self.ensure_directory_exists()

#     def ensure_directory_exists(self):
#         """Ensure the base directory for plots exists."""
#         if not os.path.exists(self.base_plot_directory):
#             os.makedirs(self.base_plot_directory)
#             print(f"Directory created at {self.base_plot_directory}")
#         else:
#             print(f"Directory already exists at {self.base_plot_directory}")

#     def generate_filename(self, lipid_name):
#         """Generate a safe filename for saving plots based on the lipid name."""
#         safe_lipid_name = lipid_name.replace("/", "-").replace(" ", "_").replace(":", "-")
#         return f"{self.base_plot_directory}/{safe_lipid_name}_OzON.png"

#     def scatter(self):
#         """Create a scatter plot of Peak Area by Lipid."""
#         grouped = self.dataframe.groupby('Lipid')
#         fig, ax = plt.subplots(figsize=(10, 6))
#         for (key, group), color in zip(grouped, plt.cm.tab20.colors):
#             ax.scatter(group['Lipid'], group['Peak_Area'], label=str(key), color=color, s=100)

#         ax.set_xlabel('Lipid')
#         ax.set_ylabel('Peak Area')
#         ax.set_title('Peak Area by Lipid')
#         ax.set_xticks([i for i, _ in enumerate(grouped.groups.keys())])
#         ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
#         ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
#         plt.show()

#     def plot_bar(self):
#         """Create a bar plot of Peak Area by Lipid."""
#         grouped = self.dataframe.groupby('Lipid')
#         fig, ax = plt.subplots(figsize=(10, 6))
#         bar_width = 0.35
#         num_groups = len(grouped.groups.keys())
#         indices = list(range(num_groups))

#         for i, (key, group) in enumerate(grouped):
#             positions = [x + bar_width * i for x in indices]
#             normalized_peak_areas = [group[group['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
#             ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

#         ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
#         ax.set_xticks(ticks_positions)
#         ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
#         ax.set_xlabel('Lipid')
#         ax.set_ylabel('Peak Area')
#         ax.set_title('Peak Area by Lipid')
#         ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
#         plt.show()

#     def plot_bar_by_group_sample(self):
#         """Create bar plots of Peak Area by Lipid for each combination of Biology, Genotype, Mouse, and Cage."""
#         group_cols = ['Biology', 'Genotype', 'Mouse', 'Cage']
#         unique_groups = self.dataframe[group_cols].drop_duplicates()

#         for _, group in unique_groups.iterrows():
#             # Filter the dataframe for the current group combination
#             group_filter = (self.dataframe['Biology'] == group['Biology']) & \
#                            (self.dataframe['Genotype'] == group['Genotype']) & \
#                            (self.dataframe['Mouse'] == group['Mouse']) & \
#                            (self.dataframe['Cage'] == group['Cage'])
#             group_df = self.dataframe[group_filter]

#             # Group by 'Lipid'
#             grouped = group_df.groupby('Lipid')
#             fig, ax = plt.subplots(figsize=(10, 6))
#             bar_width = 0.35
#             num_groups = len(grouped.groups.keys())
#             indices = list(range(num_groups))

#             for i, (key, group_data) in enumerate(grouped):
#                 positions = [x + bar_width * i for x in indices]
#                 # Sum the peak areas for each lipid
#                 normalized_peak_areas = [group_data[group_data['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
#                 ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

#             ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
#             ax.set_xticks(ticks_positions)
#             ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
#             ax.set_xlabel('Lipid')
#             ax.set_ylabel('Peak Area')
#             title = f'Peak Area by Lipid for Biology {group["Biology"]}, Genotype {group["Genotype"]}, Mouse {group["Mouse"]}, Cage {group["Cage"]}'
#             ax.set_title(title)
#             ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
#             plt.show()


######## updated CSV code below so that the pathing is correct PATH FIXED
import os
import pandas as pd
import matplotlib.pyplot as plt

class Plot:
    def __init__(self, csv_data_folder, peak_analysis_csv_on, peak_analysis_csv_off, ozone_compare_csv, plot_directory):
        self.csv_data_folder = csv_data_folder
        self.peak_analysis_csv_on = peak_analysis_csv_on
        self.peak_analysis_csv_off = peak_analysis_csv_off
        self.ozone_compare_csv = ozone_compare_csv
        self.plot_directory = plot_directory

    def load_data(self, use_ozone_compare=True, mode='OFF'):
        """Load data from the appropriate CSV file based on user input."""
        if use_ozone_compare:
            csv_file_path = os.path.join(self.csv_data_folder, self.ozone_compare_csv)
        else:
            csv_file_path = os.path.join(self.csv_data_folder, self.peak_analysis_csv_on if mode == 'ON' else self.peak_analysis_csv_off)

        print(f"Looking for CSV file at: {csv_file_path}")  # Debugging print statement

        if os.path.exists(csv_file_path):
            df = pd.read_csv(csv_file_path)
            if isinstance(df, pd.DataFrame):
                print(f"Loaded DataFrame with {len(df)} rows and {len(df.columns)} columns.")
                print(f"DataFrame columns: {df.columns}")
                return df, csv_file_path
            else:
                raise TypeError("Loaded data is not a DataFrame.")
        else:
            raise FileNotFoundError(f"File not found: {csv_file_path}")

    def ensure_directory_exists(self):
        """Ensure the base directory for plots exists."""
        if not os.path.exists(self.plot_directory):
            os.makedirs(self.plot_directory)
            print(f"Directory created at {self.plot_directory}")
        else:
            print(f"Directory already exists at {self.plot_directory}")

    # def plot_bar_by_group_sample(self, use_ozone_compare=True, mode='OFF', group_type='Group_Sample'):
    #     """Create bar plots of Peak Area by Lipid for each combination of Biology, Genotype, Mouse, and Cage or Match_Group."""
    #     self.dataframe, csv_file_path = self.load_data(use_ozone_compare, mode)
    #     self.ensure_directory_exists()

    #     if group_type == 'Match_Group':
    #         group_cols = ['Match_Group']
    #     else:
    #         group_cols = ['Biology', 'Genotype', 'Mouse', 'Cage']

    #     unique_groups = self.dataframe[group_cols].drop_duplicates()

    #     for _, group in unique_groups.iterrows():
    #         # Filter the dataframe for the current group combination
    #         if group_type == 'Match_Group':
    #             group_filter = (self.dataframe['Match_Group'] == group['Match_Group'])
    #         else:
    #             group_filter = (self.dataframe['Biology'] == group['Biology']) & \
    #                         (self.dataframe['Genotype'] == group['Genotype']) & \
    #                         (self.dataframe['Mouse'] == group['Mouse']) & \
    #                         (self.dataframe['Cage'] == group['Cage'])
            
    #         group_df = self.dataframe[group_filter]

    #         # Group by 'Lipid'
    #         grouped = group_df.groupby('Lipid')
    #         fig, ax = plt.subplots(figsize=(10, 6))
    #         bar_width = 0.35
    #         num_groups = len(grouped.groups.keys())
    #         indices = list(range(num_groups))

    #         for i, (key, group_data) in enumerate(grouped):
    #             positions = [x + bar_width * i for x in indices]
    #             # Sum the peak areas for each lipid
    #             normalized_peak_areas = [group_data[group_data['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
    #             ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

    #         ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
    #         ax.set_xticks(ticks_positions)
    #         ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
    #         ax.set_xlabel('Lipid')
    #         ax.set_ylabel('Peak Area')
            
    #         if group_type == 'Match_Group':
    #             title_parts = [f'Peak Area by Lipid for Match_Group {group["Match_Group"]} ({os.path.basename(csv_file_path)})']
    #         else:
    #             title_parts = [f'Peak Area by Lipid for Biology {group["Biology"]}, Genotype {group["Genotype"]}, Mouse {group["Mouse"]}, Cage {group["Cage"]} ({os.path.basename(csv_file_path)})']

    #         # Update lipid names if FAC_OFF column exists
    #         if 'FAC_OFF' in self.dataframe.columns:
    #             group_df['Lipid'] = group_df.apply(lambda row: f"{row['Lipid']} + {row['FAC_OFF']}" if not pd.isna(row['FAC_OFF']) else row['Lipid'], axis=1)
    #             title_parts.append('with FAC_OFF')

    #         ax.set_title(' '.join(title_parts))
    #         ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
    #         plt.show()

    def scatter(self, use_ozone_compare=True, mode='OFF'):
        """Create a scatter plot of Intensity by Lipid."""
        self.dataframe, csv_file_path = self.load_data(use_ozone_compare, mode)
        self.ensure_directory_exists()

        # Filter out rows where OzESI_Intensity is zero or missing
        filtered_df = self.dataframe[self.dataframe['OzESI_Intensity'] > 0].dropna(subset=['OzESI_Intensity'])

        if filtered_df.empty:
            print("No lipids with intensities found. Scatter plot will not be displayed.")
            return

        grouped = filtered_df.groupby('Lipid')
        
        # Debug: Print grouped data
        for lipid, group in grouped:
            print(f"Lipid: {lipid}, Data Points: {len(group)}, Intensities: {group['OzESI_Intensity'].tolist()}")

        fig, ax = plt.subplots(figsize=(10, 6))
        colors = plt.cm.tab20.colors * (len(grouped) // len(plt.cm.tab20.colors) + 1)

        for (key, group), color in zip(grouped, colors):
            ax.scatter(group['Lipid'], group['OzESI_Intensity'], label=str(key), color=color, s=100)

        ax.set_xlabel('Lipid')
        ax.set_ylabel('Intensity')
        ax.set_title(f'Intensity by Lipid ({os.path.basename(csv_file_path)})')
        ax.set_xticks([i for i, _ in enumerate(grouped.groups.keys())])
        ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
        ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()


    # def plot_bar(self, use_ozone_compare=True, mode='OFF'):
    #     """Create a bar plot of Intensity by Lipid."""
    #     self.dataframe, csv_file_path = self.load_data(use_ozone_compare, mode)
    #     self.ensure_directory_exists()

    #     # Filter out rows where OzESI_Intensity is zero or missing
    #     filtered_df = self.dataframe[self.dataframe['OzESI_Intensity'] > 0].dropna(subset=['OzESI_Intensity'])

    #     if filtered_df.empty:
    #         print("No lipids with intensities found. Bar plot will not be displayed.")
    #         return

    #     grouped = filtered_df.groupby('Lipid')['OzESI_Intensity'].max().reset_index()
        
    #     fig, ax = plt.subplots(figsize=(10, 6))
    #     ax.bar(grouped['Lipid'], grouped['OzESI_Intensity'], width=0.35)

    #     ax.set_xlabel('Lipid')
    #     ax.set_ylabel('Intensity')
    #     ax.set_title(f'Intensity by Lipid ({os.path.basename(csv_file_path)})')
    #     ax.set_xticks([i for i, _ in enumerate(grouped['Lipid'])])
    #     ax.set_xticklabels(grouped['Lipid'], rotation='vertical')
    #     plt.show()

    def plot_bar(self, use_ozone_compare=True, mode='OFF'):
        """Create a bar plot of Intensity by Lipid."""
        self.dataframe, csv_file_path = self.load_data(use_ozone_compare, mode)
        self.ensure_directory_exists()

        # Filter out rows where OzESI_Intensity is zero or missing
        filtered_df = self.dataframe[self.dataframe['OzESI_Intensity'] > 0].dropna(subset=['OzESI_Intensity'])

        grouped = filtered_df.groupby('Lipid')['OzESI_Intensity'].max().reset_index()
        
        # Debug: Print grouped data
        print(f"Grouped Data: {grouped}")

        fig, ax = plt.subplots(figsize=(10, 6))
        bar_width = 0.35
        colors = plt.cm.tab20.colors * (len(grouped) // len(plt.cm.tab20.colors) + 1)

        for i, (lipid, intensity) in enumerate(zip(grouped['Lipid'], grouped['OzESI_Intensity'])):
            ax.bar(lipid, intensity, width=bar_width, color=colors[i], label=lipid)

        ax.set_xlabel('Lipid')
        ax.set_ylabel('Intensity')
        ax.set_title(f'Intensity by Lipid ({os.path.basename(csv_file_path)})')
        ax.set_xticks(range(len(grouped)))
        ax.set_xticklabels(grouped['Lipid'], rotation='vertical')
        ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()
            
    def plot_bar_by_group_sample(self, use_ozone_compare=True, mode='OFF', group_type='Group_Sample'):
        """Create bar plots of Intensity by Lipid for each combination of Biology, Genotype, Mouse, and Cage or Match_Group."""
        self.dataframe, csv_file_path = self.load_data(use_ozone_compare, mode)
        self.ensure_directory_exists()

        if group_type == 'Match_Group':
            group_cols = ['Match_Group']
        else:
            group_cols = ['Biology', 'Genotype', 'Mouse', 'Cage']

        unique_groups = self.dataframe[group_cols].drop_duplicates()

        for _, group in unique_groups.iterrows():
            # Filter the dataframe for the current group combination
            if group_type == 'Match_Group':
                group_filter = (self.dataframe['Match_Group'] == group['Match_Group'])
            else:
                group_filter = (self.dataframe['Biology'] == group['Biology']) & \
                            (self.dataframe['Genotype'] == group['Genotype']) & \
                            (self.dataframe['Mouse'] == group['Mouse']) & \
                            (self.dataframe['Cage'] == group['Cage'])
            
            group_df = self.dataframe[group_filter]

            # Filter out rows where Intensity is zero or missing
            filtered_df = group_df[group_df['OzESI_Intensity'] > 0].dropna(subset=['OzESI_Intensity'])

            # Group by 'Lipid' and calculate the max intensity
            grouped = filtered_df.groupby('Lipid')['OzESI_Intensity'].max().reset_index()
            
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 0.35
            colors = plt.cm.tab20.colors * (len(grouped) // len(plt.cm.tab20.colors) + 1)

            for i, (lipid, intensity) in enumerate(zip(grouped['Lipid'], grouped['OzESI_Intensity'])):
                ax.bar(lipid, intensity, width=bar_width, color=colors[i], label=lipid)

            ax.set_xlabel('Lipid')
            ax.set_ylabel('Intensity')
            
            if group_type == 'Match_Group':
                title = f'Intensity by Lipid for Match_Group {group["Match_Group"]} ({os.path.basename(csv_file_path)})'
            else:
                title = f'Intensity by Lipid for Biology {group["Biology"]}, Genotype {group["Genotype"]}, Mouse {group["Mouse"]}, Cage {group["Cage"]} ({os.path.basename(csv_file_path)})'

            ax.set_title(title)
            ax.set_xticks(range(len(grouped)))
            ax.set_xticklabels(grouped['Lipid'], rotation='vertical')
            ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.show()


