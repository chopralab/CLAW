import os
import pandas as pd
import matplotlib.pyplot as plt

class Plot:
    def __init__(self, dataframe, plot_directory):
        self.dataframe = dataframe
        self.base_plot_directory = plot_directory
        self.ensure_directory_exists()

    def ensure_directory_exists(self):
        """Ensure the base directory for plots exists."""
        if not os.path.exists(self.base_plot_directory):
            os.makedirs(self.base_plot_directory)
            print(f"Directory created at {self.base_plot_directory}")
        else:
            print(f"Directory already exists at {self.base_plot_directory}")

    def generate_filename(self, lipid_name):
        """Generate a safe filename for saving plots based on the lipid name."""
        safe_lipid_name = lipid_name.replace("/", "-").replace(" ", "_").replace(":", "-")
        return f"{self.base_plot_directory}/{safe_lipid_name}_OzON.png"

    def scatter(self):
        """Create a scatter plot of Peak Area by Lipid."""
        grouped = self.dataframe.groupby('Lipid')
        fig, ax = plt.subplots(figsize=(10, 6))
        for (key, group), color in zip(grouped, plt.cm.tab20.colors):
            ax.scatter(group['Lipid'], group['Peak_Area'], label=str(key), color=color, s=100)

        ax.set_xlabel('Lipid')
        ax.set_ylabel('Peak Area')
        ax.set_title('Peak Area by Lipid')
        ax.set_xticks([i for i, _ in enumerate(grouped.groups.keys())])
        ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
        ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()

    def plot_bar(self):
        """Create a bar plot of Peak Area by Lipid."""
        grouped = self.dataframe.groupby('Lipid')
        fig, ax = plt.subplots(figsize=(10, 6))
        bar_width = 0.35
        num_groups = len(grouped.groups.keys())
        indices = list(range(num_groups))

        for i, (key, group) in enumerate(grouped):
            positions = [x + bar_width * i for x in indices]
            normalized_peak_areas = [group[group['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
            ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

        ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
        ax.set_xticks(ticks_positions)
        ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
        ax.set_xlabel('Lipid')
        ax.set_ylabel('Peak Area')
        ax.set_title('Peak Area by Lipid')
        ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.show()

    def plot_bar_by_group_sample(self):
        """Create bar plots of Peak Area by Lipid for each combination of Biology, Genotype, Mouse, and Cage."""
        group_cols = ['Biology', 'Genotype', 'Mouse', 'Cage']
        unique_groups = self.dataframe[group_cols].drop_duplicates()

        for _, group in unique_groups.iterrows():
            # Filter the dataframe for the current group combination
            group_filter = (self.dataframe['Biology'] == group['Biology']) & \
                           (self.dataframe['Genotype'] == group['Genotype']) & \
                           (self.dataframe['Mouse'] == group['Mouse']) & \
                           (self.dataframe['Cage'] == group['Cage'])
            group_df = self.dataframe[group_filter]

            # Group by 'Lipid'
            grouped = group_df.groupby('Lipid')
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 0.35
            num_groups = len(grouped.groups.keys())
            indices = list(range(num_groups))

            for i, (key, group_data) in enumerate(grouped):
                positions = [x + bar_width * i for x in indices]
                # Sum the peak areas for each lipid
                normalized_peak_areas = [group_data[group_data['Lipid'] == lipid]['Peak_Area'].sum() for lipid in grouped.groups.keys()]
                ax.bar(positions, normalized_peak_areas, width=bar_width, label=str(key))

            ticks_positions = [i + bar_width * (num_groups / 2) - bar_width / 2 for i in indices]
            ax.set_xticks(ticks_positions)
            ax.set_xticklabels(grouped.groups.keys(), rotation='vertical')
            ax.set_xlabel('Lipid')
            ax.set_ylabel('Peak Area')
            title = f'Peak Area by Lipid for Biology {group["Biology"]}, Genotype {group["Genotype"]}, Mouse {group["Mouse"]}, Cage {group["Cage"]}'
            ax.set_title(title)
            ax.legend(title='Lipid', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.show()
