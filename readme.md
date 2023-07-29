# CLAW - Comprehensive Lipidome Automation Workflow

CLAW (Comprehensive Lipidome Automation Workflow) is a powerful lipidomics workflow designed to automate and standardize lipid data analysis. It provides a set of tools and scripts that streamline various tasks such as data parsing, matching, statistical analysis, and visualization. This workflow is particularly helpful for researchers in the field of lipidomics as it ensures consistency in data processing and enables efficient exploration and interpretation of lipid expression patterns.

## Getting Started

To get started with CLAW, follow these steps:

1. **Requirements**: Install the necessary Python libraries for `pymzml` by creating a virtual environment. Use the provided `pymzml.yml` file and run the following command:
    ```
    conda env create -f pymzml.yml
    ```

2. **Project Folder**: Organize your lipidomics project by creating a project folder. This directory will serve as the central location for all project-related files, including raw data, processed results, plots, and other data files.

3. **Workflow**: Utilize the main components of CLAW, such as the Python notebook `Lipid_MRM_parser.ipynb` and the R script `edgeR.R`, to perform data analysis, visualization, and statistical tests on lipidomics datasets. The notebook and script are designed to work together cohesively to achieve comprehensive and reproducible results.

## Folder Structure

The project folder should have the following structure:

```
Project_Folder/
    ├── mzml/                  # Folder for raw data files (mzML format)
    ├── Pre_EdgeR/             # Folder for intermediate processed data
    ├── Plots/                 # Folder for storing generated plots
    ├── Labels/                # Folder for sample labels and metadata
    ├── Results/               # Folder for storing final analysis results
    ├── pymzml.yml             # YAML file for conda environment and dependencies
    ├── Lipid_MRM_parser_v2.ipynb   # Python notebook for data parsing and analysis
    ├── main_code.R            # R script for edgeR analysis and visualization
    ├── Figures/               # Folder for supporting figures (if any)
    ├── requirements/          # Folder containing other necessary files
    └── README.md              # Documentation and project overview
```



# Lipid Data Analysis and Visualization

This repository contains R and Python code for analyzing and visualizing lipid data obtained from mass spectrometry experiments. The code is designed to process raw lipid data, perform statistical analysis using edgeR, and generate informative visualizations to gain insights into lipid expression patterns.

## Python Code (Lipid_MRM_parser_v2.ipynb)

The `Lipid_MRM_parser.ipynb` Jupyter Notebook presents a Python-based pipeline for processing and analyzing lipid data. The notebook imports various libraries for data processing and visualization. The key steps in the notebook are as follows:

1. Loading and filtering data: The code imports data from a lipid database and mzML files, and then filters and organizes the data for further analysis.

2. Parsing and matching: The code uses a custom parser (`NO_AVERAGE_SCRIPTS`) to match data with specific lipid classes and extract relevant information.

3. Data visualization: The code creates various visualizations such as pie charts, bar plots, and edge plots to analyze and compare lipid classes.

4. Custom color scheme: The user can specify a custom color scheme for lipid classes in the visualizations.

## R Code (main_code.R)

The `edgeR.R` file contains R scripts for lipid data analysis and visualization. The code is organized into the following sections:

1. Libraries and Initial Setup: The code loads required R libraries and sets up the working directory.

2. Loop through Preprocessed Data Files: The code processes preprocessed lipid data files using ridge plots, PCA plots, and heatmap plots.

3. Ridge Plot Generation: Ridge plots are created to visualize the distribution of log-fold changes across different lipid classes.

4. PCA Plot Generation: PCA is performed on the numeric columns of the data to create PCA plots.

5. Heatmap Plot Generation: Heatmap plots are generated to visualize the correlation between different lipid features.

6. Summary and Full Results: The code generates summary and full result files based on edgeR analysis.

## Conclusion

This repository provides a comprehensive set of tools and scripts for analyzing and visualizing lipid data. Researchers can use the Python notebook for data preprocessing and visualization, while the R scripts offer in-depth analysis using edgeR and additional visualization options. By combining these tools, users can gain valuable insights into lipid expression patterns and make meaningful interpretations from mass spectrometry data.


## Contributing

If you wish to contribute to CLAW, feel free to submit issues, bug reports, or pull requests on the GitHub repository.

## Acknowledgments

We acknowledge the developers of `pymzml`, `plotly`,`edgeR` and other open-source libraries used in this workflow, which greatly facilitate lipidomics data analysis.

**Happy lipidomics research with CLAW!**