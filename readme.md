# CLAW-MRM - Comprehensive MRM Lipidome Automation Workflow

CLAW-MRM (Comprehensive MRM Lipidome Automation Workflow) is a powerful lipidomics workflow designed to automate and standardize lipid data analysis. It provides a set of tools and scripts that streamline various tasks such as data parsing, matching, statistical analysis, and visualization. This workflow is particularly helpful for researchers in the field of lipidomics as it ensures consistency in data processing and enables efficient exploration and interpretation of lipid expression patterns.

## Getting Started

To get started with CLAW-MRM, follow these steps:

1. **Requirements**: Install the necessary Python libraries by creating a virtual environment. Use the provided requirements or install the following packages:
    ```
    conda install pandas=1.5.2 
    pip install pymzml==2.5.2
    conda install matplotlib
    conda install ipywidgets
    pip install plotly
    conda install ipykernel
    pip install jupyterhub jupyterlab jupyter jsonschema jupyterlab-server
    conda install openpyxl
    pip install -U kaleido
    pip install langchain==0.0.345
    pip install openai==1.3.7
    ```

2. **Project Folder**: Organize your lipidomics project by creating a project folder. This directory will serve as the central location for all project-related files, including raw data, processed results, plots, and other data files.

3. **Workflow**: Utilize the main components of CLAW, such as the Python notebook `Lipid_MRM_parser.ipynb` and the R script `edgeR.R`, to perform data analysis, visualization, and statistical tests on lipidomics datasets. The notebook and script are designed to work together cohesively to achieve comprehensive and reproducible results.

## Repository Structure

```
CLAW/
├── demo/
│   ├── demo_data/
│   │   ├── JSON/
│   │   ├── Labels/
│   │   ├── mzml/
│   │   ├── plots/
│   │   ├── Plots/
│   │   ├── Pre_EdgeR/
│   │   ├── Processed Results/
│   │   └── results/
│   └── jsons/
├── lipid_platform/
│   ├── Figures/
│   ├── lipid_database/
│   ├── Projects/
│   │   ├── BRAIN_5XFAD/
│   │   └── Lipid_Load/
│   ├── tools/
│   └── Variable_Storage/
├── requirements/
└── worklist_generator/
    └── worklists/
```

## Project Folder Structure

When creating a new lipidomics project, use the following folder structure:

```
Project_Folder/
    ├── mzml/                  # Folder for raw data files (mzML format)
    ├── Pre_EdgeR/             # Folder for intermediate processed data
    ├── Plots/                 # Folder for storing generated plots
    ├── Labels/                # Folder for sample labels and metadata
    ├── Results/               # Folder for storing final analysis results
    └── README.md              # Documentation and project overview
```

## Lipid Data Analysis and Visualization

CLAW combines Python and R code for analyzing and visualizing lipid data obtained from mass spectrometry experiments.

### Python Code (Lipid_MRM_parser.ipynb)

The `Lipid_MRM_parser.ipynb` Jupyter Notebook presents a Python-based pipeline for processing and analyzing lipid data:

1. **Loading and filtering data**: Imports data from a lipid database and mzML files, then filters and organizes the data for further analysis.

2. **Parsing and matching**: Uses custom parser to match data with specific lipid classes and extract relevant information.

3. **Data visualization**: Creates various visualizations such as pie charts, bar plots, and edge plots to analyze and compare lipid classes.

4. **Custom color scheme**: Allows specification of a custom color scheme for lipid classes in the visualizations.

### R Code (edgeR.R)

The `edgeR.R` file contains R scripts for lipid data analysis and visualization:

1. **Libraries and Initial Setup**: Loads required R libraries and sets up the working directory.

2. **Loop through Preprocessed Data Files**: Processes preprocessed lipid data files using ridge plots, PCA plots, and heatmap plots.

3. **Ridge Plot Generation**: Creates ridge plots to visualize the distribution of log-fold changes across different lipid classes.

4. **PCA Plot Generation**: Performs PCA on the numeric columns of the data to create PCA plots.

5. **Heatmap Plot Generation**: Generates heatmap plots to visualize the correlation between different lipid features.

6. **Summary and Full Results**: Generates summary and full result files based on edgeR analysis.

## Contributing

If you wish to contribute to CLAW, feel free to submit issues, bug reports, or pull requests on the GitHub repository.

## Acknowledgments

We acknowledge the developers of `pymzml`, `plotly`,`edgeR` and other open-source libraries used in this workflow, which greatly facilitate lipidomics data analysis.

**Happy lipidomics research with CLAW-MRM!**