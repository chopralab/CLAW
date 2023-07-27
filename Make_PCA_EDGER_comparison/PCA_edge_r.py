

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def make_pca_plot(df, file_path,title):


    ###Load The dataframe

    ###Filter the dataframe

    #####



    # transpose dataframe so that each column becomes a data point

    df_transposed = df.transpose()

    # standardize the data (mean=0, std=1)
    standardized_data = StandardScaler().fit_transform(df_transposed)

    # PCA transformation
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(standardized_data)

    # create a dataframe with the principal components
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])


    ##Need to update this

    # add labels based on column names
    pca_df['label'] = ['C' if name.startswith('C') else 'R' for name in df_transposed.index]

    # get percentage of variance explained by each PC
    explained_var = pca.explained_variance_ratio_

    # plot
    plt.figure(figsize=(10,10))
    for label, color in [('C', 'red'), ('R', 'blue')]:
        mask = pca_df['label'] == label
        plt.scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], c=color, edgecolors='black', label=label)
    
    plt.title(title)
    plt.xlabel(f'PC1 ({explained_var[0]*100:.2f}%)')  # include % variance explained in label
    plt.ylabel(f'PC2 ({explained_var[1]*100:.2f}%)')  # include % variance explained in label
    plt.legend()

    # save plot
    plt.savefig(file_path+title + '.png')






####Full PCA

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.graph_objs as go





lst=os.listdir('../data_files_name_changed')
print("df")

for ij in lst:
    if ".csv" in ij:
        # try:
        print("Hi")
        pi_title = ij[:-4]
            ##Gets path to open dataframe
        full_file_path = '../data_files_name_changed/'+ij


        df = pd.read_csv(full_file_path)

        old_list = list(df)
        print(old_list)
        # exit()
        lipid_names = list(df['lipid'])
        PVALUE = list(df['PValue'])
        logFC = list(df['logFC'])
        FDR = list(df['FDR'])

        lipid_types_old = list(df['type'])
        transition = list(df['Transition'])

        df = df.drop(['lipid', 'Transition','lipid','PValue','logFC','FDR','type',"logCPM","LR"], axis=1)


###PCA Information




def read_data():
    # Read the EAE_part2.xlsx file
    df1 = pd.read_excel('EAE_part2.xlsx')

    # Read the sample_names_for_PCA.csv file
    df2 = pd.read_csv('sample_names_for_PCA.csv')

    # Rename the 'Name' column in df2 to match the value columns in df1
    df2 = df2.rename(columns={'Name': 'value_column'})

    # Convert the wide format of df1 to long format
    df1_long = df1.melt(id_vars=['lipid', 'Transition', 'type', 'Blank'], var_name='value_column', value_name='value')

    # Merge the two DataFrames on the value_column
    merged_df = df1_long.merge(df2, on='value_column')
    
    return merged_df

def blank_subtraction(df):
    df['value'] = df['value'] - df['Blank']
    df.loc[df['value'] < 0, 'value'] = 0
    df = df.drop(columns=['Blank'])  # Drop the 'Blank' column
    return df

def filter_dataframe(df, cell_types=None, diseases=None, sample_types=None):
    if cell_types:
        df = df[df['Cell Type'].isin(cell_types)]
    if diseases:
        df = df[df['Disease'].isin(diseases)]
    if sample_types:
        df = df[df['Sample Type'].isin(sample_types)]
        
    return df




def normalize_data(df):
    data = df.pivot_table(index=['lipid', 'Transition', 'type'], columns='value_column', values='value').reset_index()
    data = data.set_index(['lipid', 'Transition', 'type']).T
    normalized_data = StandardScaler().fit_transform(data)
    return normalized_data, data



def perform_pca(normalized_data, data, df, n_components=2):
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(normalized_data)
    principal_components_df = pd.DataFrame(data=principal_components, columns=[f'Principal Component {i + 1}' for i in range(n_components)], index=data.index)
    
    # Retrieve 'Sample Type', 'Disease', and 'Cell Type' from the original DataFrame
    metadata = df[['value_column', 'Sample Type', 'Disease', 'Cell Type']].drop_duplicates().set_index('value_column')
    
    # Concatenate the principal components and metadata DataFrames
    final_df = pd.concat([metadata.loc[principal_components_df.index], principal_components_df], axis=1)
    # Add PC1 Variance and PC2 Variance columns to the final_df
    final_df['PC1 Variance'] = pca.explained_variance_ratio_[0]
    final_df['PC2 Variance'] = pca.explained_variance_ratio_[1]


    return final_df, pca




###PCA plot with specific and general legend
def plot_pca_with_all_legend(final_df, cell_type_colors, disease_shape_dict, title="sample", file_name="test", marker_size=10):
    # Create a scatter plot of the principal component values
    fig = go.Figure()
    for sample_type, color in cell_type_colors.items():
        data = final_df[final_df['Cell Type'] == sample_type]
        for disease, shape in disease_shape_dict.items():
            subdata = data[data['Disease'] == disease]
            fig.add_trace(go.Scatter(x=subdata['Principal Component 1'], y=subdata['Principal Component 2'],
                                     mode='markers', marker=dict(color=color, symbol=shape, size=marker_size), 
                                     name=f"{sample_type} - {disease}"))

    # Map diseases to shapes in final_df
    final_df['Disease Shape'] = final_df['Disease'].map(disease_shape_dict)

    # Set the marker shapes
    for i, disease in enumerate(disease_shape_dict.keys()):
        fig.update_traces(marker=dict(symbol=disease_shape_dict[disease]), selector=dict(name=f".*{disease}.*"))

    # Set the title and axis labels
    fig.update_layout(title=title, xaxis_title='Principal Component 1 (PC1)', yaxis_title='Principal Component 2 (PC2)')

    # Add the PC variance to the axis titles
    variances = final_df[['PC1 Variance', 'PC2 Variance']].iloc[0].tolist()
    variances = [f"{v:.1%}" for v in variances]
    fig.update_layout(xaxis_title=f"Principal Component 1 (PC1)<br>({variances[0]} variance explained)",
                      yaxis_title=f"Principal Component 2 (PC2)<br>({variances[1]} variance explained)")

    # Add the color and shape legends
    fig.update_traces(showlegend=True)
    fig.update_layout(legend=dict(traceorder='normal', itemsizing='constant', font=dict(size=10)),
                      legend_title_text='', legend_title_font=dict(size=12),
                      colorway=[color for sample_type, color in cell_type_colors.items()])
    for sample_type, color in cell_type_colors.items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color=color),
                                 legendgroup='Color Legend', showlegend=True, name=sample_type))
    for disease, shape in disease_shape_dict.items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(symbol=shape, color='black'),
                                 legendgroup='Shape Legend', showlegend=True, name=disease))

    # Save the plot as an HTML file and an SVG file
    pio.write_html(fig, file_name + ".html")
    pio.write_image(fig, file_name + ".svg")

    # Show the plot in the Jupyter notebook
    fig.show()

# ###Correct PCA separate legends!!!! ##No hovering
def plot_pca(final_df, cell_type_colors, disease_shape_dict, title="sample", file_name="test", marker_size=10):
    # Create a scatter plot of the principal component values
    fig = go.Figure()

    for sample_type, color in cell_type_colors.items():
        data = final_df[final_df['Cell Type'] == sample_type]
        for disease, shape in disease_shape_dict.items():
            subdata = data[data['Disease'] == disease]
            hover_text = subdata['Disease'] + " - " + subdata['Cell Type']
            fig.add_trace(go.Scatter(x=subdata['Principal Component 1'], y=subdata['Principal Component 2'],
                                     mode='markers', marker=dict(color=color, symbol=shape, size=marker_size),
                                     showlegend=False, hovertemplate="%{text}<extra></extra>", text=hover_text))

    # Set the title and axis labels
    fig.update_layout(title=title, xaxis_title='Principal Component 1 (PC1)', yaxis_title='Principal Component 2 (PC2)')

    # Add the PC variance to the axis titles
    variances = final_df[['PC1 Variance', 'PC2 Variance']].iloc[0].tolist()
    variances = [f"{v:.1%}" for v in variances]
    fig.update_layout(xaxis_title=f"Principal Component 1 (PC1)<br>({variances[0]} variance explained)",
                      yaxis_title=f"Principal Component 2 (PC2)<br>({variances[1]} variance explained)")

    # Add the color and shape legends
    color_legend_items = []
    shape_legend_items = []
    for sample_type, color in cell_type_colors.items():
        color_legend_items.append((sample_type, color))
    for disease, shape in disease_shape_dict.items():
        shape_legend_items.append((disease, shape))
    fig.update_layout(legend=dict(traceorder='normal', itemsizing='constant', font=dict(size=10)),
                      legend_title_text='', legend_title_font=dict(size=12),
                      colorway=[color for sample_type, color in cell_type_colors.items()])
    for sample_type, color in cell_type_colors.items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color=color),
                                 legendgroup='Color Legend', showlegend=True, name=sample_type))
    for disease, shape in disease_shape_dict.items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(symbol=shape, color='black'),
                                 legendgroup='Shape Legend', showlegend=True, name=disease))
    # Add dummy traces for the color and shape legends
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color='black', symbol='circle'),
                             showlegend=True, name='Color Legend'))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color='black', symbol='square'),
                             showlegend=True, name='Shape Legend'))

    # Save the plot as an HTML file and an SVG file
    pio.write_html(fig, file_name + ".html")
    pio.write_image(fig, file_name + ".svg")

    # Show the plot in the Jupyter notebook
    fig.show()



merged_df = read_data()

merged_df.to_csv("Test_merged.csv")
# Perform blank subtraction
merged_df = blank_subtraction(merged_df)



# Get all unique cell types and diseases
unique_cell_types = list(merged_df['Cell Type'].unique())
unique_diseases = list(merged_df['Disease'].unique())


# Define the color map for cell types and shape sequence for diseases this will cause consistent through different filtering
cell_type_colors = {cell_type: color for cell_type, color in zip(unique_cell_types, px.colors.qualitative.Plotly)}

###working and hard coded
# disease_shapes = ['circle', 'square', 'diamond', 'cross', 'x', 'star']
# # Create a dictionary mapping each unique disease to a shape
# disease_shape_dict = dict(zip(unique_diseases, disease_shapes))

available_shapes = ['circle', 'square', 'diamond', 'cross', 'x', 'star', 'triangle-up', 'triangle-down', 'triangle-left', 'triangle-right']

# Create a dictionary mapping each unique disease to a shape
disease_shape_dict = {}
for i, disease in enumerate(unique_diseases):
    shape_index = i % len(available_shapes) # Use modulo to cycle through available shapes
    shape = available_shapes[shape_index]
    disease_shape_dict[disease] = shape

# exit()

##All All
filtered_df = filter_dataframe(merged_df, cell_types=unique_cell_types, diseases=unique_diseases, sample_types=["Meniges","Spinal Cord Tissue"])
# Normalize the data, perform PCA, and create the PCA plot
normalized_data, data = normalize_data(filtered_df)
final_df, pca = perform_pca(normalized_data, data, merged_df)
plot_pca(final_df,cell_type_colors,disease_shape_dict, title="All Samples", file_name="All_samples")


####All Spinal Cord Tissue
# Filter the DataFrame using the lists
filtered_df = filter_dataframe(merged_df, cell_types=unique_cell_types, diseases=unique_diseases, sample_types=["Spinal Cord Tissue"])
# Normalize the data, perform PCA, and create the PCA plot
normalized_data, data = normalize_data(filtered_df)
final_df, pca = perform_pca(normalized_data, data, merged_df)
plot_pca(final_df,cell_type_colors,disease_shape_dict, title="All Spinal Cord Tissue Samples", file_name="All_spinal_cord_tissue_samples")

####All Meniges Samples
# Filter the DataFrame using the lists
filtered_df = filter_dataframe(merged_df, cell_types=unique_cell_types, diseases=unique_diseases, sample_types=["Meniges"])
# Normalize the data, perform PCA, and create the PCA plot
normalized_data, data = normalize_data(filtered_df)
final_df, pca = perform_pca(normalized_data, data, merged_df)
plot_pca(final_df,cell_type_colors,disease_shape_dict, title="All Meniges Tissue Samples", file_name="All_Meniges_samples")


####All Spinal Cord Tissue
# Filter the DataFrame using the lists
filtered_df = filter_dataframe(merged_df, cell_types=["CD45high","CD45low","CD45"], diseases=unique_diseases, sample_types=["Spinal Cord Tissue"])
# Normalize the data, perform PCA, and create the PCA plot
normalized_data, data = normalize_data(filtered_df)
final_df, pca = perform_pca(normalized_data, data, merged_df)
plot_pca(final_df,cell_type_colors,disease_shape_dict, title="All Spinal Cord CD45 Tissue Samples", file_name="All_spinal_cord_CD45_tissue_samples")





























# # Read the Excel file
# df1 = pd.read_excel('EAE_part2.xlsx')

# # Get the list of value columns to be updated (excluding 'lipid', 'Transition', 'type', and 'Blank')
# value_columns = [col for col in df.columns if col not in ['lipid', 'Transition', 'type', 'Blank']]

# # Subtract the 'Blank' column values from each value column and set negative values to 0
# for col in value_columns:
#     df1[col] = df1[col] - df1['Blank']
#     df1.loc[df1[col] < 0, col] = 0

# # Check the updated dataframe
# print(df1.head(None))

# df2 = pd.read_csv('sample_names_for_PCA.csv')

# # Rename the 'Name' column in df2 to match the value columns in df1
# df2 = df2.rename(columns={'Name': 'value_column'})

# # Convert the wide format of df1 to long format
# df1_long = df1.melt(id_vars=['lipid', 'Transition', 'type', 'Blank'], var_name='value_column', value_name='value')

# # Merge the two DataFrames on the value_column
# merged_df = df1_long.merge(df2, on='value_column')

# # Define your criteria for filtering the DataFrame
# criteria = [
#     {'Cell Type': 'T-cell', 'Disease': 'EAE', 'Sample Type': 'Serum'},
#     {'Disease': 'B-cell', 'Disease': 'MS', 'Sample Type': 'Plasma'},
#     {'Sample Type': 'B-cell', 'Disease': 'MS', 'Sample Type': 'Plasma'}
#     # Add more criteria as needed
# ]

# # Create a list of filtered DataFrames based on the criteria
# filtered_dfs = []
# for criterion in criteria:
#     filtered_df = merged_df[(merged_df['Cell Type'] == criterion['Cell Type']) &
#                             (merged_df['Disease'] == criterion['Disease']) &
#                             (merged_df['Sample Type'] == criterion['Sample Type'])]
#     filtered_dfs.append(filtered_df)

# # Check the filtered DataFrames
# for df in filtered_dfs:
#     print(df.head())