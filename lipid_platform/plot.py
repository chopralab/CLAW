import pandas as pd
import plotly.express as px

def plot_ratio(df, color_mapping, output_directory, ratio_threshold=None):
    """
    Plots the ratio of lipids for each unique Sample_ID in the given DataFrame.

    Parameters:
        df (pd.DataFrame): Input DataFrame with columns 'Sample_ID', 'Lipid', and 'ratio'.
        color_mapping (dict): Mapping of patterns to colors for the Lipid values.
        output_directory (str): Directory where to save the plot images.
        ratio_threshold (float, optional): Minimum ratio value for plotting. If provided, rows with a ratio 
                                           value below this threshold will be excluded from plotting.

    Returns:
        None
    """
    import plotly.express as px
    import os

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Apply the ratio threshold filter if provided
    if ratio_threshold is not None:
        df = df[df['Ratio'] >= ratio_threshold]

    # Get the unique Sample_IDs
    sample_ids = df['Sample_ID'].unique()

    # Loop over the unique Sample_IDs
    for sample_id in sample_ids:

        # Filter the dataframe for the current Sample_ID
        df_sample = df[df['Sample_ID'] == sample_id]

        # Assign colors to Lipids based on patterns
        lipid_colors = []
        for lipid in df_sample['Lipid']:
            color = 'gray'  # Default color
            for pattern, pattern_color in color_mapping.items():
                if pattern in lipid:
                    color = pattern_color
                    break
            lipid_colors.append(color)

        # Create the bar plot
        fig = px.bar(df_sample, x='Lipid', y='Ratio', text='Ratio', title=f'Bar Plot for Sample_ID: {sample_id}')

        # Apply colors to the bars
        fig.update_traces(
            marker_color=lipid_colors,
            texttemplate='%{text:.2f}',
            textposition='auto',
            marker_line_width=0
        )

        # Customize the layout
        fig.update_layout(
            uniformtext_minsize=18,
            uniformtext_mode='hide',
            xaxis=dict(
                title='Lipid',
                titlefont=dict(size=16)
            ),
            yaxis=dict(
                title='Ratio',
                titlefont=dict(size=16),
                tickfont=dict(size=16)  # Set the font size of y-axis labels
            ),
            title=dict(
                text=f'Sample_ID: {sample_id}',
                font=dict(size=20)  # Set the title font size
            )
        )
        
        # Save the plot as an image
        file_name = os.path.join(output_directory, f"plot_{sample_id}.png")

        # Check if the file already exists
        index = 1
        while os.path.exists(file_name):
            file_name = os.path.join(output_directory, f"plot_{sample_id}_{index}.png")
            index += 1

        fig.write_image(file_name)




def printed_ratio(df_OzESI_ratio_sort):
    """
    Prints the Lipid, Sample_ID, db_pos, and ratio for each row in the given DataFrame.

    Parameters:
        df_OzESI_ratio_sort (pd.DataFrame): DataFrame with columns 'Lipid', 'Sample_ID', 'db_pos', and 'ratio'.

    Returns:
        None
    """
    # Iterate through each row in the DataFrame
    for index, row in df_OzESI_ratio_sort.iterrows():
        # Extract Lipid, Sample_ID, Labels and ratio from the row
        lipid = row['Lipid']
        sample_id = row['Sample_ID']
        db_pos = row['db_pos']
        ratio = row['Ratio']

        # Check if ratio is not NaN
        if not pd.isna(ratio):
            # Print out the values
            print(f'Lipid: {lipid}, Sample_ID: {sample_id}, db_pos: {db_pos}, Ratio: {ratio}')
            



import os
import pandas as pd
import plotly.graph_objects as go

def plot_chromatogram(file_path, plot_path, plot_name, x_range=None):
    # Read the CSV file
    data = pd.read_csv(file_path, skiprows=1) # Skipping the first row with metadata

    # If an x_range is provided, filter the data accordingly
    if x_range:
        start_time, end_time = x_range
        data = data[(data['X(Minutes)'] >= start_time) & (data['X(Minutes)'] <= end_time)]

    # Create the plot
    fig = go.Figure()

    # Add the line trace
    fig.add_trace(go.Scatter(x=data['X(Minutes)'], y=data['Y(Counts)'], mode='lines', name='Intensity'))

    # Set the title and axis labels
    fig.update_layout(
        title="Chromatogram of Canola Oil Crude Sample",
        xaxis_title="Time (minutes)",
        yaxis_title="Intensity (Counts)",
        font=dict(
            family="Arial, sans-serif",
            size=18
        ),
        showlegend=False
    )

    # Customize the plot appearance
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey', linecolor='Black', linewidth=2, mirror=True)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey', linecolor='Black', linewidth=2, mirror=True,
                    exponentformat='e', showexponent='all')


    # Check if the plot path exists, if not then create it
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    # Check if the file exists, if it does add 1 and try again
    count = 1
    filename = plot_name
    while os.path.exists(f"{plot_path}{filename}"):
        filename = f"{plot_name[:-4]}_{count}{plot_name[-4:]}"
        count += 1

    # Save the image
    fig.write_image(f"{plot_path}{filename}", scale=2)

    # Show the plot
    fig.show()



