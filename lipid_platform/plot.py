def plot_ratios(df, color_mapping, output_directory):
    """
    Plots the ratios of lipids for each unique Sample_ID in the given DataFrame.

    Parameters:
        df (pd.DataFrame): Input DataFrame with columns 'Sample_ID', 'Lipid', and 'Ratios'.
        color_mapping (dict): Mapping of patterns to colors for the Lipid values.
        output_directory (str): Directory where to save the plot images.

    Returns:
        None
    """
    import plotly.express as px
    import os

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

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
        fig = px.bar(df_sample, x='Lipid', y='Ratios', text='Ratios', title=f'Bar Plot for Sample_ID: {sample_id}')

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
                title='Ratios',
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
