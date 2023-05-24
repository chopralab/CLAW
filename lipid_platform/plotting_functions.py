#Import all the necessary libraries
import pymzml
import csv
import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
import re
import plotly.express as px
from collections import defaultdict

import plotly.io as pio
import json
import plotly.graph_objs as go
import matplotlib.colors as mcolors

import json
import ipywidgets as widgets
from IPython.display import display

import warnings
from IPython.display import display, Image, clear_output
import time



# Function to convert HEX colors to RGBA with specified alpha value, and then back to HEX
def hex_to_rgba_hex(hex_color, alpha=1):
    rgba_color = list(mcolors.hex2color(hex_color))
    rgba_color.append(alpha)
    return mcolors.to_hex(rgba_color, keep_alpha=True)


lipid_classes = ["CAR", "CE", "Cer", "FA", "PC", "PE", "PG", "PI", "PS", "SM", "TAG",'DAG','TAG | DAG','DAG | CE','TAG | DAG | CE']
lipid_colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#808080", "#cab2d6", "#6a3d9a",'#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3']

lipid_class_colors = dict(zip([lipid.upper() for lipid in lipid_classes], lipid_colors))


# Update colors with 0.5 alpha
alpha = 0.5
transparent_colors = [hex_to_rgba_hex(color, alpha) for color in lipid_colors]
lipid_class_colors_alpha = dict(zip([lipid.upper() for lipid in lipid_classes], transparent_colors))
    


def filter_dataframe(df, json_filter):
    for column, values in json_filter.items():
        df = df[df[column].isin(values)]
    return df


# Function to convert HEX colors to RGBA with specified alpha value, and then back to HEX
def hex_to_rgba_hex(hex_color, alpha=1):
    rgba_color = list(mcolors.hex2color(hex_color))
    rgba_color.append(alpha)
    return mcolors.to_hex(rgba_color, keep_alpha=True)

def json_to_string(json_dict):
    result = []
    for key, values in json_dict.items():
        values_str = ' '.join(values)
        result.append(f"{key}: {values_str}")
    return ' | '.join(result)


def make_pie_chart_no_replicates(merged_df,save_path,json_list_singles,labels_list,blank_name,show_percentages=False):

    ##Maybe filter?
    for i in range(len(json_list_singles)):
        plot_title = 'Sum of Intensity by Class for' + json_to_string(json_list_singles[i])
        save_name = save_path+json_to_string(json_list_singles[i]).replace(" | ","__")
        filtered_df = filter_dataframe(merged_df,json_list_singles[i])
        filtered_df = filtered_df[filtered_df["Sample Name"] != blank_name]

        groupby_columns = labels_list

        # Aggregate functions to apply
        aggregations = {
            'Intensity': 'mean',
            'Blank Subtraction': 'mean'
        }

        # Group the DataFrame and apply the aggregation functions

        filtered_df_grouped = filtered_df.groupby(groupby_columns).agg(aggregations).reset_index()
        grouped_df_sum_avg = filtered_df_grouped.groupby('Class')['Blank Subtraction'].sum().reset_index()


        textinfo = 'label+percent' if show_percentages else "none"




        fig1 = px.pie(grouped_df_sum_avg, values='Blank Subtraction', names='Class', title=plot_title)
        fig1.update_traces(marker=dict(colors=[lipid_class_colors_alpha[lipid.upper()] for lipid in grouped_df_sum_avg['Class']]),
                        hovertemplate='%{label}: %{percent:.2%}', textinfo=textinfo)
        pio.write_html(fig1,  save_name + "Sum Pie.html")

        # Save the plot as plot.png
        pio.write_image(fig1, save_name + "Sum Pie.svg")
    return



def average_pie_chart_no_repeats(merged_df,save_path,json_list_singles,labels_list,blank_name,show_percentages=False):

    for i in range(len(json_list_singles)):
        plot_title = 'Mean of Intensity by Class for' + json_to_string(json_list_singles[i])
        save_name = save_path+json_to_string(json_list_singles[i]).replace(" | ","__")
        filtered_df = filter_dataframe(merged_df,json_list_singles[i])
        filtered_df = filtered_df[filtered_df["Sample Name"] != blank_name]

        groupby_columns = labels_list

        # Aggregate functions to apply
        aggregations = {
            'Intensity': 'mean',
            'Blank Subtraction': 'mean'
        }

        # Group the DataFrame and apply the aggregation functions

        filtered_df_grouped = filtered_df.groupby(groupby_columns).agg(aggregations).reset_index()
        grouped_df_avg = filtered_df_grouped.groupby('Class')['Blank Subtraction'].mean().reset_index()


        textinfo = 'label+percent' if show_percentages else "none"




        fig1 = px.pie(grouped_df_avg, values='Blank Subtraction', names='Class', title=plot_title)
        fig1.update_traces(marker=dict(colors=[lipid_class_colors_alpha[lipid.upper()] for lipid in grouped_df_avg['Class']]),
                        hovertemplate='%{label}: %{percent:.2%}', textinfo=textinfo)
        pio.write_html(fig1,  save_name + "Average Pie.html")

        # Save the plot as plot.png
        pio.write_image(fig1, save_name + "Average Pie .svg")
    return


def make_bar_plot_comparisons(merged_df, save_path, json_list_pairs,labels_list,blank_name):

#     merged_df = merged_df[merged_df["Name"] != blank_name]

    ###This plot makes relative bar plots of comparison for total Intensity and average intensity
    ###need to have a better group by column method list to add perhaps with the propper added names from first parser
    ##Should be Lipid, Class and other important info from labeling dataframe
    
    for i in range(len(json_list_pairs)):
        json1 = json_list_pairs[i][0]
        json2 = json_list_pairs[i][1]
        custom_name1 = json_to_string(json1)
        custom_name2 = json_to_string(json2)
        plot_tile = custom_name1+ " vs "+custom_name2
        save_name  = save_path+plot_tile.replace(" | ","__")

#         merged_df = merged_df[merged_df["Name"] != blank_name]

        groupby_columns = labels_list

        aggregations = {
            'Intensity': 'mean',
            'Blank Subtraction': 'mean'
        }
        merged_df1 = merged_df[merged_df["Sample Name"] != blank_name]
        merged_df1 = merged_df1.groupby(groupby_columns).agg(aggregations).reset_index()

        filtered_df1 = filter_dataframe(merged_df1, json1)
        filtered_df2 = filter_dataframe(merged_df1, json2)
        
#         filtered_df1 = filtered_df1[filtered_df1["Sample Name"] != blank_name]
#         filtered_df2 = filtered_df2[filtered_df2["Sample Name"] != blank_name]
        
        filtered_df1['Class'] = filtered_df1['Class'].str.upper()
        filtered_df2['Class'] = filtered_df2['Class'].str.upper()

        filtered_df1_avg = filtered_df1.groupby('Class')['Blank Subtraction'].mean().reset_index()
        filtered_df2_avg = filtered_df2.groupby('Class')['Blank Subtraction'].mean().reset_index()



        combined_max = filtered_df1_avg.set_index('Class')[['Blank Subtraction']].combine(filtered_df2_avg.set_index('Class')[['Blank Subtraction']], np.maximum)


        normalized_df1_avg = filtered_df1_avg.set_index('Class').divide(combined_max)
        normalized_df2_avg = filtered_df2_avg.set_index('Class').divide(combined_max)

        trace1 = go.Bar(
            x=normalized_df1_avg.index,
            y=normalized_df1_avg['Blank Subtraction'],
            name=custom_name1,
            marker=dict(color='red'),
        )

        trace2 = go.Bar(
            x=normalized_df2_avg.index,
            y=normalized_df2_avg['Blank Subtraction'],
            name=custom_name2,
            marker=dict(color='blue'),
        )

        layout = go.Layout(
            title=plot_tile,
            xaxis=dict(title="Class"),
            yaxis=dict(title="Normalized Mean Intensity Blank Subtraction"),
            barmode="group",
        )

        fig = go.Figure(data=[trace1, trace2], layout=layout)

        pio.write_image(fig, save_name + "bar.svg")
        pio.write_html(fig, save_name + "bar.html")
    return

