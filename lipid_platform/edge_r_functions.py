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



def prep_edge_R2(merged_df,json_list_pairs,Pre_edge_r_path,blank_name,labels_list):
    for i in range(len(json_list_pairs)):
        json1 = json_list_pairs[i][0]
        json2 = json_list_pairs[i][1]

        json_blank = {"Sample Name":[blank_name]}

        filtered_df1 = filter_dataframe(merged_df,json1)
        filtered_df2 = filter_dataframe(merged_df,json2)
        filtered_blank = filter_dataframe(merged_df,json_blank)



        filtered_df1 = filtered_df1.groupby(labels_list)['Intensity'].mean().reset_index()
        filtered_df2 = filtered_df2.groupby(labels_list)['Intensity'].mean().reset_index()
        filtered_blank = filtered_blank.groupby(labels_list)['Intensity'].mean().reset_index()



        reformatted_df1 = filtered_df1.pivot(index=['Lipid', 'Class'], columns='Sample Name', values='Intensity').reset_index()
        reformatted_df2 = filtered_df2.pivot(index=['Lipid', 'Class'], columns='Sample Name', values='Intensity').reset_index()
        reformatted_blank = filtered_blank.pivot(index=['Lipid', 'Class'], columns='Sample Name', values='Intensity').reset_index()

        num_value_columns_df1 = reformatted_df1.shape[1] - 2  # subtract 2 for 'Lipid' and 'Class'
        num_value_columns_df2 = reformatted_df2.shape[1] - 2  # subtract 2 for 'Lipid' and 'Class'



        combined_df = reformatted_df1.merge(reformatted_df2, on=['Lipid', 'Class'], how='inner')

        combined_df = combined_df.merge(reformatted_blank, on=['Lipid', 'Class'], how='inner')




        ##Aquiring Names ##Add an extra string option

        title1 = json_to_string(json1)
        title2 = json_to_string(json2)

        title = title1 +" vs "+title2
        title = title.replace(" | ","__")
        length1 = num_value_columns_df1
        length2 = num_value_columns_df2

        combined_df["Title1"] = title1
        combined_df["Title2"] = title2
        combined_df["Title"] = title
        combined_df["length1"] = length1
        combined_df["length2"] = length2
        combined_df["Blank_name"] = blank_name




        combined_df.to_csv(Pre_edge_r_path+title+".csv",index=False)
    return combined_df

