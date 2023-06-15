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
        result.append(f"{key}_ {values_str}")
    return ' | '.join(result)


def make_pie_chart_no_replicates(merged_df,save_path,json_list_singles,labels_list,blank_name,extra_name,show_percentages=False):

    ##Maybe filter?
    for i in range(len(json_list_singles)):
        plot_title = 'Sum of Intensity by Class for' + json_to_string(json_list_singles[i])+ extra_name
        save_name = save_path+json_to_string(json_list_singles[i]).replace(" | ","__")+ extra_name
        filtered_df = filter_dataframe(merged_df,json_list_singles[i])
        filtered_df = filtered_df[filtered_df["Sample Name"] != blank_name]

        groupby_columns = labels_list

        # Aggregate functions to apply
        aggregations = {
            'Intensity': 'mean',
            'Blank Subtraction': 'mean'
        }

        # Group the DataFrame and apply the aggregation functions

        # filtered_df_grouped = filtered_df.groupby(groupby_columns).agg(aggregations).reset_index()
        filtered_df_grouped = filtered_df
        grouped_df_sum_avg = filtered_df_grouped.groupby('Class')['Blank Subtraction'].sum().reset_index()


        textinfo = 'label+percent' if show_percentages else "none"




        fig1 = px.pie(grouped_df_sum_avg, values='Blank Subtraction', names='Class', title=plot_title)
        fig1.update_traces(marker=dict(colors=[lipid_class_colors_alpha[lipid.upper()] for lipid in grouped_df_sum_avg['Class']]),
                        hovertemplate='%{label}: %{percent:.2%}', textinfo=textinfo)
        pio.write_html(fig1,  save_name + "Sum Pie.html")

        # Save the plot as plot.png
        pio.write_image(fig1, save_name + "Sum Pie.svg")
    return



def average_pie_chart_no_repeats(merged_df,save_path,json_list_singles,labels_list,blank_name,extra_name,show_percentages=False):

    for i in range(len(json_list_singles)):
        plot_title = 'Mean of Intensity by Class for' + json_to_string(json_list_singles[i])+ extra_name
        save_name = save_path+json_to_string(json_list_singles[i]).replace(" | ","__") + extra_name
        filtered_df = filter_dataframe(merged_df,json_list_singles[i])
        filtered_df = filtered_df[filtered_df["Sample Name"] != blank_name]

        groupby_columns = labels_list

        # Aggregate functions to apply
        aggregations = {
            'Intensity': 'mean',
            'Blank Subtraction': 'mean'
        }

        # Group the DataFrame and apply the aggregation functions

        # filtered_df_grouped = filtered_df.groupby(groupby_columns).agg(aggregations).reset_index()
        filtered_df_grouped = filtered_df
        grouped_df_avg = filtered_df_grouped.groupby('Class')['Blank Subtraction'].mean().reset_index()


        textinfo = 'label+percent' if show_percentages else "none"




        fig1 = px.pie(grouped_df_avg, values='Blank Subtraction', names='Class', title=plot_title)
        fig1.update_traces(marker=dict(colors=[lipid_class_colors_alpha[lipid.upper()] for lipid in grouped_df_avg['Class']]),
                        hovertemplate='%{label}: %{percent:.2%}', textinfo=textinfo)
        pio.write_html(fig1,  save_name + "Average Pie.html")

        # Save the plot as plot.png
        pio.write_image(fig1, save_name + "Average Pie .svg")
    return


def make_bar_plot_comparisons(merged_df, save_path, json_list_pairs,labels_list,blank_name, extra_name):

#     merged_df = merged_df[merged_df["Name"] != blank_name]

    ###This plot makes relative bar plots of comparison for total Intensity and average intensity
    ###need to have a better group by column method list to add perhaps with the propper added names from first parser
    ##Should be Lipid, Class and other important info from labeling dataframe
    
    for i in range(len(json_list_pairs)):
        json1 = json_list_pairs[i][0]
        json2 = json_list_pairs[i][1]
        custom_name1 = json_to_string(json1)
        custom_name2 = json_to_string(json2)
        plot_tile = custom_name1+ " vs "+custom_name2+ extra_name
        save_name  = save_path+plot_tile.replace(" | ","__")+ extra_name

#         merged_df = merged_df[merged_df["Name"] != blank_name]

        groupby_columns = labels_list

        aggregations = {
            'Intensity': 'mean',
            'Blank Subtraction': 'mean'
        }
        merged_df1 = merged_df[merged_df["Sample Name"] != blank_name]
        # merged_df1 = merged_df1.groupby(groupby_columns).agg(aggregations).reset_index()

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


#Function to read in MRM database
#Option to remove STDs from database##Not finished need option to use another database with no qualitative ACs


def read_mrm_list(filename,remove_std = True,custom_data=False):
    if custom_data==False:
        mrm_list_new = pd.read_excel(filename, sheet_name=None)
        mrm_list_new = pd.concat(mrm_list_new, ignore_index=True)
        mrm_list_offical = mrm_list_new[['Compound Name', 'Parent Ion', 'Product Ion', 'Class']]
        # Add underscore to middle of columns names
        mrm_list_offical.columns = mrm_list_offical.columns.str.replace(' ', '_')
        # Round Parent Ion and Product Ion to 1 decimal place
        mrm_list_offical['Parent_Ion'] = np.round(mrm_list_offical['Parent_Ion'],1)
        mrm_list_offical['Product_Ion'] = np.round(mrm_list_offical['Product_Ion'],1)
        # Create transition column by combining Parent Ion and Product Ion with arrow between numbers
        mrm_list_offical['Transition'] = mrm_list_offical['Parent_Ion'].astype(str) + ' -> ' + mrm_list_offical['Product_Ion'].astype(str)
        # Change column compound name to lipid
        mrm_list_offical = mrm_list_offical.rename(columns={'Compound_Name': 'Lipid'})
        # Make a column called Class match lipid column to lipid types
        
        if remove_std == True:
            lipid_class = mrm_list_offical['Class'].unique()
            lipid_class_to_keep = ['PS','PG','CE','PC', 'DAG', 'PE', 'TAG', 'FA', 'Cer', 'CAR', 'PI','SM']
            mrm_list_offical = mrm_list_offical[mrm_list_offical['Class'].isin(lipid_class_to_keep)]

    else:
        mrm_list_offical = pd.read_csv('lipid_database/Lipids_database_old.csv')

        # Round Parent Ion and Product Ion to 1 decimal place
        mrm_list_offical['Parent_Ion'] = np.round(mrm_list_offical['Parent_Ion'],1)
        mrm_list_offical['Product_Ion'] = np.round(mrm_list_offical['Product_Ion'],1)
        # Create transition column by combining Parent Ion and Product Ion with arrow between numbers
        mrm_list_offical['Transition'] = mrm_list_offical['Parent_Ion'].astype(str) + ' -> ' + mrm_list_offical['Product_Ion'].astype(str)
    return mrm_list_offical



def mzml_parser(file_name):
    df = pd.DataFrame(columns=['Lipid','Parent_Ion','Product_Ion','Intensity','Transition','Class','Sample_ID'])
    data_folder = os.listdir(file_name) #Path to the mzml files
    data_folder.sort()
    path_to_mzml_files = file_name
    for file in data_folder:
            if file.endswith('.mzML'):

                    run = pymzml.run.Reader(path_to_mzml_files+file, skip_chromatogram=False) #Load the mzml file into the run object



                    df_all = pd.DataFrame(columns=['Lipid','Parent_Ion','Product_Ion','Intensity','Transition','Class','Sample_ID']) #Create empty pandas dataframe to store the data

                    #create pandas dataframe to store the data with the columns Parent Ion, Product Ion, Intensity, Transition Lipid and Class
                   
                    q1_mz = 0 #Create empty variables to store the Q1 and Q3 m/z values
                    q3_mz = 0
                    count = 0 #Create a counter to keep track of the number of transitions
                    for spectrum in run:


                            for element in spectrum.ID.split(' '):
                                    intensity_store = np.array([])
                                    if 'Q1' in element:
                                            q1 = element.split('=')
                                            q1_mz= np.round((float(q1[1])),1)

                                    if 'Q3' in element:
                                
                                            q3 = element.split('=')
  
                                            q3_mz=np.round(float(q3[1]),1)


                                            for mz,intensity in spectrum.peaks(): #Get the m/z and intensity values from the spectrum
                                                    intensity_store = np.append(intensity_store,intensity) #Store the intensity values in an array



                                    if 'Q3' in element:
                                            # print(intensity_sum)
                                            intensity_sum = np.sum(intensity_store) #Sum the intensity values
                                            df_all.loc[count,'Parent_Ion'] = q1_mz #Store the Q1 and Q3 m/z values in the pandas dataframe
                                            df_all.loc[count,'Product_Ion'] = q3_mz
                                            #round the Q1 and Q3 m/z values to 1 decimal places
                                            df_all.loc[count,'Parent_Ion'] = np.round(df_all.loc[count,'Parent_Ion'],1)
                                            df_all.loc[count,'Product_Ion'] = np.round(df_all.loc[count,'Product_Ion'],1)
                                            df_all.loc[count,'Intensity'] = intensity_sum #Store the intensity values in the pandas dataframe
                                            df_all.loc[count,'Transition'] = str(q1_mz)+ ' -> '+ str(q3_mz) #Store the transition values in the pandas dataframe
                                            #add file name to Sample_ID column without the mzmL extension
                                            df_all.loc[count,'Sample_ID'] = file[:-5]
                                            count+=1

            #append df_all to df
            df = df.append(df_all, ignore_index=True)
    return df

# Function to create an ion dictionary from an MRM database DataFrame
def create_ion_dict(mrm_database):
    ion_dict = defaultdict(list)
    # Iterate through the rows of the MRM database DataFrame
    for index, row in mrm_database.iterrows():
        # Add a tuple with Lipid and Class to the ion dictionary using Parent_Ion and Product_Ion as the key
        ion_dict[(row['Parent_Ion'], row['Product_Ion'])].append((row['Lipid'], row['Class']))
    return ion_dict

# Function to check if the absolute difference between two values is within a given tolerance
def within_tolerance(a, b, tolerance=0.1):
    return abs(a - b) <= tolerance

# Function to match the ions in a DataFrame row with the ions in an ion dictionary
def match_ions(row, ion_dict, tolerance=0.1):
    ions = (row['Parent_Ion'], row['Product_Ion'])
    matched_lipids = []
    matched_classes = []

    # Iterate through the ion dictionary
    for key, value in ion_dict.items():
        # Check if both the Parent_Ion and Product_Ion values are within the specified tolerance
        if within_tolerance(ions[0], key[0], tolerance) and within_tolerance(ions[1], key[1], tolerance):
            # If within tolerance, extend the matched_lipids and matched_classes lists with the corresponding values
            matched_lipids.extend([match[0] for match in value])
            matched_classes.extend([match[1] for match in value])

    # If any matches were found, update the Lipid and Class columns in the row
    if matched_lipids and matched_classes:
        row['Lipid'] = ' | '.join(matched_lipids)
        row['Class'] = ' | '.join(matched_classes)

    return row

####Combined functions for Matching

def match_lipids_parser(mrm_database,df, tolerance=0.3):
    ion_dict = create_ion_dict(mrm_database)
    # Assuming you have the df DataFrame to apply the match_ions function
    df_matched = df.apply(lambda row: match_ions(row, ion_dict=ion_dict, tolerance=tolerance), axis=1)


    df_matched = df_matched.dropna()
    
    return df_matched


def save_dataframe(df, folder_name, file_name, max_attempts=5):
    folder_path = f'{folder_name}'
    os.makedirs(folder_path, exist_ok=True)

    for i in range(max_attempts):
        file_path = f'{folder_path}/{file_name}.csv'
        if not os.path.isfile(file_path):
            df.to_csv(file_path, index=False)
            print(f"Saved DataFrame to {file_path}")
            break
    else:
        print(f"Failed to save DataFrame after {max_attempts} attempts.")
        return None

##Adds labels and method type
def add_labels(labels_df,matched_df):
    for _, group_row in labels_df.iterrows():
        for index, burda_row in matched_df.iterrows():
            if group_row['Sample Name'].lower() in burda_row['Sample_ID'].lower():
                # Add group_row data to the corresponding burda_row if the condition is met
                for col in labels_df.columns:
                    if col not in matched_df.columns:
                        matched_df[col] = None
                    matched_df.at[index, col] = group_row[col]

    matched_df['method_type'] = matched_df['Sample_ID'].apply(lambda x: x.split('_')[0])

    
    return matched_df
    
# def subtract_blank(labels_df,matched_df,remove_list,blank_name):
def subtract_blank(labels_df,matched_df,blank_name):
    ###Removing Stuff
#     for i in remove_list:
#         matched_df = matched_df[matched_df['Sample Name'] != i]
    
    matched_df = matched_df.dropna()

    ###Subtracting the Blank
    blank_intensities_df = matched_df[matched_df['Sample Name'] == blank_name][['Lipid',"Transition", 'Intensity', 'method_type']]
    blank_intensities_df.columns = ['Lipid',"Transition", 'Blank_Intensity', 'method_type']


# Get the unique names in the Name column, excluding the blank_name
    unique_names = matched_df['Sample Name'].unique()
    # print(unique_names)
    ###Keep this line to drop the blank
    # unique_names = unique_names[unique_names != blank_name] ###This drops the blank
    merged_df = pd.DataFrame()  
    for name in unique_names:
        temp_df = matched_df[matched_df['Sample Name'] == name]
#         print(len(temp_df),"TEMP")
        numbers = np.array((temp_df["Intensity"] ))
        numbers1 = np.array((blank_intensities_df['Blank_Intensity'] ))
#         print(len(temp_df),"TEMP")
        numbers2 = numbers - numbers1
        numbers2[numbers2<0] = 0
        # Merge the blank intensities DataFrame with the temporary DataFrame
        temp_df["Blank Subtraction"] = numbers2
#         print(len(temp_df),"TEMP")
        merged_df = merged_df.append(temp_df)



    merged_df['Class'] = merged_df['Class'].replace({'TAG | TAG': 'TAG', 'FA | FA': 'FA'})

    return merged_df
def add_subclass_and_length(merged_df):
    merged_df = merged_df.reset_index(drop=True)
    merged_df.loc[merged_df['Class'] != 'Cer', 'Lipid'] = merged_df.loc[merged_df['Class'] != 'Cer', 'Lipid'].str.replace(',', '|')




    # merged_df = merged_df[merged_df['Class'] != 'PG']
    # merged_df = merged_df[merged_df['Class'] != 'PE']
    # merged_df = merged_df[merged_df['Class'] != 'PI']
    # merged_df = merged_df[merged_df['Class'] != 'PS']
    # merged_df = merged_df[merged_df['Class'] != 'SM']

    # print(set(list(merged_df["Class"])))
    # print(list(set(list(merged_df["Class"]))))

    subclasses = []
    chain_length = []
    saturation = []
    between_parenthese = []
    FA_Chain = []
    ##Would need to recheck for CE -O and all other classes if tolerance or MRMs are added
    ##Should add natural Multiples and replace the , with a |

    string = "DG(24:0)_C16:0"
    print(string[string.find("_C")+2:])

    print(len(merged_df))
    # exit()

    for i in range(len(merged_df)):
        jj = merged_df["Class"][i]
        xx = merged_df["Lipid"][i]
        subclass = ""

        if jj == "FA": ##Maybe add chain length and saturation 
            subclasses.append("FA")

            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]

                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append(between_parents_temp)


            else:


                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append(xx[xx.find("("):xx.find(")")+1])


        elif jj == "CAR":
            FA_Chain.append("None")
            subclasses.append("CAR")

            if "(" in xx:
                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])

            else:
                chain_length.append("0")
                saturation.append("0")    
                between_parenthese.append("None")


        ###What about O????
        ######
        elif jj == "CE":
            FA_Chain.append("None")

            if "O2" not in xx:
                subclasses.append("CE")

                # if "(" in xx:
                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])

            elif "O2" in xx:
                subclasses.append("CE-O")

                # if "(" in xx:
                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(";")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])



        elif jj == "Cer":
            FA_Chain.append("None")

            if "CerP" in xx:
                subclasses.append("CerP")
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                saturation.append("NA")

                chain_length.append("NA")

            elif "1-O-" in xx:
                subclasses.append("1-O-Cer")
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                saturation.append("NA")
                chain_length.append("NA")


            elif "omega-linoleoyloxy-GlcCer(" in xx:
                subclasses.append("omega-linoleoyloxy-GlcCer")
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                saturation.append("NA")
                chain_length.append("NA")


            elif "omega-linoleoyloxy-Cer" in xx:

                subclasses.append("omega-linoleoyloxy-Cer")
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                saturation.append("NA")

                chain_length.append("NA")

            else:
                subclasses.append("Cer")
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                saturation.append("NA")

                chain_length.append("NA")

        elif jj == "DAG" or jj== "DAG | CE" or jj=="TAG | DAG" or jj == "TAG | DAG | CE":
            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""
                FA_Chain_temp = ""
                temp_sub_class = ""
                for i in temp_lipid:
                    if "DG" in i:
                        if len(temp_string_length) >0:
                            temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                            temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                            between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                            FA_Chain_temp = FA_Chain_temp+" | "+ i[i.find("_C")+2:]
                            if "DG(dO" in i:
                                temp_sub_class = temp_sub_class +" | "+ "DG-DO"
                            elif "DG(" in i:
                                temp_sub_class = temp_sub_class +" | "+"DG"
                            elif "DG(P" in i:
                                temp_sub_class = temp_sub_class +" | "+"DG-P"
                            else:
                                temp_sub_class = temp_sub_class +" | "+ "DG-O"      
                        else:

                            temp_string_length = i[i.find("(")+1:i.find(":")]
                            temp_string_saturation = i[i.find(":")+1:i.find(")")]
                            between_parents_temp = i[i.find("("):i.find(")")+1]
                            FA_Chain_temp = i[i.find("_C")+2:]
                            if "DG(dO" in i:
                                temp_sub_class = "DG-DO"

                            elif "DG(P" in i:
                                temp_sub_class = "DG-P"
                            elif "DG(O" in i:
                                temp_sub_class = "DG-O" 
                            else:
                                temp_sub_class = "DG"

                    elif "CE" in i:                            
                        if len(temp_string_length) >0:
                            temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                            temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                            between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                            FA_Chain_temp = FA_Chain_temp+" | "+ "NA"
                        else:
                            temp_string_length =  i[i.find("(")+1:i.find(":")]
                            temp_string_saturation = i[i.find(":")+1:i.find(")")]
                            between_parents_temp = i[i.find("("):i.find(")")+1]
                            FA_Chain_temp =  "NA"

                    elif "TG" in i:
        ##Add TAG for this one after I do individual TAG
                        if len(temp_string_length) >0:
                            temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                            temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                            between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                            FA_Chain_temp = FA_Chain_temp+" | "+ i[i.find("_FA")+3:]
                            if "TG(O" in i:
                                temp_sub_class = temp_sub_class +" | "+ "TG(O"
                            else:
                                temp_sub_class = temp_sub_class +" | "+"TG"

                        else:

                            temp_string_length = i[i.find("(")+1:i.find(":")]
                            temp_string_saturation = i[i.find(":")+1:i.find(")")]
                            between_parents_temp = i[i.find("("):i.find(")")+1]
                            FA_Chain_temp = i[i.find("_FA")+3:]
                            if "TG(O" in i:
                                temp_sub_class = "TG(O"
                            else:
                                temp_sub_class = "TG"



                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append(FA_Chain_temp)
                subclasses.append(temp_sub_class)


            else:


                # exit()
                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append(xx[xx.find("_C")+2:])

                if "DG(dO" in xx:
                    temp_sub_class = "DG-DO"

                elif "DG(P" in xx:
                    temp_sub_class = "DG-P"
                elif "DG(O" in xx:
                    temp_sub_class = "DG-O" 
                else:
                    temp_sub_class = "DG"
                subclasses.append(temp_sub_class)


        elif jj == "TAG":

            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""
                FA_Chain_temp = ""
                temp_sub_class = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                        FA_Chain_temp = FA_Chain_temp+" | "+ i[i.find("_FA")+3:]
                        if "TG(O" in i:
                            temp_sub_class = temp_sub_class +" | "+ "TG(O"
                        else:
                            temp_sub_class = temp_sub_class +" | "+"TG"

                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]
                        FA_Chain_temp = i[i.find("_FA")+3:]
                        if "TG(O" in i:
                            temp_sub_class = "TG(O"
                        else:
                            temp_sub_class = "TG"

                subclasses.append(temp_sub_class)
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append(FA_Chain_temp)


            else:

                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append(xx[xx.find("_FA")+3:])

                if "TG(O" in xx:
                    temp_sub_class = "TG(O"
                else:
                    temp_sub_class = "TG"
                subclasses.append(temp_sub_class)

        elif jj == "PC":
            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""

                temp_sub_class = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                        if "LPC(O-" in i:
                            temp_sub_class =temp_sub_class +" | "+ "LPC(O-"
                        elif "LPC(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPC(P-"    
                        elif "LPC" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPC"    
                        elif "PC(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PC(P-"
                        elif "PC(O-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PC(O-"      
                        else:
                            temp_sub_class = temp_sub_class +" | "+"PC"

                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]

                        if "LPC(O-" in i:
                            temp_sub_class = "LPC(O-"
                        elif "LPC(P-" in i:
                            temp_sub_class = "LPC(P-"    
                        elif "LPC" in i:
                            temp_sub_class = "LPC"    
                        elif "PC(P-" in i:
                            temp_sub_class = "PC(P-"
                        elif "PC(O-" in i:
                            temp_sub_class = "PC(O-"      
                        else:
                            temp_sub_class = "PC"

                subclasses.append(temp_sub_class)
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append("NA")


            else:

                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append("NA")
                if "LPC(O-" in xx:
                    temp_sub_class = "LPC(O-"
                elif "LPC(P-" in xx:
                    temp_sub_class = "LPC(P-"    
                elif "LPC" in xx:
                    temp_sub_class = "LPC"    
                elif "PC(P-" in xx:
                    temp_sub_class = "PC(P-"
                elif "PC(O-" in xx:
                    temp_sub_class = "PC(O-"      
                else:
                    temp_sub_class = "PC"
                subclasses.append(temp_sub_class)

        elif jj == "PE":
            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""

                temp_sub_class = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                        if "LPE(O-" in i:
                            temp_sub_class =temp_sub_class +" | "+ "LPE(O-"
                        elif "LPE(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPE(P-"    
                        elif "LPE" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPE"    
                        elif "PE(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PE(P-"
                        elif "PE(O-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PE(O-"      
                        else:
                            temp_sub_class = temp_sub_class +" | "+"PE"

                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]

                        if "LPE(O-" in i:
                            temp_sub_class = "LPE(O-"
                        elif "LPE(P-" in i:
                            temp_sub_class = "LPE(P-"    
                        elif "LPE" in i:
                            temp_sub_class = "LPE"    
                        elif "PE(P-" in i:
                            temp_sub_class = "PE(P-"
                        elif "PE(O-" in i:
                            temp_sub_class = "PE(O-"      
                        else:
                            temp_sub_class = "PE"

                subclasses.append(temp_sub_class)
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append("NA")


            else:

                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append("NA")
                if "LPE(O-" in xx:
                    temp_sub_class = "LPE(O-"
                elif "LPE(P-" in xx:
                    temp_sub_class = "LPE(P-"    
                elif "LPE" in xx:
                    temp_sub_class = "LPE"    
                elif "PE(P-" in xx:
                    temp_sub_class = "PE(P-"
                elif "PE(O-" in xx:
                    temp_sub_class = "PE(O-"      
                else:
                    temp_sub_class = "PE"
                subclasses.append(temp_sub_class)

        elif jj == "PG":
            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""

                temp_sub_class = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                        if "LPG(O-" in i:
                            temp_sub_class =temp_sub_class +" | "+ "LPG(O-"
                        elif "LPG(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPG(P-"    
                        elif "LPG" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPG"    
                        elif "PG(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PG(P-"
                        elif "PG(O-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PG(O-"      
                        else:
                            temp_sub_class = temp_sub_class +" | "+"PG"

                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]

                        if "LPG(O-" in i:
                            temp_sub_class = "LPG(O-"
                        elif "LPG(P-" in i:
                            temp_sub_class = "LPG(P-"    
                        elif "LPG" in i:
                            temp_sub_class = "LPG"    
                        elif "PG(P-" in i:
                            temp_sub_class = "PG(P-"
                        elif "PG(O-" in i:
                            temp_sub_class = "PG(O-"      
                        else:
                            temp_sub_class = "PG"

                subclasses.append(temp_sub_class)
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append("NA")


            else:

                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append("NA")
                if "LPG(O-" in xx:
                    temp_sub_class = "LPG(O-"
                elif "LPG(P-" in xx:
                    temp_sub_class = "LPG(P-"    
                elif "LPG" in xx:
                    temp_sub_class = "LPG"    
                elif "PG(P-" in xx:
                    temp_sub_class = "PG(P-"
                elif "PG(O-" in xx:
                    temp_sub_class = "PG(O-"      
                else:
                    temp_sub_class = "PG"
                subclasses.append(temp_sub_class)

        elif jj == "PI":
            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""

                temp_sub_class = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                        if "LPI(O-" in i:
                            temp_sub_class =temp_sub_class +" | "+ "LPI(O-"
                        elif "LPI(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPI(P-"    
                        elif "LPI" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPI"    
                        elif "PI(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PI(P-"
                        elif "PI(O-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PI(O-"      
                        else:
                            temp_sub_class = temp_sub_class +" | "+"PI"

                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]

                        if "LPI(O-" in i:
                            temp_sub_class = "LPI(O-"
                        elif "LPI(P-" in i:
                            temp_sub_class = "LPI(P-"    
                        elif "LPI" in i:
                            temp_sub_class = "LPI"    
                        elif "PI(P-" in i:
                            temp_sub_class = "PI(P-"
                        elif "PI(O-" in i:
                            temp_sub_class = "PI(O-"      
                        else:
                            temp_sub_class = "PI"

                subclasses.append(temp_sub_class)
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append("NA")


            else:

                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append("NA")
                if "LPI(O-" in xx:
                    temp_sub_class = "LPI(O-"
                elif "LPI(P-" in xx:
                    temp_sub_class = "LPI(P-"    
                elif "LPI" in xx:
                    temp_sub_class = "LPI"    
                elif "PI(P-" in xx:
                    temp_sub_class = "PI(P-"
                elif "PI(O-" in xx:
                    temp_sub_class = "PI(O-"      
                else:
                    temp_sub_class = "PI"
                subclasses.append(temp_sub_class)


    #####
        elif jj == "PS":
            if " | " in xx:

                temp_lipid = xx.split("|")
                temp_string_length = ""
                temp_string_saturation = ""
                between_parents_temp = ""

                temp_sub_class = ""
                for i in temp_lipid:
                    if len(temp_string_length) >0:
                        temp_string_length = temp_string_length+" | "+ i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = temp_string_saturation+" | "+ i[i.find(":")+1:i.find(")")]
                        between_parents_temp = between_parents_temp +" | "+ i[i.find("("):i.find(")")+1]
                        if "LPS(O-" in i:
                            temp_sub_class =temp_sub_class +" | "+ "LPS(O-"
                        elif "LPS(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPS(P-"    
                        elif "LPS" in i:
                            temp_sub_class = temp_sub_class +" | "+"LPS"    
                        elif "PS(P-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PS(P-"
                        elif "PS(O-" in i:
                            temp_sub_class = temp_sub_class +" | "+"PS(O-"      
                        else:
                            temp_sub_class = temp_sub_class +" | "+"PS"

                    else:

                        temp_string_length = i[i.find("(")+1:i.find(":")]
                        temp_string_saturation = i[i.find(":")+1:i.find(")")]
                        between_parents_temp = i[i.find("("):i.find(")")+1]

                        if "LPS(O-" in i:
                            temp_sub_class = "LPS(O-"
                        elif "LPS(P-" in i:
                            temp_sub_class = "LPS(P-"    
                        elif "LPS" in i:
                            temp_sub_class = "LPS"    
                        elif "PS(P-" in i:
                            temp_sub_class = "PS(P-"
                        elif "PS(O-" in i:
                            temp_sub_class = "PS(O-"      
                        else:
                            temp_sub_class = "PS"

                subclasses.append(temp_sub_class)
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(between_parents_temp)
                FA_Chain.append("NA")


            else:

                temp_string_length = xx[xx.find("(")+1:xx.find(":")]
                temp_string_saturation = xx[xx.find(":")+1:xx.find(")")]
                chain_length.append(temp_string_length)
                saturation.append(temp_string_saturation)
                between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
                FA_Chain.append("NA")
                if "LPS(O-" in xx:
                    temp_sub_class = "LPS(O-"
                elif "LPS(P-" in xx:
                    temp_sub_class = "LPI(P-"    
                elif "LPS" in xx:
                    temp_sub_class = "LPS"    
                elif "PS(P-" in xx:
                    temp_sub_class = "PS(P-"
                elif "PS(O-" in xx:
                    temp_sub_class = "PS(O-"      
                else:
                    temp_sub_class = "PS"
                subclasses.append(temp_sub_class)



        elif jj == "SM":
            FA_Chain.append("None")

            subclasses.append("SM")
            between_parenthese.append(xx[xx.find("("):xx.find(")")+1])
            saturation.append("NA")

            chain_length.append("NA")

    merged_df["subclasses"] = subclasses
    merged_df["chain_length"] = subclasses
    merged_df["saturation"] = saturation
    merged_df["between_parenthese"] = between_parenthese
    merged_df["FA_Chain"] = FA_Chain

    return merged_df

def full_parse(data_base_name_location,mzml_folder, folder_name_to_save,labels_df, blank_name,file_name_to_save,tolerance, custom_data=False,remove_std = True,
               save_data=False):
    mrm_database = read_mrm_list(data_base_name_location,remove_std=remove_std,custom_data=custom_data)
    df = mzml_parser(mzml_folder)
    df_matched = match_lipids_parser(mrm_database,df, tolerance=tolerance)
    df_matched = add_labels(labels_df,df_matched)
    
    df_matched = subtract_blank(labels_df,df_matched,blank_name)
    if custom_data ==False:
        df_matched = add_subclass_and_length(df_matched)
    if save_data == True:
        
        save_dataframe(df_matched, folder_name_to_save, file_name_to_save)

    return df_matched

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



def prep_edge_R(merged_df,json_list_pairs,Pre_edge_r_path,blank_name,labels_list,extra_name):
    for i in range(len(json_list_pairs)):
        json1 = json_list_pairs[i][0]
        json2 = json_list_pairs[i][1]

        json_blank = {"Sample Name":[blank_name]}

        filtered_df1 = filter_dataframe(merged_df,json1)
        filtered_df2 = filter_dataframe(merged_df,json2)
        filtered_blank = filter_dataframe(merged_df,json_blank)



        # filtered_df1 = filtered_df1.groupby(labels_list)['Intensity'].mean().reset_index()
        # filtered_df2 = filtered_df2.groupby(labels_list)['Intensity'].mean().reset_index()
        # filtered_blank = filtered_blank.groupby(labels_list)['Intensity'].mean().reset_index()



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
        title = title + extra_name
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

