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




def read_mrm_list(filename):
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
    return mrm_list_offical

mrm_database = read_mrm_list('lipid_database/Lipid_Database.xlsx')
mrm_database.tail()



#Create for loop to load all mzml files from the data folder into the run object from pymzml reader function and store in pandas dataframe
#Create empty dictionary to store all the data
df_OzESI = pd.DataFrame(columns=['Lipid','Parent_Ion','Product_Ion','Intensity','Retention_Time','Transition','Class','Sample_ID'])
###
# OzESI_time = {}
###

data_folder = os.listdir('./data_mzml/liver_LD/') #Path to the mzml files
path_to_mzml_files = './data_mzml/liver_LD/'
#data_dict = {} #Empty dictionary to store all the data
df = pd.DataFrame(columns=['Lipid','Parent_Ion','Product_Ion','Intensity','Transition','Class','Sample_ID'])
#Create a similar for loop, except store all data in a single pandas dataframe
df_all = pd.DataFrame(columns=['Lipid','Parent_Ion','Product_Ion','Intensity','Transition','Class','Sample_ID']) #Create empty pandas dataframe to store the data
#df_all = pd.DataFrame(columns=['Q1','Q3','Intensity','Transition','Lipid','Class']) #Create empty pandas dataframe to store the data




##My edit
for file in data_folder:
        if file.endswith('.mzML'):
                print(file)
                run = pymzml.run.Reader(path_to_mzml_files+file, skip_chromatogram=False) #Load the mzml file into the run object
                print('Spectrum # = ',run.get_spectrum_count())
                print('Chromatogram # =',run.get_chromatogram_count())


                
                #create pandas dataframe to store the data with the columns Parent Ion, Product Ion, Intensity, Transition Lipid and Class
                #df_sample = pd.DataFrame(columns=['Parent_Ion','Product_Ion','Intensity','Transition','Lipid','Class']) #Create empty pandas dataframe to store the data
                #df_sample = pd.DataFrame(columns=['Q1','Q3','Intensity','Transition','Lipid','Class']) #Create empty pandas dataframe to store the data
                q1_mz = 0 #Create empty variables to store the Q1 and Q3 m/z values
                q3_mz = 0
                count = 0 #Create a counter to keep track of the number of transitions
                for spectrum in run:
                        
                        ###
                        # if isinstance(spectrum,pymzml.spec.Chromatogram):
                        #         for time, intensity in spectrum.peaks():
                        #                 print(time, intensity)
                        #                 OzESI_time[time] = intensity
                        #         # OzESI_time.append(time_list)
                        ###

                        for element in spectrum.ID.split(' '):
                                # print('element',element)
                                intensity_store = np.array([])
                                if 'Q1' in element:
                                        #print('Q1',element)
                                        q1 = element.split('=')
                                        #print('q1',q1[1])
                                        q1_mz= np.round((float(q1[1])),1)
                                        # print('q1',q1)
                                
                                if 'Q3' in element:
                                        #print('Q3',element)
                                        q3 = element.split('=')
                                        #print('q3',q3[1])
                                        q3_mz=np.round(float(q3[1]),1)
                                        # print('q3',q3)
                                        # df_sample.loc[count,'Q1'] = q1_mz
                                        # df_sample.loc[count,'Q3'] = q3_mz
                                        
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
        #append df_all to df_all2
        df = df.append(df_all, ignore_index=True)
# df.tail(5) 
print(len(df))