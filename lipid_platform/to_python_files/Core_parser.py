
###All functions

#Function to read in MRM database
#Option to remove STDs from database##Not finished need option to use another database with no qualitative ACs
import pandas as pd

def read_mrm_list(filename,remove_std = True):
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
    return mrm_list_offical

#OzESI
OzESI_time = {}
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
                        # if isinstance(spectrum, pymzml.spec.Chromatogram):
                        #     for time, intensity in spectrum.peaks():
                        #         OzESI_time[time] = np.round(intensity), q1_mz, q3_mz

         

                        for element in spectrum.ID.split(' '):
                                intensity_store = np.array([])
                                if 'Q1' in element:
                                        q1 = element.split('=')
                                        q1_mz= np.round((float(q1[1])),1)

                                if 'Q3' in element:
                            
                                        q3 = element.split('=')

                                        q3_mz=np.round(float(q3[1]),1)
                                    #### OzESI (Key)time, ( Values) intensity and transitions stored in OzESI_time dictionary
                                        for time, intensity in spectrum.peaks():
                                            OzESI_time[time] = np.round(intensity), q1_mz, q3_mz

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


    # df_matched = df_matched.dropna()
    
    return df_matched


def save_dataframe(df, folder_name, file_name, max_attempts=5):
    folder_path = f'data_results/data/data_matching/{folder_name}'
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

def full_parse(data_base_name_location,mzml_folder, folder_name_to_save, file_name_to_save,tolerance,remove_std = True,
               save_data=False):
    mrm_database = read_mrm_list(data_base_name_location,remove_std=remove_std)
    df = mzml_parser(mzml_folder)
    df_matched = match_lipids_parser(mrm_database,df, tolerance=tolerance)
    
    if save_data == True:
        
        save_dataframe(df_matched, folder_name_to_save, file_name_to_save)

    return df_matched