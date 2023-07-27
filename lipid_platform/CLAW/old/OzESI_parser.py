

def process_chromatogram(OzESI_time):
    # Create dataframe from OzESI_time dictionary
    OzESI_rt_df = pd.DataFrame(list(OzESI_time.items()), columns=['Retention_Time', 'intensity'])
    
    # Split intensity column into three columns intensity, Parent_Ion and Product_Ion
    OzESI_rt_df[['intensity','Parent_Ion','Product_Ion']] = pd.DataFrame(OzESI_rt_df['intensity'].tolist(), index=OzESI_rt_df.index)
    
    # Round retention Retention_Time to 1 decimal place
    OzESI_rt_df['Retention_Time'] = round(OzESI_rt_df['Retention_Time'], 2)
    
    # Create a column called Transition with the Parent_Ion and Product_Ion
    OzESI_rt_df['Transition'] = OzESI_rt_df['Parent_Ion'].astype(str) + ' -> ' + OzESI_rt_df['Product_Ion'].astype(str)
    
    ########### HARDCODED TO DROP RETENTION TIMES BELOW 7 SECONDS ############
    #drop Rention_Time below 10.5 seconds and above 15.5 seconds
    OzESI_rt_df = OzESI_rt_df[OzESI_rt_df['Retention_Time'] > 10.5]
    OzESI_rt_df = OzESI_rt_df[OzESI_rt_df['Retention_Time'] < 15.5]

    # Get the top 10 records for each 'Transition' based on 'intensity'
    OzESI_rt_df_top = OzESI_rt_df.groupby('Transition').apply(lambda x: x.nlargest(1, 'intensity')).reset_index(drop=True)

    print('OzESI 1 largest per transition: \n', OzESI_rt_df_top)

    # peaks, _ = find_peaks(OzESI_rt_df['intensity'], height=0.5e5,distance=1000)
    # plt.scatter(OzESI_rt_df['Retention_Time'] ,OzESI_rt_df['intensity'])
    # plt.scatter(OzESI_rt_df.iloc[peaks]['Retention_Time'], OzESI_rt_df.iloc[peaks]['intensity'], "x")
    # plt.ylabel('Intensity')
    # plt.xlabel('Retention Time')
    # plt.title('OzESI LC Chromatogram')
    # plt.show()
    
    return OzESI_rt_df_top


def add_rt_intensity(df, OzESI_rt_df_top):

    #### STRING OF RETENTION TIMES
    # Group by 'Transition' and get lists of retention times and intensities
    transitions_to_rt = OzESI_rt_df_top.groupby('Transition')['Retention_Time'].apply(list).to_dict()
    transitions_to_intensity = OzESI_rt_df_top.groupby('Transition')['intensity'].apply(list).to_dict()
    
    # Use the map function to add retention times and intensities to the dataframe as strings
    df['Retention_Time'] = df['Transition'].map(transitions_to_rt).apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)
    df['Intensity_OzESI'] = df['Transition'].map(transitions_to_intensity).apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)
    

    #MEAN OF RENETION TIMES
    # Convert the 'Retention_Time' column from string of lists to lists of floats
    df['Retention_Time'] = df['Retention_Time'].apply(lambda x: [float(i) for i in x.split(', ')] if isinstance(x, str) else x)

    # Calculate the mean of 'Retention_Time' for each transition
    df['Mean_Retention_Time'] = df['Retention_Time'].apply(np.mean).round(2)

    # Convert the 'Intensity_OzESI' column from string of lists to lists of floats
    df['Intensity_OzESI'] = df['Intensity_OzESI'].apply(lambda x: [float(i) for i in x.split(', ')] if isinstance(x, str) else x)

    # Calculate the mean of 'Intensity_OzESI' for each transition
    df['Mean_Intensity_OzESI'] = df['Intensity_OzESI'].apply(np.mean)

    return df

def create_aldehyde_ion_dataframe():
    # Create a pandas dataframe with columns for DB_Position and Aldehyde_Ion
    df_OzESI = pd.DataFrame(columns=['DB_Position','Aldehyde_Ion'])

    # Loop over the range of DB_Position values and calculate the corresponding Aldehyde_Ion values
    for i in range(3, 21):
        df_OzESI.loc[i,'DB_Position'] = i
        df_OzESI.loc[i,'Aldehyde_Ion'] = 26 + (14 * (i-3)) 

    # Print the first 25 rows of the dataframe
    # print(df_OzESI.head(25))

    # Return the dataframe
    return df_OzESI

#OzESI_list = [3,5,7,9,11]
OzESI_list = [7,9,12]
#def calculate_n_minus_values(df_matched, df_OzESI, OzESI_list=[3,5,7,9,11], starting_column=9, last_column=14):
def calculate_n_minus_values(df_matched, df_OzESI, OzESI_list=[7,9,12], starting_column=9, last_column=14):
    """
    Given a pandas dataframe df and a dataframe df_OzESI containing DB_Position and Aldehyde_Ion values,
    calculates the n-i values for each i in OzESI_list by subtracting the corresponding Aldehyde_Ion value
    from the Parent_Ion column in df and storing the result in new columns named 'n-i' in df.
    The starting_column and last_column parameters specify the range of columns in which the n-i values should be stored.
    """
    # Create new columns in df for n-i values
    for i in OzESI_list:
        df_matched[f"n-{i}"] = df_matched["Parent_Ion"] - df_OzESI.loc[df_OzESI["DB_Position"] == i, "Aldehyde_Ion"].values[0]
    
    return df_matched



#Add Lipid and Labels columns
def add_lipid_info(df, OzESI_list):
    df_test = df.copy()
    
    for i in OzESI_list:
        df_test['n-' + str(i)] = df_test['n-' + str(i)].astype(float)
    
    for i in range(len(df_test)):
        if pd.isna(df_test.loc[i, 'Lipid']):
            parent_ion = df_test.loc[i, 'Parent_Ion']
            
            for j in range(len(df_test)):
                if parent_ion == df_test.loc[j, 'n-7'] and isinstance(df_test.loc[j, 'Lipid'], str):
                    #df_test.loc[i, 'Lipid'] = 'n-7 ' + df_test.loc[j, 'Lipid']
                    df_test.loc[i, 'Lipid'] =  df_test.loc[j, 'Lipid']
                    df_test.loc[i, 'Labels'] = 'n-7' + df_test.loc[j, 'Labels']
                elif parent_ion == df_test.loc[j, 'n-9'] and isinstance(df_test.loc[j, 'Lipid'], str):
                    #df_test.loc[i, 'Lipid'] = 'n-9 ' + df_test.loc[j, 'Lipid']
                    df_test.loc[i, 'Lipid'] =  df_test.loc[j, 'Lipid']
                    df_test.loc[i, 'Labels'] = 'n-9' + df_test.loc[j, 'Labels']
                elif parent_ion == df_test.loc[j, 'n-12'] and isinstance(df_test.loc[j, 'Lipid'], str):
                    #df_test.loc[i, 'Lipid'] = 'n-12 ' + df_test.loc[j, 'Lipid']
                    df_test.loc[i, 'Lipid'] = df_test.loc[j, 'Lipid']
                    df_test.loc[i, 'Labels'] = 'n-12' + df_test.loc[j, 'Labels']
    
    df_test.dropna(subset=['Lipid'], inplace=True)
    return df_test

df_test = add_lipid_info(df_OzESI_processed, OzESI_list)
df_test.head(None)
