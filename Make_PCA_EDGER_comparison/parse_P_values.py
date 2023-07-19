import pandas as pd
import re
# Load the data
df = pd.read_excel('P_adj_P for.xlsx')

full_df = pd.read_csv("cl_e1_full.csv")

full_df['lipid'] = full_df['lipid'].str.replace(' ', '')

filtered_df = full_df.loc[full_df['FDR'] < 0.1]

print(list(full_df))

columns_to_keep = ["lipid",'C1', 'C2', 'C3', 'C4', 'C5', 'R1', 'R2', 'R3', 'R4', 'R5']
columns_to_keep2 = ['C1', 'C2', 'C3', 'C4', 'C5', 'R1', 'R2', 'R3', 'R4', 'R5']


subset_df = full_df[columns_to_keep]

print(list(subset_df))

# exit()

# Extract the 'Pvalue' and 'AdjustedP' columns
Pvalue_list = df['Pvalue'].dropna().tolist()
AdjustedP_list = df['AdjustedP'].dropna().tolist()

# Print the lists
print('Pvalue list:', Pvalue_list)
print('AdjustedP list:', AdjustedP_list)

feature_list_P = []
feature_list_adjP = []

for item in AdjustedP_list:
    # Extract the feature using regex
    feature = re.search('feature : (.*) <br> log2', item)
    
    if feature is not None: 
        feature = feature.group(1)
        # Remove spaces from the feature
        feature = feature.replace(" ", "")
        feature_list_adjP.append(feature)




for item in Pvalue_list:
    # Extract the feature using regex
    feature = re.search('feature : (.*) <br> log2', item)
    
    if feature is not None: 
        feature = feature.group(1)
        # Remove spaces from the feature
        feature = feature.replace(" ", "")
        feature_list_P.append(feature)

print(feature_list_P)
print(feature_list_adjP)



# Convert list to DataFrame
df_feature_list_P = pd.DataFrame(feature_list_P, columns=['lipid'])

# Merge with subset_df
merged_df_P = pd.merge(df_feature_list_P, subset_df, how='inner', on='lipid')

# Convert list to DataFrame
df_feature_list_adjP = pd.DataFrame(feature_list_adjP, columns=['lipid'])

# Merge with subset_df
merged_df_adjP = pd.merge(df_feature_list_adjP, subset_df, how='inner', on='lipid')


merged_df_adjP = merged_df_adjP[columns_to_keep2]
merged_df_P = merged_df_P[columns_to_keep2]
filtered_df = filtered_df[columns_to_keep2]



merged_df_adjP.to_csv("merged_df_adjP.csv",index=False)
merged_df_P.to_csv("merged_df_P.csv",index=False)
filtered_df.to_csv("filtered_df_FDR.csv",index=False)



print(len(merged_df_adjP),"merged_df_adjP")
print(len(merged_df_P),"merged_df_P")
print(len(feature_list_P),"feature_list_P")
print(len(feature_list_adjP),"feature_list_adjP")
print(len(full_df),"full_df")
