import pandas as pd

# Load the CSV file into a DataFrame
pooja_df = pd.read_csv('Pooja.csv')

# Assuming 'Method Path' is the second column, we filter the rows
# We can access the second column by its index which is 1 (since Python is zero-indexed)
filtered_pooja_df = pooja_df[pooja_df.iloc[:, 1].str.contains('FA_AMP_MRMs_FlowInj.m', na=False)]

# Save the filtered DataFrame to a new CSV file
filtered_pooja_df.to_csv('Filtered_Pooja.csv', index=False)
