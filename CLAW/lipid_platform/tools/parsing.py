##Remove duplicates
import pandas as pd

def add_suffix(df):
    df['Count'] = df.groupby(['Lipid', 'Class', 'Sample Name']).cumcount().add(1)
    df['Lipid'] = df.apply(lambda row: row['Lipid'] + '_' + str(row['Count']) if row['Count'] > 1 else row['Lipid'], axis=1)
    df.drop(columns=['Count'], inplace=True)
    return df
