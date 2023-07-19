import pandas as pd




merged_df_adjP = pd.read_csv("merged_df_adjP.csv")
merged_df_P = pd.read_csv("merged_df_P.csv")
filtered_df = pd.read_csv("filtered_df_FDR.csv")


full_df = pd.read_csv("cl_e1_full.csv")
print(list(full_df))

filtered_df2 = full_df.loc[full_df['PValue'] < 0.1]

print(len(filtered_df2))

columns_to_keep2 = ['C1', 'C2', 'C3', 'C4', 'C5', 'R1', 'R2', 'R3', 'R4', 'R5']


filtered_df2 = filtered_df2[columns_to_keep2]
# exit()





import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
def make_pca_plot(df, df_name):
    # transpose dataframe so that each column becomes a data point
    df_transposed = df.transpose()

    # standardize the data (mean=0, std=1)
    standardized_data = StandardScaler().fit_transform(df_transposed)

    # PCA transformation
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(standardized_data)

    # create a dataframe with the principal components
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    # add labels based on column names
    pca_df['label'] = ['C' if name.startswith('C') else 'R' for name in df_transposed.index]

    # get percentage of variance explained by each PC
    explained_var = pca.explained_variance_ratio_

    # plot
    plt.figure(figsize=(10,10))
    for label, color in [('C', 'red'), ('R', 'blue')]:
        mask = pca_df['label'] == label
        plt.scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], c=color, edgecolors='black', label=label)
    
    plt.title(df_name)
    plt.xlabel(f'PC1 ({explained_var[0]*100:.2f}%)')  # include % variance explained in label
    plt.ylabel(f'PC2 ({explained_var[1]*100:.2f}%)')  # include % variance explained in label
    plt.legend()

    # save plot
    plt.savefig(df_name + '.png')


# merged_df_adjP = pd.read_csv("merged_df_adjP.csv")
# merged_df_P = pd.read_csv("merged_df_P.csv")
# filtered_df = pd.read_csv("filtered_df.csv")




make_pca_plot(merged_df_adjP, 'merged_df_adjP')
make_pca_plot(merged_df_P, 'merged_df_P')
make_pca_plot(filtered_df, 'filtered_df_original')
make_pca_plot(filtered_df2, 'filtered_df_original P Value')





