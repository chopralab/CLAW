import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt

















full_file_path = "List_MRMs_edit.csv"
df = pd.read_csv(full_file_path, delimiter="\t",header=0)

print(df)

my_list = list(df)

print(my_list)
print(df['Lipid Name'])


lipids_orig = df['Lipid Name']
transitions_orig = df['Transition']

transitions_set = [*set(transitions_orig)]

transitions = []

print(len(transitions_orig))
print(len(lipids_orig))
print(len(transitions_set))

new_trans = []

duplicates = []

print("")
print("")
print("")
print("")
print("")
print("")
indexes = []
for j,i in enumerate(transitions_orig):
    transitions.append(i)
    print(type(i))
    if i not in new_trans:
        new_trans.append(i)
    else:
        duplicates.append(i)
        indexes.append(j)
        # print(transitions_orig.count(i))


for i in range(len(duplicates)):
    if duplicates[i] in transitions:
        print(duplicates[i],"Count:",transitions.count(duplicates[i]))
    else:
        print(i)

# print(df[0])

for i in indexes:
    print(lipids_orig[i])