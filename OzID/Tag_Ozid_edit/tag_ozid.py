import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt



df = pd.read_excel("FA MUFA 18_1.xlsx")


headings = list(df)

print(headings)

headings2 = headings[5:]

print(list(df))

print(df['loss 18:1'])



print(headings)
print(headings2)

# print(df[0])


#####Auto Generate
#14.01565 n2 - n17
#N1 - +1.9793

#N2 - -12.0264

transitions1 = df["m/z + NH4"]
transitions2 = df["loss 18:1"]
lipid_name = df["TG "]


mass_subtract = [1.9793,-12.0264,-26.052,-40.0677,-54.0833,-68.099,-82.1146,-96.1303,-110.1459,-124.1616,-138.1772,-152.1929,-166.2085,-180.2242,-194.2398,-208.2555,-222.2711]
db_position = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]

for i in range(2,17):
    number = i*-14.01565
    print(number)

for i in range(len(db_position)):
    print(db_position[i],mass_subtract[i])


print(len(transitions1))
print(len(transitions2))
print(len(lipid_name))

no_ozone_csv = []
ozone_csv = []

import csv

# name of csv file 
filename = "new_transitions.csv"
    
# writing to csv file 
csvfile= open(filename, 'w')
csvwriter = csv.writer(csvfile,delimiter="\t") 



for i in range(len(transitions1)):
    # to_add_no_ozone = str(lipid_name[i])+"\t"+str(transitions1[i])+"->"+str(transitions2[i])
    to_add_no_ozone = [str(lipid_name[i]),str(round(transitions1[i],1))+"->"+str(round(transitions2[i],1))]
    no_ozone_csv.append(to_add_no_ozone)
    print(to_add_no_ozone)
    csvwriter.writerow(to_add_no_ozone)
    for j in range(len(mass_subtract)):
        number = round(transitions1[i]-mass_subtract[j],1)
        try:
            number = str(number).replace('.0', '')
        except:
            pass
            # print("passed")
        # to_add_ozone = str(lipid_name[i])+"_OZONE_DB_Position_N"+str(db_position[j])+"\t"+str(transitions1[i])+"->"+str(number)
        to_add_ozone = [str(lipid_name[i])+"_OZONE_DB_1_Position_N"+str(db_position[j]),str(round(transitions1[i],1))+"->"+str(number)]
        print(to_add_ozone)
        ozone_csv.append(to_add_ozone)
        csvwriter.writerow(to_add_ozone)


        










