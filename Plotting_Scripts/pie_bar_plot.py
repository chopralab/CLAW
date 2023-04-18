import os
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt

# logFC
###Chooses FDR or P value for making bar plot
threshold = 'FDR'
# threshold = 'logFC'
# threshold_value = 1

# threshold = 'PValue'
threshold_value = 0.1
ncomponets = 4

##Create Directory for 
lst=os.listdir('data_files_name_changed')

for ij in lst:
    if ".csv" in ij:
        ##Directory and file name for the count of up and down regulated lipids
        new_xcel_file_name_for_palak = "sig_lipids_count/"+ij[:-4]+'_'+threshold+'_'+str(threshold_value)+'_.xlsx'

        print(ij[:-4])
        # exit()
        ##Gets path to open dataframe
        full_file_path = 'data_files_name_changed/'+ij

        ###Opens the file and makes it a dataframe and adds the 
        df = pd.read_csv(full_file_path)




        ##Get list of names for finding blank
        my_list = list(df)
        print(my_list)



        ###Finds the blank in order to use it for normilization also is the la
        for j,i in enumerate(my_list):
            if 'Blank' in i:
                end_of_list=j
        
        print(end_of_list)



        # exit()
        files_list = []



        ##Intensities start at position 8
        for i in range(8,end_of_list+1):
            files_list.append(my_list[i])






        # exit()

        # print(my_list[8])
        # print(my_list[9])
        # print(my_list[10])
        # print(my_list[11])
        # print(my_list[12])
        # # exit()
        bar_plot = ij[:-4]+' _Bar_Plot_'+threshold+'_'+str(threshold_value)+'_'

        pi_title = ij[:-4]+' Pie Chart_'+threshold+'_'+str(threshold_value)+'_'
        pi_Log2_title = ij[:-4]+' Pie Chart Log2_'+threshold+'_'+str(threshold_value)+'_'
        PCA_title = ij[:-4]+" PCA_minus_blank_Ncomponets_"+str(ncomponets)

        blank_int = np.array(df[my_list[end_of_list]])

        print(blank_int)
        # exit()
        # files_list = [my_list[8],my_list[9],my_list[10],my_list[11]]

        itense_values = []
        print(files_list[:-1])
        # exit()
        for i in files_list[:-1]:
            value = np.array(df[i]) - blank_int
            print(np.array(value))
            value = np.array(value) # /max(value)
            print(type(value))
            # value = np.log(value)
            itense_values.append(value)



        thresh_hold_values = []

        for i in range(len(df)):
            x = df[threshold][i]

            if x < threshold_value:
                thresh_hold_values.append(x)
        # exit()


        new_intensity_values_sig = []

        for i in range(len(itense_values)):
            new_intensity_values_sig.append(itense_values[i][:len(thresh_hold_values)])

        print(len(thresh_hold_values))

        print(len(new_intensity_values_sig))
        print(len(itense_values))
        print(len(new_intensity_values_sig[i]))
        print(len(itense_values[i]))
        # continue
        # exit()

        print(df)
        print(len(df))

        lipid_log2FC = []
        # list_of_sig_lipids = []
        list_of_sig_lipids = []
        up_regulated = []
        up_regulated_values = []
        down_regulated = []
        down_regulated_values = []
        
        print(math.log2(8))



        # exit()
        for i in range(len(df)):
            x = df[threshold][i]
            # print(x)
            if x < threshold_value: ##For FDR or P value
            # print(abs(x))
            # if abs(x) > abs(threshold_value): ##For LogFC
                if (df['type'][i]) == 'PCandSM':
                    xx =(df['lipid'][i])
                    print(xx)
                    if 'PC' in xx:

                        list_of_sig_lipids.append("PC")
                        lipid_log2FC.append(df['logFC'][i])

                        if df['logFC'][i] >0:
                            up_regulated.append("PC")
                            up_regulated_values.append(df['logFC'][i])

                        elif df['logFC'][i] <0:
                            down_regulated.append("PC")
                            down_regulated_values.append(df['logFC'][i])

                    else: #'SM' in xx:
                        lipid_log2FC.append(df['logFC'][i])

                        list_of_sig_lipids.append("SM")

                        if df['logFC'][i] >0:
                            up_regulated.append("SM")
                            up_regulated_values.append(df['logFC'][i])

                        elif df['logFC'][i] <0:
                            down_regulated.append("SM")
                            down_regulated_values.append(df['logFC'][i])

                elif "TAG" in (df['type'][i]):
                    list_of_sig_lipids.append("TAG")
                    lipid_log2FC.append(df['logFC'][i])

                    if df['logFC'][i] >0:
                        up_regulated.append("TAG")
                        up_regulated_values.append(df['logFC'][i])

                    elif df['logFC'][i] <0:
                        down_regulated.append("TAG")
                        down_regulated_values.append(df['logFC'][i])

                else:
                    lipid_log2FC.append(df['logFC'][i])

                    list_of_sig_lipids.append(df['type'][i])

                    if df['logFC'][i] >0:
                        up_regulated.append(df['type'][i])
                        up_regulated_values.append(df['logFC'][i])
                    elif df['logFC'][i] <0:
                        down_regulated.append(df['type'][i])

                        down_regulated_values.append(df['logFC'][i])

        lipid_classes = []
        counts = []
        for i in list_of_sig_lipids:
            # print(i)
            if i not in lipid_classes:
                lipid_classes.append(i)


        print(len(lipid_classes))
        print(len(list_of_sig_lipids))

        lipid_classes.sort()

        counts_up2 = []
        counts_down2 = []

        counts_up = []
        counts_down = []

        for i in lipid_classes:
            counts.append(list_of_sig_lipids.count(i))
            counts_up.append(up_regulated.count(i))
            counts_down.append(down_regulated.count(i))

        print(lipid_classes)
        print(list_of_sig_lipids.count("TAGs1"))
        print(list_of_sig_lipids.count("TAGs2"))
        
        len_tags1 = []
        len_tags2 = []
        for i in list_of_sig_lipids:
            if 'TAGs1' == i:
                print(i)
                len_tags1.append('TAGs1')
            if 'TAGs2' == i:
                len_tags2.append('TAGs2')


        print(len(len_tags1))
        print(len(len_tags2))
        # exit()


        sums_of_log2s = []

        ###Get sums of log2 values
        for i in lipid_classes:
            temp = 0
            for j in range (len(lipid_log2FC)):
                if list_of_sig_lipids[j] == i:
                    print(list_of_sig_lipids[j])
                    temp = temp + abs(lipid_log2FC[j])
            sums_of_log2s.append(temp)






        data_dict= {"Lipid Classes":lipid_classes,"Signifigantly changed Lipids":counts,"Up Regulated Lipids":counts_up,"Down Regulated Lipids":counts_down}

        print(lipid_classes)
        print(counts)
        print(counts_down)
        print(counts_up)
        print(up_regulated)
        print(down_regulated)

        up_regulated2 = []
        down_regulated2 = []

        for i in down_regulated:
            if i not in down_regulated2:
                down_regulated2.append(i)


        for i in up_regulated:
            if i not in up_regulated2:
                up_regulated2.append(i)





        up_regulated_values2 = []

        down_regulated_values2 = []


    ###Gets the sums of upregulated values
        for i in up_regulated2:
            temp = 0
            counts_up2.append(up_regulated.count(i))

            for j in range(len(up_regulated)):
                up_regulated.count(i)

                if up_regulated[j] ==i:
                    temp = temp +up_regulated_values[j]
            up_regulated_values2.append(abs(temp))


        for i in down_regulated2:
            temp = 0
            counts_down2.append(down_regulated.count(i))

            for j in range(len(down_regulated)):
                if down_regulated[j] ==i:

                    temp = temp +down_regulated_values[j]
            down_regulated_values2.append(abs(temp))

        print(ij)

        print(len(down_regulated_values2))
        print(len(counts_down))
        print((counts_down))
        print(len(down_regulated2))
        print(len(down_regulated))
        print(len(down_regulated))


        true_down_values = []
        true_up_values = []

        for i in lipid_classes:
            temp_up = []
            temp_down = [] 
            for j in range(len(down_regulated2)):
                if i == down_regulated2[j]:
                    temp_down.append(down_regulated_values2[j])

            for k in range(len(up_regulated2)):
                if i == up_regulated2[k]:
                    temp_up.append(up_regulated_values2[k])

            true_down_values.append(sum(temp_down))
            true_up_values.append(sum(temp_up))

        print(np.array(true_down_values)/np.array(counts_down))

        

        normalized_down_values = np.array(true_down_values)/np.array(counts_down)

        normalized_up_values = np.array(true_up_values)/np.array(counts_up)


        normalized_all_values = np.array(sums_of_log2s)/np.array(counts)

        normalized_all_values = np.nan_to_num(normalized_all_values, copy=True, nan=0.0, posinf=None, neginf=None)
        normalized_up_values = np.nan_to_num(normalized_up_values, copy=True, nan=0.0, posinf=None, neginf=None)
        normalized_down_values = np.nan_to_num(normalized_down_values, copy=True, nan=0.0, posinf=None, neginf=None)

        df2 = pd.DataFrame(data_dict)

        df2.to_excel(new_xcel_file_name_for_palak,index=None)



        print(counts)
        print(lipid_classes)
        print(sums_of_log2s)


        if sum(counts) != 0:
            plt.bar(lipid_classes, counts, color ='maroon',
                    width = 0.4)
            
            plt.xlabel("Class of Signifigant Lipid")
            plt.ylabel("Number of Lipids")
            plt.title("Signifigant Lipids")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("sig_lipids.png")
            plt.tight_layout()

            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)

            plt.close()
            plt.cla()
            plt.clf()

    ###UP_down_sig_all_in_one
        if sum(counts) != 0:

            barWidth = 0.25
            fig = plt.subplots(figsize =(12, 8))



            # Set position of bar on X axis
            br1 = np.arange(len(lipid_classes))
            br2 = [x + barWidth for x in br1]
            # br3 = [x + barWidth for x in br2]

            # Make the plot
            plt.bar(br1, counts_down, color ='black', width = barWidth,
                    edgecolor ='grey', label ='Down Regulated')
            plt.bar(br2, counts_up, color ='gray', width = barWidth,
                    edgecolor ='grey', label ='Up Regulated')
            # plt.bar(br3, counts_down, color ='b', width = barWidth,
            #         edgecolor ='grey', label ='Down Regulated')

            # Adding Xticks
            plt.xlabel('Lipid Class', fontweight ='bold', fontsize = 15)
            plt.ylabel('Number of Signifigant Lipids', fontweight ='bold', fontsize = 15)
            plt.xticks([r + barWidth for r in range(len(lipid_classes))],lipid_classes)
            plt.title(bar_plot.replace("_"," "))
            plt.tight_layout()


            plt.legend()
            plt.show()
            save_title_bar_up = "combined_bar/"+bar_plot +str("counts_combined_bar.png")

            plt.savefig(save_title_bar_up ,dpi=700)


            plt.close()
            plt.cla()
            plt.clf()

            ##LOG_FC
            fig = plt.subplots(figsize =(12, 8))



            # Set position of bar on X axis
            br1 = np.arange(len(lipid_classes))
            br2 = [x + barWidth for x in br1]
            # br3 = [x + barWidth for x in br2]

            # Make the plot
            plt.bar(br1, true_down_values, color ='black', width = barWidth,
                    edgecolor ='grey', label ='Down Regulated')
            plt.bar(br2, true_up_values, color ='gray', width = barWidth,
                    edgecolor ='grey', label ='Up Regulated')
            # plt.bar(br3, counts_down, color ='b', width = barWidth,
            #         edgecolor ='grey', label ='Down Regulated')

            # Adding Xticks
            plt.xlabel('Lipid Class', fontweight ='bold', fontsize = 15)
            plt.ylabel('Sum of LogFC of Signifigant Lipids', fontweight ='bold', fontsize = 15)
            plt.xticks([r + barWidth for r in range(len(lipid_classes))],lipid_classes)
            plt.title(bar_plot.replace("_"," "))
            plt.tight_layout()


            plt.legend()
            plt.show()
            save_title_bar_up = "combined_bar/"+bar_plot +str("LOG_FC_combined_bar.png")

            plt.savefig(save_title_bar_up ,dpi=700)


            plt.close()
            plt.cla()
            plt.clf()


            ##LOG_FC
            fig = plt.subplots(figsize =(12, 8))



            # Set position of bar on X axis
            br1 = np.arange(len(lipid_classes))
            br2 = [x + barWidth for x in br1]
            # br3 = [x + barWidth for x in br2]

            # Make the plot
            plt.bar(br1, normalized_down_values, color ='black', width = barWidth,
                    edgecolor ='grey', label ='Down Regulated')
            plt.bar(br2, normalized_up_values, color ='gray', width = barWidth,
                    edgecolor ='grey', label ='Up Regulated')
            # plt.bar(br3, counts_down, color ='b', width = barWidth,
            #         edgecolor ='grey', label ='Down Regulated')

            # Adding Xticks
            plt.xlabel('Lipid Class', fontweight ='bold', fontsize = 15)
            plt.ylabel('Average LogFC of Signifigant Lipids', fontweight ='bold', fontsize = 15)
            plt.title(bar_plot.replace("_"," "))
            plt.xticks([r + barWidth for r in range(len(lipid_classes))],lipid_classes)
            plt.tight_layout()

            plt.legend()
            plt.show()
            save_title_bar_up = "combined_bar/"+bar_plot +str("Normalized_LOG_FC_combined_bar.png")

            plt.savefig(save_title_bar_up ,dpi=700)


            plt.close()
            plt.cla()
            plt.clf()


        if sum(counts_up) != 0:
        
        # creating the bar plot
            print(up_regulated2)
            print(counts_up)
            plt.bar(lipid_classes, counts_up, color ='maroon',
                    width = 0.4)
            plt.tight_layout()
            
            plt.xlabel("Class of Signifigant Lipid")
            plt.ylabel("Number of Up Regulated Lipids")
            plt.title("Upregulated Lipids")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("_up_regulated.png")

            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)

            plt.close()
            plt.cla()
            plt.clf()



            plt.bar(up_regulated2, up_regulated_values2, color ='maroon',
                    width = 0.4)
            
            plt.xlabel("Class of Signifigant Lipid")
            plt.ylabel("LogFC of Up Regulated Lipids")
            plt.title("Upregulated Lipids Log FC")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("_LOGFC_up_regulated.png")
            plt.tight_layout()

            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)

            plt.close()
            plt.cla()
            plt.clf()





            plt.bar(lipid_classes, normalized_up_values, color ='maroon',
                    width = 0.4)
            
            plt.xlabel("Class of Signifigant Lipid")
            plt.ylabel("Normalized LogFC of Up Regulated Lipids")
            plt.title("Upregulated Lipids Normalized LOG FC")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("_LOGFC_NORMALIZED_up_regulated.png")

            plt.tight_layout()

            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)
            plt.close()
            plt.cla()
            plt.clf()



        if sum(counts_down) != 0:
            print(down_regulated2)
            print(counts_down)
            plt.bar(lipid_classes, counts_down, color ='maroon',
                    width = 0.4)
            
            plt.xlabel("Class of Down regulated Lipid")
            plt.ylabel("Number of Lipids  Down regulated")
            plt.title("Down regulated Lipids")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("_down_regulated.png")
            plt.tight_layout()

            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)

            plt.close()
            plt.cla()
            plt.clf()






            plt.bar(down_regulated2, down_regulated_values2, color ='maroon',
                    width = 0.4)
            
            plt.xlabel("Class of Signifigant Lipid")
            plt.ylabel("LogFC of Down Regulated Lipids")
            plt.title("Down Regulated Lipids")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("_NORMALIZED_down_regulated.png")

            plt.tight_layout()
            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)

            plt.close()
            plt.cla()
            plt.clf()



            plt.bar(lipid_classes, normalized_down_values, color ='maroon',
                    width = 0.4)
            
            plt.xlabel("Class of Signifigant Lipid")
            plt.ylabel("LogFC Normalized of Down Regulated Lipids")
            plt.title("Down Regulated Lipids")
            plt.show()

            save_title_bar_up = "plots/bar/"+bar_plot +str("_LOGFC_Normalized_down_regulated.png")
            plt.tight_layout()

            # show plot
            plt.show()
            plt.savefig(save_title_bar_up ,dpi=700)

            plt.close()
            plt.cla()
            plt.clf()









        from matplotlib import pyplot as plt
        import numpy as np
        # Creating plot
    ###PCA uncommented


        # if len(thresh_hold_values) > ncomponets:

        from sklearn.decomposition import PCA


        pca = PCA(n_components = ncomponets)
        components = pca.fit_transform(new_intensity_values_sig)
        print(components)

        print(components[0])
        print(components[0][0])


        X = []
        Y = []



        for i in range(len(components)):
            X.append(components[i][0])
            Y.append(components[i][1])

        print(pca.explained_variance_ratio_)

        plt.scatter(X,Y)
        xlabel = "PC1_"+ str(100*pca.explained_variance_ratio_[0])[:2]+"%"
        ylabel = "PC2_"+ str(100*pca.explained_variance_ratio_[1])[:2]+"%"
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        for i, label in enumerate(files_list[:-1]):
            plt.annotate(label, (X[i], Y[i]))




        plt.title(PCA_title)
        # show plot
        plt.show()
        save_title_PCA = "plots/"+PCA_title +str("_PCA.png")
        plt.savefig(save_title_PCA)

        plt.close()
        plt.cla()
        plt.clf()



    ####Where Pie charts are made
        ##Old Colors
        # colors = []
        # colors.append(["ACs","#4e08ff"])
        # colors.append(["AC","#4e08ff"])
        # colors.append(["CEs","#e82a60"])
        # colors.append(["CE","#e82a60"])
        # colors.append(["Cer","#4fff00"])
        # colors.append(["CER","#4fff00"])
        # colors.append(["Ceramide","#4fff00"])
        # colors.append(["Ceramides","#4fff00"])
        # colors.append(["FFA","red"])
        # colors.append(["TAG2","#1e90ff"])
        # colors.append(["TAG","#1e90ff"])
        # colors.append(["TAGs","#1e90ff"])
        # colors.append(["TAG1","#1e90ff"])
        # colors.append(["TAG","#1e90ff"])
        # colors.append(["PSs","#808080"])
        # colors.append(["PS","#808080"])
        # colors.append(["PC","#ffff00"])
        # colors.append(["PEs","#228b22"])
        # colors.append(["PE","#228b22"])
        # colors.append(["PIs","#8fbc8f"])
        # colors.append(["PI","#8fbc8f"])
        # colors.append(["SM","#ff8c00"])
        # colors.append(["PGs","#8a1c1c"])
        # colors.append(["PGs","#8a1c1c"])

###new colors
        # colors = []
        # # colors.append(["ACs","#E49B0F"])
        # colors.append(["AC","#E49B0F"])
        # # colors.append(["CEs","#000000"])
        # colors.append(["CE","#000000"])
        # colors.append(["FFA","#4596a7"])
        # # colors.append(["TAG2","#de2e95"])
        # colors.append(["TAG","#de2e95"])
        # # colors.append(["TAGs","#de2e95"])
        # # colors.append(["TAG1","#de2e95"])
        # # colors.append(["TAG","#de2e95"])
        # colors.append(["SM","red"])
        # # colors.append(["Sm","red"])
        # # colors.append(["sm","red"])
        # colors.append(["PC","#6e34a4"])
        # # colors.append(["PEs","#c3a2e2"])
        # colors.append(["PE","#c3a2e2"])
        # colors.append(["PGs","#9fc5e8"])
        # colors.append(["PGs","#9fc5e8"])
        # colors.append(["PSs","#e1a7c3"])
        # colors.append(["PS","#e1a7c3"])


        # ###My own colors
        # # colors.append(["Cer","#b6d7a8"])
        # colors.append(["CER","#b6d7a8"])
        # # colors.append(["Ceramide","#b6d7a8"])
        # # colors.append(["Ceramides","#b6d7a8"])
        # # colors.append(["PIs","blue"])
        # colors.append(["PI","blue"])
        # # colors.append(["Pi","blue"])




        # colors_2_use = []
        # print(lipid_classes)
        # print("done")

        # for i in range(len(lipid_classes)):
        #     for j in range(len(colors)):
        #         if lipid_classes[i] == colors[j][0]:
        #             colors_2_use.append(colors[j][1])
        

        # fig = plt.figure(figsize =(10, 7))
        # plt.pie(counts, labels = lipid_classes,autopct='%.1f%%',colors=colors_2_use )
        # plt.title(pi_title.replace("_"," ")+" Total "+str(sum(counts)))
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        # save_title_pi = "plots/pie/"+pi_title +str("_PIE_.png")
        # # save_title_pi = save_title_pi.replace("_"," ")
        # plt.tight_layout()

        # # show plot
        # plt.show()
        # plt.savefig(save_title_pi ,dpi=700)

        # plt.close()
        # plt.cla()
        # plt.clf()




    ####Where Pie charts are made
        
        # colors = []
        # colors.append(["ACs","#4e08ff"])
        # colors.append(["CEs","#e82a60"])
        # colors.append(["Cer","#4fff00"])
        # colors.append(["FFA","red"])
        # colors.append(["FFAs","red"])
        # colors.append(["TAG2","#1e90ff"])
        # colors.append(["TAG1","#1e90ff"])
        # colors.append(["TAG","#1e90ff"])
        # colors.append(["TAGs","#1e90ff"])
        # colors.append(["PSs","#808080"])
        # colors.append(["PC","#ffff00"])
        # colors.append(["PCs","#ffff00"])
        # colors.append(["PEs","#228b22"])
        # colors.append(["PIs","#8fbc8f"])
        # colors.append(["SM","#ff8c00"])
        # colors.append(["SMs","#ff8c00"])
        # colors.append(["PGs","#8a1c1c"])



###new colors
        colors = []
        # colors.append(["ACs","#E49B0F"])
        colors.append(["AC","#E49B0F"])
        # colors.append(["CEs","#000000"])
        colors.append(["CE","#666699"])
        colors.append(["FFA","#4596a7"])
        # colors.append(["TAG2","#de2e95"])
        colors.append(["TAG","#cc40cc"])
        # colors.append(["TAGs","#de2e95"])
        # colors.append(["TAG1","#de2e95"])
        # colors.append(["TAG","#de2e95"])
        colors.append(["SM","red"])
        # colors.append(["Sm","red"])
        # colors.append(["sm","red"])
        colors.append(["PC","#6e34a4"])
        # colors.append(["PEs","#c3a2e2"])
        colors.append(["PE","#c3a2e2"])
        colors.append(["PG","#9fc5e8"])
        colors.append(["PGs","#9fc5e8"])
        colors.append(["PSs","#e1a7c3"])
        colors.append(["PS","#e1a7c3"])


        ###My own colors
        # colors.append(["Cer","#b6d7a8"])
        colors.append(["CER","#b6d7a8"])
        # colors.append(["Ceramide","#b6d7a8"])
        # colors.append(["Ceramides","#b6d7a8"])
        # colors.append(["PIs","blue"])
        colors.append(["PI","#cabcab"])
        # colors.append(["Pi","blue"])


        colors_2_use = []

        print("done")

        # lipid_classes.sort()

        for i in range(len(lipid_classes)):
            for j in range(len(colors)):
                if lipid_classes[i] == colors[j][0]:
                    colors_2_use.append(colors[j][1])
                    print(colors[j][0],"  ",lipid_classes[i],"    ",colors[j][1])
        print(lipid_classes)
        print(counts)     
        print(ij)
        # exit()
        print(lipid_classes)
        print(colors_2_use)
        fig = plt.figure(figsize =(10, 7))
        plt.pie(counts, labels = lipid_classes,autopct='%.1f%%',colors=colors_2_use )
        plt.title(pi_title.replace("_"," ")+" Total "+str(sum(counts)))
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        save_title_pi = "plots/pie/"+pi_title +str("_ALL_PIE_.png")
        # save_title_pi = save_title_pi.replace("_"," ")

        # show plot
        plt.show()
        plt.tight_layout()

        plt.savefig(save_title_pi, dpi=600)

        plt.close()
        plt.cla()
        plt.clf()


        # exit()


        colors_2_use2 = []

        print("done")

        # lipid_classes.sort()

        for i in range(len(up_regulated2 )):
            for j in range(len(colors)):
                if up_regulated2 [i] == colors[j][0]:
                    colors_2_use2.append(colors[j][1])
                    print(colors[j][0],"  ",up_regulated2 [i],"    ",colors[j][1])
        print(lipid_classes)
        print(counts)     
        print(ij)
        # exit()
        print(lipid_classes)
        print(colors_2_use)
        fig = plt.figure(figsize =(10, 7))
        plt.pie(counts_up2  , labels = up_regulated2 ,autopct='%.1f%%',colors=colors_2_use2 )
        plt.title(pi_title.replace("_"," ")+" Up Regulated Total "+str(sum(counts_up2 )))
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        save_title_pi = "plots/pie/"+pi_title +str("_UP_Regulated_PIE_.png")
        # save_title_pi = save_title_pi.replace("_"," ")

        # show plot
        plt.show()
        plt.tight_layout()

        plt.savefig(save_title_pi, dpi=600)

        plt.close()
        plt.cla()
        plt.clf()



        colors_2_use3 = []

        print("done")

        # lipid_classes.sort()

        for i in range(len(down_regulated2  )):
            for j in range(len(colors)):
                if down_regulated2  [i] == colors[j][0]:
                    colors_2_use3.append(colors[j][1])
                    print(colors[j][0],"  ",down_regulated2[i],"    ",colors[j][1])
        print(lipid_classes)
        print(counts)     
        print(ij)
        # exit()
        print(lipid_classes)
        print(colors_2_use)
        fig = plt.figure(figsize =(10, 7))
        plt.pie(counts_down2   , labels = down_regulated2 ,autopct='%.1f%%',colors=colors_2_use3 )
        plt.title(pi_title.replace("_"," ")+" Down Regulated Total "+str(sum(counts_down2 )))
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        save_title_pi = "plots/pie/"+pi_title +str("_down_Regulated_PIE_.png")
        # save_title_pi = save_title_pi.replace("_"," ")

        # show plot
        plt.show()
        plt.tight_layout()

        plt.savefig(save_title_pi, dpi=600)

        plt.close()
        plt.cla()
        plt.clf()



        ###LOG_FC_Sums

#         fig = plt.figure(figsize =(10, 7))
#         plt.pie(sums_of_log2s, labels = lipid_classes,autopct='%.1f%%',colors=colors_2_use )
#         plt.title(pi_Log2_title.replace("_"," "))
#         # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

#         save_title_pi_log2 = "plots/pie/"+pi_Log2_title +str("_PIE.png")
#         # save_title_pi_log2 = save_title_pi_log2.replace("_"," ")

#         # show plot
#         plt.show()
#         plt.savefig(save_title_pi_log2 ,dpi=700)
#         plt.tight_layout()

#         plt.close()
#         plt.cla()
#         plt.clf()


#         fig = plt.figure(figsize =(10, 7))
#         plt.pie(normalized_all_values, labels = lipid_classes,autopct='%.1f%%',colors=colors_2_use )
#         plt.title(pi_Log2_title.replace("_"," ")+"Normalized")
#         # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
# # .replace("_"," ")
#         save_title_pi_log2 = "plots/pie/"+pi_Log2_title +str("Normalized_PIE.png")
#         # save_title_pi_log2 = save_title_pi_log2
#         # show plot
#         plt.show()
#         plt.savefig(save_title_pi_log2 ,dpi=700)
#         plt.tight_layout()

#         plt.close()
#         plt.cla()
#         plt.clf()

#         print(lipid_classes)