
'''
Katrina Kalantar
December 10, 2017

Extracting all the pathogen-associated functions required to run pathogen analysis for the mBAL Study.
These were previously implemented in the .ipynb file, but due to space and visual constraints,
I am re-locating them to a separate file.
'''


import pandas as pd
import sys
import matplotlib.pyplot as plt
import sklearn
from sklearn import metrics
import numpy as np

'''
Fn: calculate_scores(df_rna, df_dna, rna, dna, output_file)
Function to calculate the scores given a matrix from rna and dna and the filenames associated with those dataframes.

input:
df_rna - dataframe containing DeRisi lab pipeline output for RNA file;
df_dna - dataframe containing DeRisi lab pipeline output for associated DNA file;
rna - string filename for rna file (for debugging purposes);
dna - string filename for dna file (for debugging purposes);
output_file - filename for the per-patient top 10 microbes to be written to

output: returns three dataframes, one with just bacterial data / scores, one with just viral, and one with the combined (exception for RNA viruses)
[bacterial_df, viral_df, combined]
'''

def calculate_scores(df_rna, df_dna, rna, dna, output_file, output_directory):

    keep_kingdoms = ['Bacteria','Fungi','Viruses','Escherichia coli']

    ### prepare the RNA dataframe
    df_rna = df_rna[df_rna['Category'].isin(keep_kingdoms)]  #remove kingdoms not in the list of kingdoms

    # RNA must have non-zero expression for all microbes, including DNA viruses;
    # require concordance on NR and NT at genus level
    df_rna = df_rna[df_rna['NR Genus rM'] > 0]
    df_rna = df_rna[df_rna['NT Genus rM'] > 0]

    df_rna = df_rna[['Category','Genus','NT Genus rM','NT Genus Z','Species', 'NT Species rM','NR Genus rM']]   # subset the original dataframe to contain only the relevant columns
    idx = df_rna.groupby(['Genus'])['NT Species rM'].transform(max) == df_rna['NT Species rM']                  # group by Genus, collapsing to take the species with the greatest species-level rM
    df_rna = df_rna[['Category','Genus','NT Genus rM','NT Genus Z','Species','NR Genus rM']][idx]               # keep the species with max rM for each genus
    df_rna.columns = ['Category','Genus','NT Genus rM RNA','NT Genus Z RNA','Species','NR Genus rM RNA']        # re-name the columns

    # write to file for manual verification
    df_rna.sort_values(by='NT Genus rM RNA', inplace=True, ascending=False)   # sort by Genus rM

    #f = open('./output/' + output_directory + '/TableS4.csv','a')
    #f.write("\nRNA\n")
    #f.close()
    #df_rna[['Category','Genus','Species','NT Genus rM RNA','NR Genus rM RNA']].head(n=20).to_csv(output_file, sep='\t', mode='a',index=False)   # write top 20 microbes to the output file

    ### prepare the DNA dataframe
    df_dna = df_dna[df_dna['Category'].isin(keep_kingdoms)]

    df_dna = df_dna[['Category','Genus','NT Genus rM','NT Genus Z','Species', 'NT Species rM','NR Genus rM']]   # subset the original dataframe to contain only the relevant columns
    idx = df_dna.groupby(['Genus'])['NT Species rM'].transform(max) == df_dna['NT Species rM']                  # group by Genus, collapsing to take the species with the greatest species-level rM
    df_dna = df_dna[['Category','Genus','NT Genus rM','NT Genus Z','Species','NR Genus rM']][idx]               # keep the species with max rM for each genus
    df_dna.columns = ['Category','Genus','NT Genus rM DNA','NT Genus Z DNA','Species','NR Genus rM DNA']        # re-name the columns

    # write to file for manual verification
    df_dna.sort_values(by='NT Genus rM DNA', inplace=True, ascending=False)   # sort by Genus rM

    #f = open('./output/' + output_directory + '/TableS4.csv','a')
    #f.write("\nDNA\n")
    #f.close()
    #df_dna[['Category','Genus','Species','NT Genus rM DNA','NR Genus rM DNA']].head(n=20).to_csv(output_file, sep='\t', mode='a',index=False)   # write top 20 microbes to the output file

    # merge DNA and RNA into one matrix with DNA rM and RNA rM values
    master_df = pd.merge(df_dna, df_rna, how='outer', on=['Category','Genus'])   # merge on Genus ID
    master_df.fillna(value=0, inplace=True)    # convert NA values -> 0

    master_df.drop_duplicates(subset=['Genus'],inplace=True)

    # Filter for presence in both RNA and DNA-seq
    no_vir = master_df[master_df['Category']!='Viruses']           # get the subset of microbes that are not viruses (bacterial/fungal rows)
    no_vir_rna = no_vir[no_vir['NT Genus rM RNA'] > 0]             # for bacterial/fungal rows, remove microbes with NT Genus rM RNA = 0
    no_vir_dna = no_vir_rna[no_vir_rna['NT Genus rM DNA'] > 0]     # for bacterial/fungal rows, remove microbes with NT Genus rM DNA = 0
    no_vir_dna = no_vir_dna[no_vir_dna['NR Genus rM RNA'] > 0]             # added 2/13 - for bacterial/fungal rows, remove microbes with NT Genus rM RNA = 0
    no_vir_dna = no_vir_dna[no_vir_dna['NR Genus rM DNA'] > 0]     # added 2/13 - for bacterial/fungal rows, remove microbes with NT Genus rM DNA = 0
    v = master_df[master_df['Category'] == 'Viruses']              # get the subset of microbes that are viruses (no restrictions on RNA and DNA presence)

    dna_viruses_list = ['Mastadenovirus ( 10509 )','Iridovirus ( 10487 )','Roseolovirus ( 40272 )','Betapapillomavirus ( 333922 )','Molluscipoxvirus ( 10278 )','Gammapapillomavirus ( 325455 )',
                        'Alphatorquevirus ( 687331 )','Cytomegalovirus ( 10358 )','Simplexvirus ( 10294 )','Lymphocryptovirus ( 10375 )','Polyomavirus ( 10624 )','Varicellovirus ( 10319 )']
    dna_virs = [i for i in range(len(v.index)) if v.iloc[i]['Genus'] in dna_viruses_list ]
    not_dna_virs = [i for i in range(len(v.index)) if v.iloc[i]['Genus'] not in dna_viruses_list ]
    v_dna_vir = v.iloc[dna_virs]
    v_rna_vir = v.iloc[not_dna_virs]
    #print(v_dna_vir)
    #print(v_rna_vir)
    v_dna_vir = v_dna_vir[v_dna_vir['NT Genus rM DNA'] > 0]
    v_dna_vir = v_dna_vir[v_dna_vir['NT Genus rM RNA'] > 0]
    v_dna_vir = v_dna_vir[v_dna_vir['NR Genus rM DNA'] > 0]   # added 2/13
    v_dna_vir = v_dna_vir[v_dna_vir['NR Genus rM RNA'] > 0]   # added 2/13
    #print(v_dna_vir)

    #temp = pd.concat([no_vir_dna,v], axis=0)                       # concatonate filtered bacterial/fungal with viral reads
    temp = pd.concat([no_vir_dna,v_dna_vir, v_rna_vir])
    master_df = temp

    # create separate output DataFrames for each bacterial, viral, and combined
    bacterial_df = master_df[master_df['Category'] != 'Viruses']
    viral_df = master_df[master_df['Category'] == 'Viruses']
    combined = master_df

    combined.sort_values(by='NT Genus rM RNA', inplace=True, ascending=False)  # sort combined data by NT rM RNA value
    f = open('./output/' + output_directory + '/TableS4.csv','a')
    f.write("\nCombined\n")
    f.close()


    assigned_species_names = []
    for i in combined.index:
        curr = combined.loc[i]['Species_x'] # default is DNA
        if(combined.loc[i]['NT Genus rM DNA'] < combined.loc[i]['NT Genus rM RNA']): #if DNA < RNA
            curr = combined.loc[i]['Species_y']
        assigned_species_names.append(curr)
    combined['Species_Assignment'] = assigned_species_names
    combined[['Category','Genus','Species_Assignment','NT Genus rM DNA','NR Genus rM DNA','NT Genus Z DNA','NT Genus rM RNA','NR Genus rM RNA','NT Genus Z RNA']].head(n=20).to_csv(output_file, sep='\t', mode='a',index=False)
    #combined.head(n=10).to_csv(output_file, sep='\t', mode='a')                # write top 20 microbes to the output file


    #####
    #f = open('./output/' + output_directory + '/TableS4_lowcount.csv','a')
    #f.write("\nRNA\n")
    #f.close()
    #t1 = df_rna[df_rna['NT Genus rM RNA'] < 1]
    #t2 = t1[t1['Species'].isin(list(combined.head(n=20)['Species_Assignment']))]
    #t2[['Category','Genus','Species','NT Genus rM RNA','NR Genus rM RNA']].head(n=20).to_csv("./output/" + output_directory + "/TableS4_lowcount.csv", sep='\t', mode='a',index=False)
    #f = open('./output/' + output_directory + '/TableS4_lowcount.csv','a')
    #f.write("\nDNA\n")
    #f.close()
    #t1 = df_dna[df_dna['NT Genus rM DNA'] < 1]
    #t2 = t1[t1['Species'].isin(list(combined.head(n=20)['Species_Assignment']))]
    #t2[['Category','Genus','Species','NT Genus rM DNA','NR Genus rM DNA']].head(n=20).to_csv("./output/" + output_directory + "/TableS4_lowcount.csv", sep='\t', mode='a',index=False)
    #####


    return [bacterial_df, viral_df, combined]


'''

Fn: main(BMs, investigative = False)
Function to iterate over all RNA/DNA file pairs, for each possible background models, and calculate the scores,
input:
BMs - comma-separated list of background models (ie BM_4,BM_16)
output: returns a dataframe of background models with the corresponding scores for each file

'''

def main(BMs, file_pairs_dict, output_directory, metadata, investigative=False):

    #print("UPDATED")
    #print(file_pairs_dict)
    list_of_BMs = BMs.split(',')
    #print(list_of_BMs)

    full_bacterial_scores = {}
    full_viral_scores = {}
    full_scores = {}

    #for each background model:
    for bm in list_of_BMs:

        #create a list of RNA- and DNA-specific files from the pairs dictionary
        dna_csv_files = ['./' + bm + '/DNA/' + list(file_pairs_dict.keys())[i] +'.report.csv' for i in range(len(file_pairs_dict))]
        rna_csv_files = ['./' + bm + '/RNA/' + list(file_pairs_dict.values())[i] +'.report.csv' for i in range(len(file_pairs_dict))]

        #for each file:
        for c in range(len(rna_csv_files)):  #iterate through all csv files in the directory
            rna = rna_csv_files[c]
            dna = dna_csv_files[c]

            #print('\n\n')
            #print(rna)
            #print(dna)

            #reset these dataframe values to avoid accidentally comparing two different residual frames
            df_rna = pd.DataFrame()
            df_dna = pd.DataFrame()

            #Try to read in the files, if they are not present then write the error and move on.
            try:
                df_rna = pd.read_csv(rna, error_bad_lines=False)
            except:
                print('failed to read RNA file: ' + rna)
                print(sys.exc_info()[0])
                continue
            try:
                df_dna = pd.read_csv(dna, error_bad_lines=False)
            except:
                print('failed to read DNA file: ' + dna)
                print(sys.exc_info()[0])
                continue

            #print(df_dna.head())
            try:
                df_dna = df_dna[df_dna.Genus != 'No genus ( -1 )'] #remove this genus!
                df_rna = df_rna[df_rna.Genus != 'No genus ( -1 )'] #remove this genus!
            except:
                print("DNA or RNA file was malformed - check .csv report file")
                print(sys.exc_info()[0])
                continue

            #print(df_dna.sum(axis=0)['NT Genus rM'])
            #print(df_rna.sum(axis=0)['NT Genus rM'])

            f = open('./output/' + output_directory + '/TableS4.csv','a')   # write the RNA filename to the output file
            f.write("\n\n")
            f.write('TA-'+str(metadata.loc[metadata['RNAfilename'] == rna.split('/')[-1].split('.')[0]].index[0]) + "\n")
            f.close()

            ##### 3/4
            #f = open('./output/' + output_directory + '/TableS4_lowcount.csv','a')   # write the RNA filename to the output file
            #f.write("\n\n")
            #f.write('TA-'+str(metadata.loc[metadata['RNAfilename'] == rna.split('/')[-1].split('.')[0]].index[0]) + "\n")
            #f.close()
            #####

            x = calculate_scores(df_rna, df_dna, rna, dna, './output/' + output_directory + '/TableS4.csv', output_directory)  #calculate scores for the rna and dna files

            f = open('./output/' + output_directory + '/TableS4.csv','a')   # create a break in the output file
            f.write("\n\n")
            f.close()

            # get the scores from the list output by calculate_scores() function
            bacterial_scores = x[0]
            viral_scores = x[1]
            combined_scores = x[2]

            # append the results for each file to the output
            full_bacterial_scores[rna] = bacterial_scores
            full_viral_scores[rna] = viral_scores
            full_scores[rna] = combined_scores

    return [full_bacterial_scores, full_viral_scores, full_scores]






'''
Fn: apply_rule_based_method_cluster(input_data, input_res, groupID, output_filename, title=None)
Function to run the Rules-Based Method for pathogen identification
input:
input_rbm - a dataframe containing the reference information for each microbe/patient (this was the training data for the model(s))
input_res - these are the true classifications corresponding to the training data
groupID - group that is being trained / tested on
output_filename - the desired output filename (for the .pdf)
title - title for the output plot of predictions
output: returns the predicted microbes per patient (input to heatmap functions)
'''


#THIS IS THE ONE TO USE!!!!

def get_split(list_of_values):  #returns the index
    splits = {}
    max_split = 0
    max_split_id = 0
    split_list = []
    for i in range(len(list_of_values)):
        if i == 0:
            splits[i] = 0
            split_list.append(0)
        else:
            splits[i] = -1 * (list_of_values[i] - list_of_values[i-1] )
            split_list.append(splits[i])

        if splits[i] > max_split:
            max_split = splits[i]
            max_split_id = i

    return max_split_id


def apply_rule_based_method_cluster(input_data, input_res, groupID, output_filename, title=None, annotate_plots=False,rbm_method_plot=False):

    rbm1_patient = []
    rbm1_genus = []
    rbm1_microbe = []

    s = int(np.round(len(set(input_data['patient']))/5)+1)

    if rbm_method_plot:
        fig, axarr = plt.subplots(s,5, figsize=(20, s*4), sharex=True, sharey=True)
        row = 0
        col = 0

    for i in list(set(input_data['patient'])):#[0:3]:
        new_adds = 0
        #print(i)

        sub = input_data[input_data['patient']==i]
        sub['rna+dna'] = sub['RNAvalue'] + sub['DNAvalue']
        sub.sort_values(by='rna+dna', ascending=False, inplace=True)

        #ADD POTENTIAL BACTERIAL MICROBES
        subB = sub[sub['pathogenic_green'] == False] #cannot be a virus

        size_of_top_clust = get_split(list(subB['rna+dna']))
        # this was added on 12/28 to deal with cases in which only 1 bacterial microbe is present
        if(size_of_top_clust < 1):
            size_of_top_clust = 1 # there should always be at least 1 microbe in top cluster (if only 1 microbe present)


        ## NEW PLOTTING FUNCTION 2/1
        #plt.figure(figsize=[3,3])
        if rbm_method_plot:
            axarr[row,col].scatter([i for i in range(len(subB['rna+dna']))],subB['rna+dna'],color=[['darkblue','red'][int(i)] for i in subB['pathogenic_red']])
            axarr[row,col].axvline(x=size_of_top_clust - .5, color='r', linestyle='-',lw=.4)
            axarr[row,col].set_title(i)
            axarr[row,col].set_xlim(-1,15)
            axarr[row,col].set_ylim(-1,12)
            axarr[row,col].set_xticks(np.arange(0,15,1)) #MODS 2/1
            axarr[row,col].set_yticks(np.arange(0,12,1)) #MODS 2/1
            for l in range(len(subB.index)):
                axarr[row,col].annotate(subB.iloc[l]['microbe'],(l + .1,subB.iloc[l]['rna+dna']-.1),
                                 fontsize=8, color='grey')
            if col == 4:
                row += 1
                col = 0
            else:
                col += 1


        #print(size_of_top_clust)
        top = subB.head(n=size_of_top_clust)
        top = top[top['pathogenic_red'] == True]
        #print(top)
        for j in range(len(top.index)):
            rbm1_genus.append(str(top['microbe_genus'].iloc[j]))
            rbm1_microbe.append(str(top['microbe'].iloc[j]))
            new_adds += 1

        possible_viruses = sub[[(sub.loc[j]['pathogenic_green'] == True and sub.loc[j]['pathogenic_red'] == True) for j in sub.index]]
        possible_viruses = possible_viruses[possible_viruses['RNAvalue'] > np.log10(.1 + 1)] # RNAvalue is log10(RNA rpm + 1), and we want to keep viruses with rpm > .1
        if(len(possible_viruses.index) > 0):
            res = possible_viruses.head(n=1)
            for j in range(len(res.index)):
                rbm1_genus.append(str(res['microbe_genus'].iloc[j]))
                rbm1_microbe.append(str(res['microbe'].iloc[j]))
                new_adds += 1

        rbm1_patient = rbm1_patient + [i for n in range(new_adds)]

    rbm_g1 = pd.DataFrame.from_dict({'microbe': rbm1_microbe,'patient':rbm1_patient,'microbe_genus':rbm1_genus})
    rbm_g1.drop_duplicates(inplace=True)

    if groupID == 2:
        c = ['blue','blue']
    else:
        c = ['blue','red']

    if rbm_method_plot:
        plt.savefig(output_filename+"_methodplot.pdf")
        plt.show()

    #print([c[int(i)] for i in input_res])

    plt.figure(figsize=[8,8])
    a = [rbm_g1.iloc[i]['microbe'] + "-" + str(rbm_g1.iloc[i]['patient']) for i in range(len(rbm_g1))]
    match = [int(input_data.iloc[i]['microbe'] + "-" + input_data.iloc[i]['patient'] in a) for i in range(len(input_data.index))]
    plt.scatter(input_data['RNAvalue'],input_data['DNAvalue'],
                edgecolor=[c[int(input_data.loc[i]['positive'])]
                       if input_data[input_data['patient'] == input_data.loc[i]['patient']]['positive'].sum() < 2
                       or c[int(input_data.loc[i]['positive'])] == 'blue'
                       or input_data.loc[i]['RNAvalue'] == max(input_data[(input_data.patient == input_data.loc[i]['patient']) & (input_data.positive == True)]['RNAvalue'])
                       else 'orange' for i in input_data.index],
                s=[85 if match[i] == 1 else 50 for i in range(len(match)) ],linewidth=1.2,
                facecolor = [c[input_res[i]] if match[i] == 1 else 'white' for i in range(len(match)) ], alpha=.7)
                #facecolor=[c[int(Y_g2[i])] if predicted[:,1][i] > probability_threshold else 'white' for i in range(len(Y_g2)) ],alpha=.6) #s=predicted[:,1]*200
    if(annotate_plots):
        for i in range(len(input_res)):
            #if(predicted[:,1][i] > probability_threshold):
            if input_res[i] or match[i] == 1: #if the microbe was True in input_res
                try:
                    sp = input_data.iloc[i]['microbe'].split(' ')
                except:
                    sp = 'NA'
                pa = input_data.iloc[i]['patient'].split('-')[1]
                plt.annotate(sp[0][0] + '. ' + sp[1] + ' - ' + pa,(input_data.iloc[i]['RNAvalue'] + .1,input_data.iloc[i]['DNAvalue']-.1),
                             fontsize=5, color='grey')

    plt.xlim((-.2, 4.5))
    plt.ylim((-.2, 4.5))
    if not title == None:
        plt.title(title)
    else:
        plt.title("Group " + str(groupID) + " - Rule Based Method")
    #plt.axis((-1,12,-1,12))
    plt.xlabel('log( RNA rpM )' )
    plt.ylabel('log( DNA rpM )' )
    plt.savefig(output_filename + ".pdf")

    return rbm_g1


#rbm_g1 = apply_rule_based_method_cluster(X_g1, Y_g1, 1, "rbm_1")
#rbm_g1g4 = apply_rule_based_method_cluster(pd.concat([X_g1,X_g4]), Y_g1+Y_g4, 1, './output/' + output_directory + '/Figure2B', title="Rules-Based Model Predictions\n for LRTI+C/M and LRTI-NEG microbes")


# Function: apply_euclidean_RBM
#
# This is a separate parallel implementation of the RBM with euclidean distance in RNA x DNA space (as opposed to the change in RNA+DNA rpm)
#

from scipy.spatial import distance

def get_euclidean_split(df):  #returns the index
    dmat = scipy.spatial.distance.cdist(10**df[['RNAvalue','DNAvalue']],10**df[['RNAvalue','DNAvalue']],metric='euclidean')
    newDF = pd.DataFrame(dmat,columns=[i for i in range(len(df.index))])
    #newDF.reindex(columns=[i for i in range(len(newDF.index))])
    #print(newDF)
    a = [newDF.iloc[i][i+1] for i in range(len(newDF.index)-1)]
    #print(a)
    keep_this_many = 0 # default
    if(len(a) > 0):
        keep_this_many = max( (v, i) for i, v in enumerate(a) )[1]+1  #adjust default if there is more than 1 value in the list

    return keep_this_many

def apply_euclidean_RBM(input_data, input_res, groupID, output_filename, title=None, annotate_plots=False,rbm_method_plot=False):

    rbm1_patient = []
    rbm1_genus = []
    rbm1_microbe = []

    s = int(np.round(len(set(input_data['patient']))/5)+1)

    if rbm_method_plot:
        fig, axarr = plt.subplots(s,5, figsize=(20, s*4), sharex=True, sharey=True)
        row = 0
        col = 0

    for i in list(set(input_data['patient'])):#[0:3]:
        new_adds = 0
        print(i)

        sub = input_data[input_data['patient']==i]
        sub['rna+dna'] = sub['RNAvalue'] + sub['DNAvalue']
        sub.sort_values(by='rna+dna', ascending=False, inplace=True)

        #ADD POTENTIAL BACTERIAL MICROBES
        subB = sub[sub['pathogenic_green'] == False] #cannot be a virus

        size_of_top_clust = get_euclidean_split(subB)
        # this was added on 12/28 to deal with cases in which only 1 bacterial microbe is present
        if(size_of_top_clust < 1):
            size_of_top_clust = 1 # there should always be at least 1 microbe in top cluster (if only 1 microbe present)


        ## NEW PLOTTING FUNCTION 2/1
        #plt.figure(figsize=[3,3])
        if rbm_method_plot:
            axarr[row,col].scatter([i for i in range(len(subB['rna+dna']))],subB['rna+dna'],color=[['darkblue','red'][int(i)] for i in subB['pathogenic_red']])
            axarr[row,col].axvline(x=size_of_top_clust - .5, color='r', linestyle='-',lw=.4)
            axarr[row,col].set_title(i)
            axarr[row,col].set_xlim(-1,15)
            axarr[row,col].set_ylim(-1,12)
            axarr[row,col].set_xticks(np.arange(0,15,1)) #MODS 2/1
            axarr[row,col].set_yticks(np.arange(0,12,1)) #MODS 2/1
            for l in range(len(subB.index)):
                axarr[row,col].annotate(subB.iloc[l]['microbe'],(l + .1,subB.iloc[l]['rna+dna']-.1),
                                 fontsize=8, color='grey')
            if col == 4:
                row += 1
                col = 0
            else:
                col += 1


        #print(size_of_top_clust)
        top = subB.head(n=size_of_top_clust)
        top = top[top['pathogenic_red'] == True]
        #print(top)
        for j in range(len(top.index)):
            rbm1_genus.append(str(top['microbe_genus'].iloc[j]))
            rbm1_microbe.append(str(top['microbe'].iloc[j]))
            new_adds += 1

        possible_viruses = sub[[(sub.loc[j]['pathogenic_green'] == True and sub.loc[j]['pathogenic_red'] == True) for j in sub.index]]
        possible_viruses = possible_viruses[possible_viruses['RNAvalue'] > np.log10(.1 + 1)] # RNAvalue is log10(RNA rpm + 1), and we want to keep viruses with rpm > .1
        if(len(possible_viruses.index) > 0):
            res = possible_viruses.head(n=1)
            for j in range(len(res.index)):
                rbm1_genus.append(str(res['microbe_genus'].iloc[j]))
                rbm1_microbe.append(str(res['microbe'].iloc[j]))
                new_adds += 1

        rbm1_patient = rbm1_patient + [i for n in range(new_adds)]

    rbm_g1 = pd.DataFrame.from_dict({'microbe': rbm1_microbe,'patient':rbm1_patient,'microbe_genus':rbm1_genus})
    rbm_g1.drop_duplicates(inplace=True)

    if groupID == 2:
        c = ['blue','blue']
    else:
        c = ['blue','red']

    if rbm_method_plot:
        plt.show()

    #print([c[int(i)] for i in input_res])

    plt.figure(figsize=[8,8])
    a = [rbm_g1.iloc[i]['microbe'] + "-" + str(rbm_g1.iloc[i]['patient']) for i in range(len(rbm_g1))]
    match = [int(input_data.iloc[i]['microbe'] + "-" + input_data.iloc[i]['patient'] in a) for i in range(len(input_data.index))]
    plt.scatter(input_data['RNAvalue'],input_data['DNAvalue'],
                edgecolor=[c[int(input_data.loc[i]['positive'])]
                       if input_data[input_data['patient'] == input_data.loc[i]['patient']]['positive'].sum() < 2
                       or c[int(input_data.loc[i]['positive'])] == 'blue'
                       or input_data.loc[i]['RNAvalue'] == max(input_data[(input_data.patient == input_data.loc[i]['patient']) & (input_data.positive == True)]['RNAvalue'])
                       else 'orange' for i in input_data.index],
                s=[85 if match[i] == 1 else 50 for i in range(len(match)) ],linewidth=1.2,
                facecolor = [c[input_res[i]] if match[i] == 1 else 'white' for i in range(len(match)) ], alpha=.7)
                #facecolor=[c[int(Y_g2[i])] if predicted[:,1][i] > probability_threshold else 'white' for i in range(len(Y_g2)) ],alpha=.6) #s=predicted[:,1]*200
    if(annotate_plots):
        for i in range(len(input_res)):
            #if(predicted[:,1][i] > probability_threshold):
            if input_res[i] or match[i] == 1: #if the microbe was True in input_res
                try:
                    sp = input_data.iloc[i]['microbe'].split(' ')
                except:
                    sp = 'NA'
                pa = input_data.iloc[i]['patient'].split('-')[1]
                plt.annotate(sp[0][0] + '. ' + sp[1] + ' - ' + pa,(input_data.iloc[i]['RNAvalue'] + .1,input_data.iloc[i]['DNAvalue']-.1),
                             fontsize=5, color='grey')

    plt.xlim((-.2, 4.5))
    plt.ylim((-.2, 4.5))
    if not title == None:
        plt.title(title)
    else:
        plt.title("Group " + str(groupID) + " - Rule Based Method")
    #plt.axis((-1,12,-1,12))
    plt.xlabel('log( RNA rpM )' )
    plt.ylabel('log( DNA rpM )' )
    plt.savefig(output_filename + ".pdf")

    return rbm_g1







'''
Fn: perf_measure( truth_data, predicted_data )
Function to count true positives, false positives, true negatives, and false negatives for a given prediction
output: returns an array of [TP, FP, TN, FN]
'''

# basic fn to get the components of sensitivity and specificity
def perf_measure(y_actual, y_hat):
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    for i in range(len(y_hat)):
        if y_actual[i]==y_hat[i]==1:
           TP += 1
    for i in range(len(y_hat)):
        if y_hat[i]==1 and y_actual[i]!=y_hat[i]:
           FN += 1
    for i in range(len(y_hat)):
        if y_actual[i]==y_hat[i]==0:
           TN += 1
    for i in range(len(y_hat)):
        if y_hat[i]==0 and y_actual[i]!=y_hat[i]:
           FP += 1

    return(TP, FP, TN, FN)



def Find_Optimal_Cutoff(target, predicted):
    """ Find the optimal probability cutoff point for a classification model related to event rate
    Parameters
    ----------
    target : Matrix with dependent or target data, where rows are observations

    predicted : Matrix with predicted data, where rows are observations

    Returns
    -------
    list type, with optimal cutoff value

    """
    fpr, tpr, threshold = metrics.roc_curve(target, predicted,pos_label=1)
    i = np.arange(len(tpr))
    roc = pd.DataFrame({'tf' : pd.Series(tpr-(1-fpr), index=i), 'threshold' : pd.Series(threshold, index=i), 'tpr':pd.Series(tpr), '1-fpr':pd.Series(1-fpr)})
    print(roc)
    roc_t = roc.ix[(roc.tf-0).abs().argsort()[:1]]
    print(roc_t)
    return list(roc_t['threshold'])




# METHOD To get 95% CI - adapted from stack exchange: https://stackoverflow.com/questions/19124239/scikit-learn-roc-curve-with-confidence-intervals

import numpy as np
from scipy.stats import sem
from sklearn.metrics import roc_auc_score

def get_CI(y_true, y_pred, n_bootstraps):

    y_pred = np.array(y_pred)
    y_true = np.array(y_true)

    #print("Original ROC area: {:0.3f}".format(roc_auc_score(y_true, y_pred)))

    rng_seed = 42  # control reproducibility
    bootstrapped_scores = []

    rng = np.random.RandomState(rng_seed)
    for i in range(n_bootstraps):
        # bootstrap by sampling with replacement on the prediction indices
        indices = rng.random_integers(0, len(y_pred) - 1, len(y_pred))
        if len(np.unique(y_true[indices])) < 2:
            # We need at least one positive and one negative sample for ROC AUC
            # to be defined: reject the sample
            continue

        score = roc_auc_score(y_true[indices], y_pred[indices])
        bootstrapped_scores.append(score)
        #print("Bootstrap #{} ROC area: {:0.3f}".format(i + 1, score))

    #calculate CI based on the bootstrapped_scores
    sorted_scores = np.array(bootstrapped_scores)
    sorted_scores.sort()

    # Computing the lower and upper bound of the 90% confidence interval
    # You can change the bounds percentiles to 0.025 and 0.975 to get
    # a 95% confidence interval instead.
    confidence_lower = sorted_scores[int(0.05 * len(sorted_scores))]
    confidence_upper = sorted_scores[int(0.95 * len(sorted_scores))]
    print("Original ROC area: {:0.3f}".format(roc_auc_score(y_true, y_pred)) + ", [{:0.3f} - {:0.3}]".format(
        confidence_lower, confidence_upper))
    return([confidence_lower, confidence_upper])
