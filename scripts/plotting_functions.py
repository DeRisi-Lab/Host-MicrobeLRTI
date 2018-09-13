'''
Katrina Kalantar
December 10, 2017

Extracting all the plotting functions required to plot results from the pathogen analysis in the mBAL Study.
These were previously implemented in the .ipynb file, but due to space and visual constraints,
I am re-locating them to a separate file.
'''

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns

'''
Fn: evaluate_sensitivity(input_df)
Function to determine the sensitivity and specificity of a classification method
input:
input_df - dataframe containing predicted values for each patient
output: returns nothing, but runs a series of print statements to give information on the sensitivity / specificity
'''

### Get the sensitivity/specificity for the microbes identified by the model [THIS IS ACTUALLY VERY SIMILAR TO THE CODE ABOVE...STILL IN PROGRESS]

def evaluate_sensitivity(input_df, list_of_files, metadata, true_classification, full, output_directory):

    #print(input_df.head())

    sensitivity_specificity_df = {}
    #sensitivity_specificity_df[microbe = {'ID in culture'}]

    had_cultured_organism = 0
    at_least_one_organism_identified_by_ngs = 0
    all_organisms_identified_by_ngs = 0
    patients_where_all_NGS_microbes_were_positive_by_model = 0
    #print(input_df.head())

    for i in list_of_files:
        sample_id = i.split('/')[-1].split('.')[-3]

        #if '209' in sample_id:
        #    continue
        metadata_id = metadata.loc[metadata['RNAfilename'] == sample_id].index[0]
        tc = true_classification.loc[[sample_id],:]['effective_group'][0]

        bac = full[i]

        try:
            organism_list = metadata.loc[metadata_id]['organism'].split(',') #will fail if nan, so since it doesn't fail there was a cultured organism
            ID_by_MODEL = list(input_df[input_df['patient']==i.split('/')[-1]]['microbe_genus'])
            ID_by_NGS = [bac.loc[i]['Genus'] for i in bac.index if bac.loc[i]['Genus'] in organism_list]

            for m in [set(organism_list + ID_by_MODEL)]:
                for n in m:
                    sensitivity_specificity_df[sample_id + " / " + n] = {'ID in culture': (n in organism_list), 'ID in NGS':(n in ID_by_NGS or n in ID_by_MODEL), 'ID in Model': (n in ID_by_MODEL)}


        except:
            if tc == 1:
                print(metadata.loc[metadata_id]['organism'])  #if tc is 1 and no organism is listed, we want to know...otherwise skip

    #print("# of patients that had a cultured organism: " + str(had_cultured_organism))
    #print('# of patients that had at least one cultured organism identified in NGS: ' + str(at_least_one_organism_identified_by_ngs))
    #print('# of patients that had all cultured organisms identified in NGS: ' + str(all_organisms_identified_by_ngs))

    #print(pd.DataFrame(sensitivity_specificity_df).transpose())#.head())

    SESP_DF = pd.DataFrame(sensitivity_specificity_df).transpose()
    SESP_DF.to_csv('./output/' + output_directory + '/sensitivity_specificity_df.csv')

    print("Total # of microbes identified in culture")
    print(SESP_DF[SESP_DF['ID in culture']].shape[0])
    print("Number of microbes identified in culture that were also identified by RBM:")
    print(sum(SESP_DF[SESP_DF['ID in culture'] == True]['ID in Model']))

    print("Number of patients with at least one microbe identified by RBM:")
    print(len(set([i.split('/')[0] for i in SESP_DF[SESP_DF['ID in culture']].index])))

    culture_pos_microbes = SESP_DF[SESP_DF['ID in culture']==True]
    print(culture_pos_microbes.sum())

    #number of individuals with at least 1 microbe identified in culture and identified by the model...
    print(len(set([i.split('/')[0] for i in culture_pos_microbes[culture_pos_microbes['ID in Model']==True].index])))


    model_pos_microbes = SESP_DF[SESP_DF['ID in Model']==True]
    print(model_pos_microbes[model_pos_microbes['ID in culture'] == False])
    print(model_pos_microbes[model_pos_microbes['ID in culture'] == False].shape)
    print(model_pos_microbes.sum())






'''
Fn: create_prediction_heatmap_combo(input_data, original_matrix, filename, groupID)
Function to create a heamap of predicted microbes by patient ID, where 2 or more groups are included in the patients to be plotted (useful for comparing when training on both group 1 and group 4)
input:
input_data - a dataframe containing the patient, microbe genus, and microbe species for every microbe predicted by one of the two pathogen v. commensal models
original_matrix - a dataframe containing the reference information for each microbe/patient (this was the training data for the model)
filename - the desired output filename (for the .pdf)
groupID - the group (1,2,3,4 indicating LRTI+C/M, LRTI+C, UNK, LRTI-NEG respectively) for which classification a particular patient came from
output: returns nothing, but outputs a .pdf heatmap of microbes by patient ID
'''

def create_prediction_heatmap_combo(input_data, original_matrix, filename, groupID, true_classification):

    CustomCmap = matplotlib.colors.ListedColormap(['white','darkblue','red'])
    yay_v = 0

    d = pd.DataFrame(0, index=np.arange(len(set(input_data['microbe']))), columns=set(original_matrix['patient']))
    d.index = set(input_data['microbe'])

    for i in input_data.index:
        d[input_data.iloc[i]['patient']][input_data.iloc[i]['microbe']] = 1
        orgs = metadata[metadata['RNAfilename'] == input_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]

        if(type(orgs) == type('string')):
            if input_data.iloc[i]['microbe_genus'] in orgs.split(','):
                yay_v += 1
                d[input_data.iloc[i]['patient']][input_data.iloc[i]['microbe']] = 2

    #loop through every patient
    #add all group 1s
    g1s = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 1:
                g1s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]

    #add all group 4s
    g4s = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 4:
                g4s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]

    gcombos = g1s + g4s
    fig_height = d.shape[0]/4.8
    fig_width = d.shape[1]/4.8

    d_1 = d

    d = d[gcombos] # this sorts the columns in order of [group1] [group4]

    #reset names to remove .report.csv notation
    new_col_names = []
    for i in d.columns:
        new_col_names.append(i.split('.')[0])
    d.columns = new_col_names

    #d.to_csv('LRcombo.csv')

    plt.figure(figsize=(fig_width,fig_height))
    #plt.figure(figsize=(5.5,5))
    ax = sns.heatmap(d,cmap=CustomCmap,vmin=0,vmax=2,cbar=False)
    ax.hlines([i for i in range(len(d.index))], *ax.get_xlim(),lw=.4)
    ax.vlines([i for i in range(len(d.columns))], *ax.get_ylim(),lw=.4)
    ax.tick_params(labelsize=10)
    plt.title( filename + ", by patient ID")
    plt.savefig(filename + ".pdf", bbox_inches='tight')





'''

Fn: create_prediction_heatmap_combo_overlap(input_rbm_data, input_lr_data, original_matrix, filename, groupID, print_output = False)
Function to create a basic heamap of predicted microbes by patient ID
input:
input_rbm - a dataframe containing the patient, microbe genus, and microbe species for every microbe predicted by the RBM
input_lr_data - a dataframe containing the predicted microbes output by the LRM
original_matrix - a dataframe containing the reference information for each microbe/patient (this was the training data for the model(s))
filename - the desired output filename (for the .pdf)
groupID - PARAMETER NO LONGER IN USE
print_output - default = False, but set to True if you want to print out the output
output: returns nothing, but outputs a .pdf heatmap of microbes by patient ID, combining the predictions from both RBM and LRM

'''



# THIS WILL BE USED AS HEATMAP FOR FINAL COMBINED MATRIX!!!

# THIS IS THE VERSION OF THE ABOVE CODE THAT I AM MODIFYING ON 12/08!!!!!!!

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

# TRYING TO CREATE COMBO HEATMAP THAT WILL TAKE IN THE INPUTS FROM ALL G1, G4, G2/3 and output the color-coded heatmap
# NOVEMBER 27

def create_prediction_heatmap_combo_overlap(input_rbm_data, input_lr_data, original_matrix, filename, groupID, true_classification, output_directory, metadata, training_sample_names, output_filename, print_output = False, collapse=True):

    input_rbm_data.to_csv('./output/' + output_directory + '/TableS6A.csv')
    input_lr_data.to_csv('./output/' + output_directory + '/TableS6B.csv')

    # my original colors - using color to indicate group, shade to indicate clin +
    #colors = ['white','deepskyblue','lightgrey','blue','yellow','orange','gold','red']
    #colors = ['white','#2FD2FF','lightgrey','#2666C7','#FFFC6F','#FF7E70','#FFD700','#DA2727']

    # new colors - using shade to indicate group, color to indicate clin +
    # WHITE, RB only (light blue), clin only (light grey), RB + Clin (light red), LR only (med blue), LR + RM (dark blue), LR + Clin (med red), LR + RB + Clin (dark red)
    #colors = ['white','#819BFC','lightgrey','#FFABAB','#2752F3','#0026B2','#F30B0B','#AE0000']
    colors = ['white','#AFECFF','lightgrey','#FFA2A2','#0BAEFF','#0026B2','#FF3232','#C00000'] # trying to make the above colors more distinct

    CustomCmap = matplotlib.colors.ListedColormap(colors)
    yay_v = 0

    d = pd.DataFrame(0, index=np.arange(len(set(input_rbm_data['microbe']))), columns=set(original_matrix['patient']))
    d.index = set(input_rbm_data['microbe'])


    #loop through every patient and create a new column for each patient
    #add all group 1s
    g1s_train = []
    g1s_test = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 1:
                if(i+'.report.csv' in training_sample_names):
                    g1s_train.append(i+'.report.csv')
                else:
                    g1s_test.append(i+'.report.csv')
                #g1s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]
    #add all group 4s
    g4s_train = []
    g4s_test = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 4:
                if(i+'.report.csv' in training_sample_names):
                    g4s_train.append(i+'.report.csv')
                else:
                    g4s_test.append(i+'.report.csv')
                #g4s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]
    #add all group 2s
    g2s = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 2:
                g2s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]
    #add all group 3s
    g3s = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 3:
                g3s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]



    # ADD VALUES FOR RB method
    for i in input_rbm_data.index:
        d[input_rbm_data.iloc[i]['patient']][input_rbm_data.iloc[i]['microbe']] = 1
        orgs = metadata[metadata['RNAfilename'] == input_rbm_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]
        #print(orgs)
        if(type(orgs) == type('string')):
            if input_rbm_data.iloc[i]['microbe_genus'] in orgs.split(','):
                yay_v += 1
                d[input_rbm_data.iloc[i]['patient']][input_rbm_data.iloc[i]['microbe']] = d[input_rbm_data.iloc[i]['patient']][input_rbm_data.iloc[i]['microbe']] + 2

    # ADD VALUES FOR LR method
    for i in input_lr_data.index:

        try:
            d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] = d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] + 4   #try to just add to existing values
        except:
            #NEED TO ADD THE ROW
            new_index = [ v for v in d.index] + [input_lr_data.iloc[i]['microbe']]
            d.loc[len(d)]=[0 for i in range(len(d.columns))]
            d.index = new_index
            d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] = 4   # if the value has not existed before, add to the dataframe

        orgs = metadata[metadata['RNAfilename'] == input_lr_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]   # loop through all the organisms
        if(type(orgs) == type('string')):
            if input_lr_data.iloc[i]['microbe_genus'] in orgs.split(','):
                yay_v += 1
                if d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] not in [3,7]:  # if the organism is confirmed, but has not previously been clinically confirmed then...add clinical to it
                    d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] = d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] + 2

    # LOOP THROUGH AND ADD CLINICAL GENERA
    for i in input_lr_data.index:
        orgs = metadata[metadata['RNAfilename'] == input_lr_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]
        if(type(orgs) == type('string')):
            #print(input_lr_data.iloc[i]['patient'])
            #print(orgs)
            for o in orgs.split(','):
                list_of_species = list(set(input_lr_data.loc[input_lr_data['microbe_genus']==o]['microbe']).union(set(input_lr_data.loc[input_lr_data['microbe_genus']==o]['microbe'])))
                a_species_was_found = False
                for s in list_of_species:
                    if not a_species_was_found:
                        try:
                            if d[input_lr_data.iloc[i]['patient']][s] > 0:
                                a_species_was_found = True
                        except:
                            continue

                #print(a_species_was_found)

                if not a_species_was_found:
                    if sum([int(o == v ) for v in  d.index]) > 0:
                        d[input_lr_data.iloc[i]['patient']][o] = 2
                    else:
                        new_index = [ v for v in d.index] + [o]
                        d.loc[len(d)]=[0 for i in range(len(d.columns))]

                        d.index = new_index
                        d[input_lr_data.iloc[i]['patient']][o] = 2

    gcombos = g1s_train + g1s_test + g4s_train +  g4s_test + g2s + g3s

    d = d[gcombos] # this sorts the columns in order of [group1] [group4]

    #reset names to remove .report.csv notation
    new_col_names = []
    for i in d.columns:
        new_col_names.append(i.split('.')[0])
    d.columns = new_col_names




    # Apply some collapsing for the indices
    if(collapse):
        strep_viridans_group = ['Streptococcus pseudopneumoniae ( 257758 )','Streptococcus sp. VT 162 ( 1419814 )',
                                'Streptococcus salivarius ( 1304 )','Streptococcus mitis ( 28037 )',
                                'Streptococcus pseudopneumoniae ( 257758 )','Streptococcus anginosus ( 1328 )']
        strep_preprocessed_names = []
        for i in d.index:
            if(i in strep_viridans_group):
                strep_preprocessed_names.append('Streptococci viridans Group ( XXX )')
            elif('coronavirus' in i):
                strep_preprocessed_names.append('Human coronavirus species ( XXX )')
            else:
                strep_preprocessed_names.append(i)  # strep pneumo and strep genera (clinical) will be merged

        d.index = strep_preprocessed_names
        first_values = [' '.join(i.split(' ')[0:-4]) if i.count(' ') > 3 else ' '.join(i.split(' ')[0:-3]) for i in d.index]
        d.index = first_values


        for j in d.columns:
            for k in set(first_values):
                #print(k)
                l = list(d[j][d.index==k])
                #print(l)
                if 2 in l and 3 in l:
                    # this is bad - it will collapse incorrectly and appear the BM and LR matched (score = 5),
                    # when in reality it was a collapse of RB+Clin (3) and Clin+ (3) again - need to modify the cases post-hoc.
                    #print([3] + [0 for i in range(len(l)-1)])
                    d[j][d.index==k] = [3] + [0 for i in range(len(l)-1)]





        d = d.groupby(d.index).sum()

        post_mod_index = []  # MODIFICATIONS TO UPDATE NAME OF STREP AND STAPH TO BE CLEAR
        for i in d.index:
            if(i == 'Streptococcus'):
                post_mod_index.append('Streptococcus pneumoniae')
            elif(i == 'Staphylococcus'):
                post_mod_index.append('Staphylococcus aureus')
            elif(i == 'Human respiratory syncytial'):
                post_mod_index.append('Human RSV')
            else:
                post_mod_index.append(i)

        d.index = post_mod_index




    # sort the indices by alphabetical order
    new_d = d.sort_index(axis=0)
    d = new_d

    #new_columns = ['-'.join(c.split('-')[0:2]) for c in d.columns]
    #new_cols_temp = ['TA-' + str(c.split('-')[1]) for c in d.columns]  #convert names to TA-XXX notation
    #new_columns = [i[0:6] for i in new_cols_temp]  #get rid of TRNA from the one sample that has weird naming convention
    #print(new_columns)
    #d.columns = new_columns




    fig_height = (d.shape[0]/4.8)*2
    #fig_width = d.shape[1]/4.8
    fig_width = (d.shape[1]/6.4)*2

    #d.to_csv('LRcombo.csv')

    if print_output:
        print(d.head())

    #colors = ['white','#4AA2B6','lightgrey','#04819E','#69DE56','#F93846','#1DD300','#F90012']
    #ordered_colors = ['white','lightgrey','deepskyblue','blue','yellow','gold','orange','red']

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))


    plt.figure(figsize=(fig_width,fig_height))
    #ax = sns.heatmap(d,cmap=CustomCmap,vmin=0,vmax=7,cbar=False)
    ax = sns.heatmap(d,cmap=cmap,vmin=0,vmax=7)#,cbar=False)
    ax.hlines([i for i in range(len(d.index))], *ax.get_xlim(),lw=.4)
    ax.vlines([i for i in range(len(d.columns))], *ax.get_ylim(),lw=.4)
    ax.tick_params(labelsize=10)
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([i for i in range(len(colors))])
    colorbar.set_ticks([0,2,1,3,4,6,5,7])
    colorbar.set_ticklabels(['NA','Clin+','RB+','RB+,Clin+','LR+','LR+,Clin+','LR+,RB+','LR+,RB+,Clin+'])
    #colorbar.set_ticklabels(['NA','RB+','Clin+','RB+,Clin+','LR+','LR+,RB+','LR+,Clin+','LR+,RB+,Clin+'])

    plt.title( filename + ", by patient ID")
    plt.savefig('./output/' + output_directory + '/' + output_filename , bbox_inches='tight')

    return d

#for testing purposes
#d = create_prediction_heatmap_combo_overlap(c,lr_results,X_full, './output/' + output_directory + '/Combined Predictions -\n Leave Out One Patient per Iteration',[1,4,2,3],print_output = True)



# this plots the heatmap group by group, separately showing clinical only v. RBM/LR results - ulimately goes into FigS2

def create_prediction_heatmap_combo_overlap_update2(input_rbm_data, input_lr_data, original_matrix, filename, groupID, true_classification, output_directory, metadata, training_sample_names, output_filename, print_output = False, collapse=True):

    input_rbm_data.sort_values(by='patient',inplace=True)
    input_rbm_data.to_csv('./output/' + output_directory + '/TableS6A.csv')
    input_lr_data.sort_values(by='patient',inplace=True)
    input_lr_data.to_csv('./output/' + output_directory + '/TableS6B.csv')

    # my original colors - using color to indicate group, shade to indicate clin +
    #colors = ['white','deepskyblue','lightgrey','blue','yellow','orange','gold','red']
    #colors = ['white','#2FD2FF','lightgrey','#2666C7','#FFFC6F','#FF7E70','#FFD700','#DA2727']

    # new colors - using shade to indicate group, color to indicate clin +
    # WHITE, RB only (light blue), clin only (light grey), RB + Clin (light red), LR only (med blue), LR + RM (dark blue), LR + Clin (med red), LR + RB + Clin (dark red)
    #THESE are the colors we have used for a long time, prior to 2x plot vvv
    #colors = ['white','#AFECFF','lightgrey','#FFA2A2','#0BAEFF','#0026B2','#FF3232','#C00000'] # trying to make the above colors more distinct
    colors = ['white','#FFA2A2','lightgrey','#FFA2A2','#FF3232','#C00000','#FF3232','#C00000']  # trying to make monochromatic FigB

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    CustomCmap = matplotlib.colors.ListedColormap(colors)
    yay_v = 0

    d = pd.DataFrame(0, index=np.arange(len(set(input_rbm_data['microbe']))), columns=set(original_matrix['patient']))
    d.index = set(input_rbm_data['microbe'])


    #loop through every patient and create a new column for each patient
    #add all group 1s
    g1s_train = []
    g1s_test = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 1:
                if(i+'.report.csv' in training_sample_names):
                    g1s_train.append(i+'.report.csv')
                else:
                    g1s_test.append(i+'.report.csv')
                #g1s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]
    #add all group 4s
    g4s_train = []
    g4s_test = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 4:
                if(i+'.report.csv' in training_sample_names):
                    g4s_train.append(i+'.report.csv')
                else:
                    g4s_test.append(i+'.report.csv')
                #g4s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]
    #add all group 2s
    g2s = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 2:
                g2s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]
    #add all group 3s
    g3s = []
    for i in true_classification.index:
        if 'RNA' in i:
            if true_classification.loc[i]['effective_group'] == 3:
                g3s.append(i+'.report.csv')
                if i + '.report.csv' not in list(d.columns):   #if patient not in heatmap, add it.
                    d[i + '.report.csv'] = [0 for i in range(len(d.index))]


    # ADD VALUES FOR RB method
    for i in input_rbm_data.index:
        d[input_rbm_data.iloc[i]['patient']][input_rbm_data.iloc[i]['microbe']] = 1
        orgs = metadata[metadata['RNAfilename'] == input_rbm_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]
        #print(orgs)
        if(type(orgs) == type('string')):
            if input_rbm_data.iloc[i]['microbe_genus'] in orgs.split(','):
                yay_v += 1
                d[input_rbm_data.iloc[i]['patient']][input_rbm_data.iloc[i]['microbe']] = d[input_rbm_data.iloc[i]['patient']][input_rbm_data.iloc[i]['microbe']] + 2

    # ADD VALUES FOR LR method
    for i in input_lr_data.index:

        try:
            d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] = d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] + 4   #try to just add to existing values
        except:
            #NEED TO ADD THE ROW
            new_index = [ v for v in d.index] + [input_lr_data.iloc[i]['microbe']]
            d.loc[len(d)]=[0 for i in range(len(d.columns))]
            d.index = new_index
            d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] = 4   # if the value has not existed before, add to the dataframe

        orgs = metadata[metadata['RNAfilename'] == input_lr_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]   # loop through all the organisms
        if(type(orgs) == type('string')):
            if input_lr_data.iloc[i]['microbe_genus'] in orgs.split(','):
                yay_v += 1
                if d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] not in [3,7]:  # if the organism is confirmed, but has not previously been clinically confirmed then...add clinical to it
                    d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] = d[input_lr_data.iloc[i]['patient']][input_lr_data.iloc[i]['microbe']] + 2



    # LOOP THROUGH AND ADD CLINICAL GENERA
    for i in input_lr_data.index:
        orgs = metadata[metadata['RNAfilename'] == input_lr_data.iloc[i]['patient'].split('.')[0]]['organism'].iloc[0]
        if(type(orgs) == type('string')):
            #print(input_lr_data.iloc[i]['patient'])
            #print(orgs)
            for o in orgs.split(','):
                list_of_species = list(set(input_lr_data.loc[input_lr_data['microbe_genus']==o]['microbe']).union(set(input_lr_data.loc[input_lr_data['microbe_genus']==o]['microbe'])))
                a_species_was_found = False
                for s in list_of_species:
                    if not a_species_was_found:
                        try:
                            if d[input_lr_data.iloc[i]['patient']][s] > 0:
                                a_species_was_found = True
                        except:
                            continue

                #print(a_species_was_found)

                if not a_species_was_found:
                    if sum([int(o == v ) for v in  d.index]) > 0:
                        d[input_lr_data.iloc[i]['patient']][o] = 2
                    else:
                        new_index = [ v for v in d.index] + [o]
                        d.loc[len(d)]=[0 for i in range(len(d.columns))]

                        d.index = new_index
                        d[input_lr_data.iloc[i]['patient']][o] = 2

    gcombos = g1s_train + g1s_test + g4s_train +  g4s_test + g2s + g3s

    d = d[gcombos] # this sorts the columns in order of [group1] [group4]


    # Apply some collapsing for the indices
    if(collapse):
        strep_viridans_group = ['Streptococcus pseudopneumoniae ( 257758 )','Streptococcus sp. VT 162 ( 1419814 )',
                                'Streptococcus salivarius ( 1304 )','Streptococcus mitis ( 28037 )',
                                'Streptococcus pseudopneumoniae ( 257758 )','Streptococcus anginosus ( 1328 )']
        strep_preprocessed_names = []
        for i in d.index:
            if(i in strep_viridans_group):
                strep_preprocessed_names.append('Streptococci viridans Group ( XXX )')
            elif('coronavirus' in i):
                strep_preprocessed_names.append('Human coronavirus species ( XXX )')
            elif("Rhino" in i):
                strep_preprocessed_names.append(' '.join(i.split(' ')[0:2]) + " virus ( XXX )")
                #print(strep_preprocessed_names)
            else:
                strep_preprocessed_names.append(i)  # strep pneumo and strep genera (clinical) will be merged

        d.index = strep_preprocessed_names
        first_values = [' '.join(i.split(' ')[0:-4]) if i.count(' ') > 3 else ' '.join(i.split(' ')[0:-3]) for i in d.index]
        d.index = first_values


        for j in d.columns:
            for k in set(first_values):
                #print(k)
                l = list(d[j][d.index==k])
                if 2 in l and 3 in l:
                    # this is bad - it will collapse incorrectly and appear the BM and LR matched (score = 5),
                    # when in reality it was a collapse of RB+Clin (3) and Clin+ (3) again - need to modify the cases post-hoc.
                    # print([3] + [0 for i in range(len(l)-1)])
                    d[j][d.index==k] = [3] + [0 for i in range(len(l)-1)]


        d = d.groupby(d.index).sum()

        post_mod_index = []  # MODIFICATIONS TO UPDATE NAME OF STREP AND STAPH TO BE CLEAR
        for i in d.index:
            if(i == 'Streptococcus'):
                post_mod_index.append('Streptococcus pneumoniae')
            elif(i == 'Staphylococcus'):
                post_mod_index.append('Staphylococcus aureus')
            elif(i == 'Human respiratory syncytial'):
                post_mod_index.append('Human RSV')
            elif(i == 'Acinetobacter'):
                post_mod_index.append('Acinetobacter calcoaceticus')
            elif(i == 'Aspergillus'):
                post_mod_index.append('Aspergillus oryzae')
            elif(i == 'Bacteroides'):
                post_mod_index.append('Bacteroides fragilis')
            elif(i == 'Burkholderia'):
                post_mod_index.append('Burkholderia cenocepacia')
            elif(i == 'Citrobacter'):
                post_mod_index.append('Citrobacter koseri')
            elif(i == 'Enterobacter'):
                post_mod_index.append('Enterobacter spp.')
            elif(i == 'Enterococcus'):
                post_mod_index.append('Enterococcus spp.')
            elif(i == 'Escherichia'):
                post_mod_index.append('Escherichia coli')
            elif(i == 'Fusobacterium'):
                post_mod_index.append('Fusobacterium nucleatum')
            elif(i == 'Haemophilus'):
                post_mod_index.append('Haemophilus influenzae')
            elif(i == 'Klebsiella'):
                post_mod_index.append('Klebsiella pneumoniae')
            elif(i == 'Moraxella'):
                post_mod_index.append('Moraxella catarrhalis')
            elif(i == 'Morganella'):
                post_mod_index.append('Morganella morganii')
            elif(i == 'Mycobacterium'):
                post_mod_index.append('Mycobacterium abscessus')
            elif(i == 'Pasteurella'):
                post_mod_index.append('Pasteurella multocida')
            elif(i == 'Prevotella'):
                post_mod_index.append('Prevotella melaninogenica')
            elif(i == 'Pseudomonas'):
                post_mod_index.append('Pseudomonas aeruginosa')
            elif(i == 'Rothia'):
                post_mod_index.append('Rothia mucilaginosa')
            elif(i == 'Serratia'):
                post_mod_index.append('Serratia marcescens')
            elif(i == 'Stenotrophomonas'):
                post_mod_index.append('Stenotrophomonas maltophilia')
            elif(i == 'Histoplasma'):
                post_mod_index.append('Histoplasma capsulatum')
            else:
                post_mod_index.append(i)

        d.index = post_mod_index


    # sort the indices by alphabetical order
    new_d = d.sort_index(axis=0)
    d = new_d


    #print(d.sum(axis=1).sort_values(ascending=False))
    row_order = list(((d>0)*1).sum(axis=1).sort_values(ascending=False).index)  # row order is the binary sum of occurrence in all
    #print(row_order)
    cg1_train_data = d[g1s_train]
    cg1_test_data = d[g1s_test]
    cg4_train_data = d[g4s_train]
    cg4_test_data = d[g4s_test]
    cg2_data = d[g2s]
    cg3_data = d[g3s]

    output_count = 0
    for subset in [cg1_train_data, cg1_test_data, cg4_train_data, cg4_test_data, cg2_data, cg3_data]:

        output_count += 1

        #sort rows and columns
        subset = subset.reindex(row_order)
        subset_bool = (subset > 0)*1
        subset_bool = (subset_bool.transpose()*[i for i in reversed(range(len(row_order)))]).transpose()  #weight it so that top left corner has values
        col_order = list(subset_bool.sum(axis=0).sort_values(ascending=False).index)
        subset = subset[col_order]


        #reset names to remove .report.csv notation
        new_col_names = []
        for i in subset.columns:
            new_col_names.append(i.split('.')[0])
        subset.columns = new_col_names

        new_cols_temp = ['TA-' + str(c.split('-')[1]) for c in subset.columns]  #convert names to TA-XXX notation
        new_columns = [i[0:6] for i in new_cols_temp]  #get rid of TRNA from the one sample that has weird naming convention
        subset.columns = new_columns
        fig_height = (subset.shape[0]/4.8)*1.2
        fig_width = (subset.shape[1]/8)*1.8
        plt.figure(figsize=(fig_width,fig_height))
        colors = ['white','#FFA2A2','lightgrey','#FFA2A2','#FF3232','#C00000','#FF3232','#C00000']  # trying to make monochromatic FigB
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
        cg = sns.clustermap(subset, cmap = cmap, vmin = 0, vmax=7,figsize=(fig_width,fig_height),row_cluster=False,col_cluster=False)
        plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        ax = cg.ax_heatmap  # this is the important part
        ax.hlines([i for i in range(len(subset.index))], *ax.get_xlim(),lw=.4)
        ax.vlines([i for i in range(len(subset.columns))], *ax.get_ylim(),lw=.4)
        plt.savefig('./output/' + output_directory + '/Fig3B_part' + str(output_count) + '.pdf', bbox_inches='tight')

        plt.show()


        plt.figure(figsize=(fig_width,fig_height))
        colors = ['white','#FFA2A2','#4F4F4F','#FFA2A2','#FF3232','#C00000','#FF3232','#C00000']  # trying to make monochromatic FigB
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
        subset.replace([3,6,7],2,inplace=True)
        subset.replace([1,4,5,8,9],0,inplace=True)
        cg = sns.clustermap(subset, cmap = cmap, vmin = 0, vmax=7,figsize=(fig_width,fig_height),row_cluster=False,col_cluster=False)
        plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        ax = cg.ax_heatmap  # this is the important part
        ax.hlines([i for i in range(len(subset.index))], *ax.get_xlim(),lw=.4)
        ax.vlines([i for i in range(len(subset.columns))], *ax.get_ylim(),lw=.4)
        plt.savefig('./output/' + output_directory + '/Fig3A_part' + str(output_count) + '.pdf' , bbox_inches='tight')

        plt.show()

    #reset names to remove .report.csv notation
    new_col_names = []
    for i in d.columns:
        new_col_names.append(i.split('.')[0])
    d.columns = new_col_names

    new_cols_temp = ['TA-' + str(c.split('-')[1]) for c in d.columns]  #convert names to TA-XXX notation
    new_columns = [i[0:6] for i in new_cols_temp]  #get rid of TRNA from the one sample that has weird naming convention
    d.columns = new_columns


    fig_height = (d.shape[0]/4.8)*2
    fig_width = (d.shape[1]/6.4)*2

    if print_output:
        print(d.head())

    return d
