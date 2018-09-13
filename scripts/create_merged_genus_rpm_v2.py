import pandas as pd
import glob
import os
import sys

#python create_merged_genus_rpm_v2.py [directory] [--genus | --species]

keep_kingdoms = ['Bacteria','Fungi','Viruses','Escherichia coli']  #keep any genus in these kingdoms
#layout of file header = 'Category,Species,Genus,Score1,Score2,Score3,NT Species Z,NT Genus Z,NT Species rM,NT Genus rM,NR Species Z,NR Genus Z,NR Species rM,NR Genus rM'


directory = sys.argv[1]
os.chdir(directory)

csv_files = glob.glob('*.csv') #('*.csv')   #get list of all .csv files in the directory

full_dict = {}
for c in csv_files:

    '''
    file = open(c).readlines()
    f = [line.strip().split(',') for line in file]

    filename = c.replace('.report.csv','')
    full_dict[filename] = {}



    for l in f:
        print(len(l))
        if len(l) == 14 + 1 and 'Category' not in l: #this is a valid non-header line
            if l[0+ 1] in keep_kingdoms:            #keep this info
                if(sys.argv[2] == '--genus'):
                    genus = l[2+ 1]
                else:
                    genus = l[1+1]  #THIS IS ACTUALLY THE SPECIES NAME
                genus_rM = l[9+ 1]
                full_dict[filename][genus] = genus_rM
    '''

    filename = c.replace('.report.csv','')
    full_dict[filename] = {}

    try:
    	df = pd.read_csv(c, error_bad_lines=False)
    	df = df[df['Category'].isin(keep_kingdoms)]
    	df = df[['Category','Genus','NT Genus rM','NT Genus Z','Species', 'NT Species rM']]
        
    	#group by Genus, but keep the species with the greatest rM value
    	idx = df.groupby(['Genus'])['NT Species rM'].transform(max) == df['NT Species rM']
    	df = df[['Category','Genus','NT Genus rM','NT Genus Z','Species']][idx] #HOPEFULLY THIS KEEPS THE MAX SPECIES IN COLLAPSING GENUS

    	for l in df.index:
            if(sys.argv[2] == '--genus'):
                genus = df.loc[l]['Genus']
            else:
                genus = df.loc[l]['Species']
            genus_rM = df.loc[l]['NT Genus rM']
            full_dict[filename][genus] = genus_rM

    except:
        print("FAILED to add sample to merged file:" + str(c))

result_DF = pd.DataFrame(full_dict)
result_DF.fillna(value=0, inplace=True) #convert NA values -> 0
print(result_DF)

if sys.argv[2] == '--genus':
    result_DF.to_csv('merged_genusrpm.tsv',sep='\t')  #write to file
else:
    result_DF.to_csv('merged_genus_withspeciesname_rpm.tsv',sep='\t')  #write to file
