import pandas as pd
import sys
import argparse
import MySQLdb
import gzip


parser = argparse.ArgumentParser(description = 'Remove all sequences that align to taxa under vertebrata')
parser.add_argument('-i','-infile', help='input retained-NR.retained-NT file',type=str)
parser.add_argument('-R1','-R1', help='Read 1 .fastq.gz file', type=str)
parser.add_argument('-R2','-R2', help='Read 2 .fastq.gz file', type=str)
parser.add_argument('-o','-outfile', help='output filepath',type=str)
args = parser.parse_args()


tsdb = MySQLdb.connect(host='localhost',user='mysql_user',passwd='balamuthia',db='taxa_scoring')
ncbi_db = MySQLdb.connect(host='localhost',user='mysql_user',passwd='balamuthia',db='NCBI_Taxonomy')

# read in the .fastq R1 file

#with gzip.open(args.R1, 'rb') as f:
#    file_content = f.read()


print("Step 1. importing retained file")

count_vertebrate = 0
keep_seqids = []
r = open(args.i, 'r')
print(args.i)
for line in r:
    l = line.strip()
    if l[0] == '>':
        taxids = [l.split(':')[1], l.split(':')[3]]
        is_vertebrate = False
        for t in taxids:
            Vpath2root = 'call vedas_path_to_root_node(' + str(t) + ')'
            Vp2r = pd.read_sql(Vpath2root,ncbi_db)['@path_to_root'][0].split(';')
            if '7742' in Vp2r:
                is_vertebrate = True
                count_vertebrate += 1

        if not is_vertebrate:
            keep_seqids.append('@'+(':'.join(l.split(':')[4:]).split('/')[0]))

print("Found " + str(count_vertebrate) + " sequences mapping to vertebrate.")
print("Keeping " + str(len(keep_seqids)) + ", " + str(float(len(keep_seqids))/float((len(keep_seqids) + count_vertebrate))*100) + "% of total sequenes.")




def parse_fastq(filename, keep):

    f = gzip.open(filename, 'rb')
    count = 0
    dict_df = {}
    seqid=''
    seq = ''
    line2 = ''
    qual = ''
    seqid_map = ''

    keep_this_seq = False

    new_file = ''

    for line in f.readlines():
        l = line.strip()
        if count % 4 == 0:
            if count > 3 and keep_this_seq:
                #dict_df[seqid_map] = {'seqid':seqid, 'qual':qual,'line2':line2,'seq':seq}
            	new_file = new_file + ('\n'.join([seqid,seq,line2,qual]))
            	keep_this_seq = False

            if l.split(' ')[0].replace('@',''):
            	keep_this_seq = True
	            seqid = l
	            #seqid_map = seqid.split(' ')[0]
        elif count % 4 == 1 and keep_this_seq:
            seq = l
        elif count % 4 == 2 and keep_this_seq:
            line2 = l
        elif count % 4 == 3 and keep_this_seq:
            qual = l

        count += 1


    f.close()

    #df = pd.DataFrame.from_dict(dict_df, orient='index')
    return new_file

print("Step 2. Parsing R1 .fastq file")
dfR1 = parse_fastq(args.R1, keep_seqids)

print("Step 3. Writing R1 output")
f = open(args.o+'_R1.fastq','w')
f.write(dfR1)
f.close()

dfR1 = None

print("Step 4. Parsing R2 .fastq file")
dfR2 = parse_fastq(args.R2, keep_seqids)

print("Step 5. Writing R2 output")
f = open(args.o+'_R2.fastq','w')
f.write(dfR2)
f.close()

dfR2 = None

'''
keep_set = list(set(keep_seqids).intersection(set(dfR1.index)))
R1output = dfR1.loc[keep_set]

keep_set = list(set(keep_seqids).intersection(set(dfR2.index)))
R2output = dfR2.loc[keep_set]

writeR1 = []
for index, row in R1output.iterrows():
    writeR1 = writeR1 + [row["seqid"],row["seq"],row["line2"],row["qual"]]

f = open(args.o+'_R1.fastq','w')
f.write('\n'.join(writeR1))
f.close()

writeR2 = []
for index, row in R2output.iterrows():
    writeR2 = writeR2 + [row["seqid"],row["seq"],row["line2"],row["qual"]]

f = open(args.o+'_R2.fastq','w')
f.write('\n'.join(writeR2))
f.close()
'''
