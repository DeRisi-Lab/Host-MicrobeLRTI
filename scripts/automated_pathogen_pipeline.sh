#!/bin/bash

# Katrina Kalantar
# December 10, 2017
# Automate the downloading of all pathogen files associated with a particular metadata
# Input: Requires a data directory set up containing the metadata file with headers "RNAfilename" and "DNAfilename"
# Output: Will populate the data directory with BM_X background-model-specific downloads of the pathogen data,
#         combined genus- and species-level rpM files, and a filenames.txt file

# $1 is the data directory name containing the metadata file
# $2 is the metadata filename
# $3 is the scripts directory

# move into the data directory to execute all commands
cd $1

# run get_all_filenames to extract associated filenames from the metadata file and output to filenames.txt
bash $3get_all_filenames.sh $1 $2
echo "completed getting filenames"

# run generate_reports_by_background_model.v2.sh to use CURL command implemented by Joe to download all reports
# associated with the filenames found in filenames.txt
source $3generate_reports_by_background_model.v2.sh background_models.txt filenames.txt
echo "completed getting files"

# run python script to merge the NT rpM per microbe genus / species
python $3create_merged_genus_rpm_v2.py BM_4 --genus
python $3create_merged_genus_rpm_v2.py BM_4 --species

# move files into RNA/ and DNA/ -specific directories for downstream pathogen analysis
for i in `cat background_models.txt`; do mkdir BM_$i/RNA; mv BM_$i/*RNA*report* BM_$i/RNA; mkdir BM_$i/DNA; mv BM_$i/*DNA*report* BM_$i/DNA; done;
echo "COMPLETED SCRIPT"
