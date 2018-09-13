# Host-MicrobeLRTI

This directory contains all the required data and scripts required to run the analysis outlined in Integrating Host Response and Unbiased Microbe Detection for Lower Respiratory Tract Infection Diagnosis in Critically Ill Adults. The Pathogen v. Commensal model is implemented in a jupyter notebook, while the host and microbe combined model is implemented in R markdown, with all data pre-processing steps pre-set in the mBALPkg.

## Citation

# Pathogen v. Commensal Model

To differentiate pathogens from respiratory commensals, we developed rules-based and logistic regression models (RBM, LRM) that are implemented in the jupyter notebook: **mBAL_PathogenvCommensal.ipynb**. 

### Input Files

Running the jupyter notebook requires the following input files: 

#### Metadata

The input metadata MUST contain the following fields: 
+ **StudyId:** This is unique sample identifier; If using the host functionality downstream, the StudyId in the metadata should match with the column names in the genecounts matrix.
+ **effective_group:** This indicates which "Infection Group" the sample is from - for the purposes of the mBAL study, Group 1 = LRTI+C+M, Group 4 = No-LRTI, Group 2 = LRTI+C, and Group 3 = LRTI-UNK. If using different input, the groups should be coded such that "1" indicates the *positive* cases and "4" indicates *negative* cases.
+ **organism:** The genus-level name of the organism identified in standard clinical microbiology for true positive cases. If more than one pathogen was identified, then this should be a comma separated list containing no spaces between organisms.
+ **DNAfilename:** The base name of the DNA-seq microbe report file.
+ **RNAfilename:** The base name of the RNA-seq microbe report file.

| MBALStudyId       | effective_group | organism  | DNAfilename | RNAfilename|
| :-----: |:---:|:----------:|:--------:|:--------:
| 205      | 1 | Enterovirus ( 12059 ) | mBAL-205-DNA-TA1-B8 | mBAL-205-RNA-TA1-QIA-61917 |
| 212      | 1      |   Escherichia ( 561 ),Klebsiella ( 570 ) | mBAL-213-DNA-TA1-B8 | mBAL-213-RNA-TA1-QIA-62317 |
| 215 | 4      |  | mBAL-215-DNA-TA1-B10 | mBAL-215-RNA-TA1-B10 |
| ... | ... | ... | ... | ... |

**Note:** The metadata file used in the published analysis contained several additional columns that were used for analysis of clinical covariates and may be referenced in the **mBAL_LRTIvNoLRTI_HostPathogenCombined.Rmd** script, but are not integral to replicating the main functionality.

#### Microbe report files
The **mBAL_PathogenvCommensal.ipynb** analysis requires report files output after alignment to NCBI NR and NT databases. The current version supports report files downloaded from the DeRisi lab pipeline. Future versions may support report files downloaded from IDSeq.

#### .reference/
The **mBAL_PathogenvCommensal.ipynb** analysis references two input files. 
1. **viruses.txt** contains a list of all viruses identified in the study - this list is used to generate the boolean variable "is_virus" that is used in the logistic regression model.
2. **known_respiratory_pathogens.txt** contains a list of known pathogens curated a priori based on reference studies of pathogen prevalence in ICU settings.

#### Merged Microbial file
The Host Pathogen Combined model implemented in **mBAL_LRTIvNoLRTI_HostPathogenCombined.Rmd** requires a merged_genusrpm.tsv file which contains a file merging all the genera from all samples. This is used for calculation of diversity metrics.


## Running the Pathogen v. Commensal Model on Published Study Input Files

1. Download the github directory including mBAL_PathogenvCommensal.ipynb script and the data, scripts, and reference sub-directories in the same directory structure as in this github repo. 
2. Create an output directory (./data/XXX) and modify the output variable accordingly - this will dictate where all .pdf and associated files are written to.
3. Run the notebook. 


## Overriding the Pathogen v. Commensal Model for use on Other Datasets

**At this point, the Pathogen v. Commensal Model is only set up to run with output report files from the DeRisi lab pathogen discovery pipeline. This limits the functionality greatly.** Integration with IDSeq pipeline output files is coming soon.

1. Download the .ipynb file into a directory: analysis_directory
2. Create a new data directory: analysis_directory/data/my_data_dir
3. Add the "metadata" file to the data directory
4. Create a "background model" file (ie background_models.txt) containing the integer corresponding the the background model for which to generate data.
```
4
```

5. Run the automated_pathogen_pipeline.sh script to download the report files with the corresponding background model:
```
bash ./scripts/automated_pathogen_pipeline.sh [/path/to/my_data_dir] [metadata filename (without path)] [scripts directory (ie ./scripts)]
```
This will create a sub-directory BM_X (ie BM_4 for background model 4) containing all RNA-seq and DNA-seq report files associated with the filenames in your metadata `RNAfilename` and `DNAfilename` fields.

6. Move the RNA files to an RNA sub-directory and DNA files to DNA sub-directory:
```
cd /path/to/my_data_dir
mkdir BM_4/RNA
for i in `ls BM_4/*RNA*`; do mv $i BM_4/RNA/; done;
mkdir BM_4/DNA
for i in `ls BM_4/*DNA*`; do mv $i BM_4/DNA/; done;
```

7. Set the input variables in the **mBAL_PathogenvCommensal.ipynb** notebook
```
data_directory =  XXX
output_directory = XXX
metadata = XXX
```

8. Run the notebook.




# Combined Host and Microbe

To identify LRTI-positive patients and differentiate them from controls with non-infectious acute respiratory illnesses, we developed pathogen, microbiome diversity, and host gene expression metrics. Here, I provide scripts to 
1. Replicate the full analysis using the published datasets
2. Apply the host gene expression metric to other transcriptomic data

## 1. Running the LRTI v. No LRTI Host Pathogen Combined Script

The host pathogen combined analysis is impelented in the R package mBALPkg. It can be installed via this command:


``` devtools::install_github("katrinakalantar/miniBAL_study", subdir="mBALPkg") ```


Once installed, the library can be loaded using command: ``` library(mBALPkg) ```

This package includes the processed host expression data and associated microbe data from the Pathogen v. Commensal model. The pre-processing steps can be viewed in the **mBALPkg/data-raw/preprocess_data.R** script. Once loaded, you can modify the "output" variable in the R markdown script **mBAL_LRTIvNoLRTI_HostPathogenCombined.Rmd** to ensure that the output files (.pdf, .csv, etc.) are written to an appropriate directory. Finally, run the script.


## 2. Run the Host classifier on your new gene counts

```
Rscript run_host.R -c counts_matrix.tsv -m host_model.txt [-o output.txt]
```

#### Input:

**counts_matrix.tsv:** File path to a gene counts matrix with ENSEMBL gene IDs as row names and unique sample names in the header. Each column corresponds to the transcriptional profile of a single sample.

**host_model.txt:** File path to a tab-delimited matrix containing gene IDs in the first column and model weights in the second column.

For example, the **./reference/mBAL_HostClassifier.txt** contains the following:
```
geneID  weight
ENSG00000117523 1
ENSG00000058673 1
ENSG00000134851 1
ENSG00000132424 1
ENSG00000135218 -1
ENSG00000092964 -1
ENSG00000106991 -1
ENSG00000107223 -1
ENSG00000157106 1
ENSG00000102908 1
ENSG00000104979 -1
ENSG00000090013 -1
```

**output.txt:** *Optional* output filename indicating where to write the host classifier predicted scores. If no filepath is indicated, default is "./output.txt"



It is possible to replicate the results for the test set in the paper and .Rmd file using this script via the following command:

```
Rscript run_host.R -c mBALPkg/data-raw/NEB-TA-ONLY-genecounts-03.29.18.csv -m reference/mBAL_HostClassifier.txt -t reference/TESTnames.txt
```
