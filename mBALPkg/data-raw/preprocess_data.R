source("../R/utility.R")  # run this SOURCE command only when generating the data files for the package
library(vegan)
library(DESeq2)
library(biomaRt)

# some steps require that the gene counts matrix have HGNC gene identifiers
# instead of ENSEMBL, so I generated genemap.rda to contain all the mappings to 
# do the conversion without having to access biomaRt

# Ran this code once to generate the genemap.rda file, which was loaded 
# listMarts(host = "www.ensembl.org")
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                dataset = "hsapiens_gene_ensembl",
#                host = "www.ensembl.org")
# genemap <-
#   getBM(
#     attributes = c(
#       "ensembl_gene_id",
#       "entrezgene",
#       "hgnc_symbol",
#       "gene_biotype",
#       "description"
#     ),
#     filters = "ensembl_gene_id",
#     values = genenames,
#     mart
#   )

load("genemap.rda")  # now just load the existing map from Ensembl -> HGNC
devtools::use_data(genemap, overwrite = TRUE)


####### Raw Data #######

#
# Read in raw data
#

# Read in the raw data from gene counts .csv file
raw_genecounts <-
  read.csv('NEB-TA-ONLY-genecounts-03.29.18.csv',
           header = TRUE,
           row.names = 1)
devtools::use_data(raw_genecounts, overwrite = TRUE)  # save to raw_genecounts

# Read in raw metadata from .csv file
raw_metadata <-
  read.csv('NEB-TA-ONLY-metadata-03.29.18withFilenames.KKadaptation_v3.csv',
           row.names = 2)
devtools::use_data(raw_metadata, overwrite = TRUE)  # save to raw_metadata




####### Format Metadata #######

#
# Format microbe and diversity data to be appended as columns to metadata
#

# Read in microbe scores (output from LR model implemented in python)
# Format: sampleID / microbe_score (top LR probability per patient) / effective_group (patient group by LRTI status)
microbe_scores <-
  format_microbe_features3(read.csv('patient_microbe_scores.csv', row.names = 1), raw_metadata)

# merge microbe data with metadata matrix
metadata_with_microbe <- merge(raw_metadata, microbe_scores, by = 0, all = TRUE)
# clean up row names after merging
rownames(metadata_with_microbe) <- metadata_with_microbe$Row.names 
# remove Row.names column
metadata_with_microbe <- metadata_with_microbe[, !(colnames(metadata_with_microbe) %in% c("Row.names"))]  


# import the table of microbes x patient - 
# this matrix was generated separately using .report.csv files from /data/BM_4 and this script:
# python create_merged_genus_rpm_v2.py /data/BM_4 --genus
merged_microbes <-
  read.csv(
    'merged_genusrpm.tsv',
    header = TRUE,
    sep = '\t',
    row.names = 1
  )   
devtools::use_data(merged_microbes, overwrite = TRUE)

# compute the Shannon's Diversity Index
# returns one diversity value per filename (RNA and DNA)
shannon_div <-
  vegan::diversity(
    t(merged_microbes),
    index = "shannon",
    MARGIN = 1,
    base = exp(1)
  )

# format diversity features 
# DNAfilename / RNA filename / effective_group / DNAdiv / RNAdiv
diversity_data <-
  format_diversity_features(shannon_div, raw_metadata)

# merge diversity data with metadata+microbe matrix
metadata_with_microbe_with_diversity <-
  merge(metadata_with_microbe,
        diversity_data,
        by = 0,
        all = TRUE)
# clean up row names after merging
rownames(metadata_with_microbe_with_diversity) <- metadata_with_microbe_with_diversity$Row.names  
# remove Row.names column
metadata_with_microbe_with_diversity <-
  metadata_with_microbe_with_diversity[, !(colnames(metadata_with_microbe_with_diversity) %in% c("Row.names"))]  


# update the names to avoid duplicates
metadata_with_microbe_with_diversity <-
  metadata_with_microbe_with_diversity[,!(
    colnames(metadata_with_microbe_with_diversity) %in% c(
      "effective_group.x.1",
      "effective_group.y",
      "effective_group.y.1"
    )
  )]
names(metadata_with_microbe_with_diversity)[names(metadata_with_microbe_with_diversity) == 'effective_group.x'] <-
  'effective_group'
metadata_with_microbe_with_diversity <-
  metadata_with_microbe_with_diversity[,!(
    colnames(metadata_with_microbe_with_diversity) %in% c("RNAfilename.y", "DNAfilename.y")
  )]
names(metadata_with_microbe_with_diversity)[names(metadata_with_microbe_with_diversity) == 'RNAfilename.x'] <-
  'RNAfilename'
names(metadata_with_microbe_with_diversity)[names(metadata_with_microbe_with_diversity) == 'DNAfilename.x'] <-
  'DNAfilename'

devtools::use_data(metadata_with_microbe_with_diversity, overwrite = TRUE)  # save to metadata_with_microbe_diversity





# generate an ExpresionSet object from the raw data (to be used in diversity analysis, where scaling cannot be included)
raw_eset <-
  basic.setup(
    raw_genecounts,
    metadata_with_microbe_with_diversity,
    genemap,
    read_filter = 0,
    nonzero_filter = .5,
    equalize = FALSE
  )
devtools::use_data(raw_eset, overwrite = TRUE)  # save to raw_eset

# filter the full transcriptome data matrix based on minimum number of reads, protein-coding, and non-zero values
filtered_eset <-
  basic.setup(
    raw_genecounts,
    metadata_with_microbe_with_diversity,
    #_with_microbevectors
    genemap,
    read_filter = 90000, #previously 10,000  but 90,000 makes no difference
    nonzero_filter = .5,
    equalize = FALSE
  )



#
# Estimate cell-type proportions
#

library("EpiDISH")

# CIBERsORT requires that the gene counts matrix have HGNC gene identifiers
# instead of ENSEMBL, so convert to HGNC notation

genenames <- rownames(exprs(filtered_eset))

# create an index in which you are rearranging ensmbl gene ids from genemap into same order as they are in genenames
idx <- match(genenames, genemap$ensembl_gene_id)  
hgnc_symbol <- genemap$hgnc_symbol[idx]

# convert original genecounts matrix, filtered for PC-genes, to have HGNC rownames (necessary to run CIBERSORT)
pc_hgnc_genecounts = exprs(filtered_eset)
rownames(pc_hgnc_genecounts) = hgnc_symbol

# load the LM22 matrix for CIBERSORT
LM22 <-
  read.csv(
    'LM22.txt',
    sep = '\t',
    row.names = 1
  )

# run CIBERSORT to estimate proportions
epidish_cbs <-
  epidish(as.matrix(pc_hgnc_genecounts), as.matrix(LM22), method = "CBS")
cell_proportions <- epidish_cbs$estF

devtools::use_data(cell_proportions, overwrite = TRUE)  # save cell type proportion for each patient in cell_proportions

# append cell type proportions to existing metadata
new_filtered_eset_pData <-
  cbind(pData(filtered_eset), cell_proportions)
pData(filtered_eset) <- new_filtered_eset_pData

devtools::use_data(filtered_eset, overwrite = TRUE)   # save filtered_eset with cell type proportions to filtered_eset




####### Format Gene Counts #######

#
# Create training and test datasets
#


# select group 1 and group 4 patients (known)
classification_input <-
  filtered_eset[, filtered_eset$effective_group %in% c(1, 4)]

# select group 2 and group 3 patients (unknown)
unknown_input <-
  filtered_eset[, filtered_eset$effective_group %in% c(2, 3)]

# format the classification data (normalize, log-transform, scale phenotype data)
classification_eset <-
  format_host_features(
    traditionalclassification.setup(
      classification_input,
      normalization_method = "DESeq",
      rld = FALSE
    ),
    zscore = TRUE
  )
pData(classification_eset) <- smart_scale(pData(classification_eset))
devtools::use_data(classification_eset, overwrite = TRUE)   # save the classification data

# format the classification data (normalize, log-transform, NO scaling on phenotype data)
classification_eset_notscaled <-
  format_host_features(
    traditionalclassification.setup(
      classification_input,
      normalization_method = "DESeq",
      rld = FALSE
    ),
    zscore = FALSE
  )
pData(classification_eset_notscaled) <- smart_scale(pData(classification_eset_notscaled))
devtools::use_data(classification_eset_notscaled, overwrite = TRUE)   # save the non-scaled classification data


# format the unknown data using the scaling and centering values from the test set
values_for_unknowns <- format_host_features_get_z_values(
  traditionalclassification.setup(
    classification_input,
    normalization_method = "DESeq",
    rld = FALSE
  ),
  zscore = TRUE
)
formatted_unknowns <- traditionalclassification.setup(unknown_input,
                                                      normalization_method = "DESeq",
                                                      rld = FALSE)
ft <- t(exprs(formatted_unknowns))
cft <- sweep(ft, 2, values_for_unknowns$centering_values)
sft <- sweep(cft, 2, values_for_unknowns$scaling_values, '/')
unknown_eset <-
  ExpressionSet(t(sft),  phenoData = AnnotatedDataFrame(pData(unknown_input)))
pData(unknown_eset) <- smart_scale(pData(unknown_eset))
devtools::use_data(unknown_eset, overwrite = TRUE)   # save the test dataset


# create a full dataset that is normalized within itself across all groups
fulldata_eset <-
  format_host_features(
    traditionalclassification.setup(
      filtered_eset,
      normalization_method = "DESeq",
      rld = FALSE
    ),
    zscore = TRUE
  )
devtools::use_data(fulldata_eset, overwrite = TRUE)  # save the full normalized, scaled dataset to fulldata_eset variable


### 12/15/17 CREATE A TRAINING AND TEST SET
TRAINING_NAMES <-
  c(
    "TA.212",
    "TA.225",
    "TA.298",
    "TA.304",
    "TA.314",
    "TA.315",
    "TA.335",
    "TA.337",
    "TA.343",
    "TA.350",
    "TA.349",
    "TA.273",
    "TA.331",
    "TA.221",
    "TA.220",
    "TA.215",
    "TA.270",
    "TA.241",
    "TA.211",
    "TA.218"
  )

# do all normalization, scaling, and  centering for the training data in isolation from the test data
TRAINING_eset <-
  format_host_features(
    traditionalclassification.setup(
      basic.setup(
        raw_genecounts[, TRAINING_NAMES],
        metadata_with_microbe_with_diversity,
        genemap,
        read_filter = 0,
        nonzero_filter = .5,
        equalize = FALSE
      ),
      normalization_method = "DESeq",
      rld = FALSE
    ),
    zscore = FALSE
  )

devtools::use_data(TRAINING_eset,
                   overwrite = TRUE)

# get the scaling values from the training data to be used to scale the test data
scaling_values_from_training_data <- format_host_features_get_z_values(
  traditionalclassification.setup(
    basic.setup(
      raw_genecounts[, TRAINING_NAMES],
      metadata_with_microbe_with_diversity,
      genemap,
      read_filter = 0,
      nonzero_filter = .5,
      equalize = FALSE
    ),
    normalization_method = "DESeq",
    rld = FALSE
  ),
  zscore = TRUE
)
devtools::use_data(scaling_values_from_training_data,
                   overwrite = TRUE)

# apply normalization, scaling, and  centering for all samples not in the training set (groups 1,2,3,4)
NOT_TRAINING_eset <-
  format_host_features(
    traditionalclassification.setup(
      filtered_eset[, !(colnames(filtered_eset) %in% TRAINING_NAMES)],
      normalization_method = "DESeq",
      rld = FALSE
    ),
    zscore = FALSE
  )

devtools::use_data(NOT_TRAINING_eset,
                   overwrite = TRUE)


# apply normalization, scaling, and  centering for test set (groups 1,4 not in training list)
TEST_NAMES <-
  colnames(classification_eset[, !(colnames(classification_eset) %in% TRAINING_NAMES)])

TEST_eset <- format_host_features(
  traditionalclassification.setup(
    ExpressionSet(
      as.matrix(raw_genecounts[rownames(exprs(TRAINING_eset)), TEST_NAMES]),
      phenoData = AnnotatedDataFrame(metadata_with_microbe_with_diversity[TEST_NAMES, ])
    ),
    normalization_method = "DESeq",
    rld = FALSE
  ),
  zscore = FALSE
)

devtools::use_data(TEST_eset,
                   overwrite = TRUE)


# bind the normalized, scaled training data and the normalized test data
TRAINING_AND_TEST_eset <- format_host_features(ExpressionSet(
  as.matrix(cbind(exprs(TRAINING_eset), exprs(TEST_eset))),
  phenoData = AnnotatedDataFrame(metadata_with_microbe_with_diversity[c(TRAINING_NAMES, TEST_NAMES), ])), zscore = FALSE)
devtools::use_data(TRAINING_AND_TEST_eset,
                   overwrite = TRUE)



#
# Create a static reference file of gene names ENSEMBL to HGNC
#

# create the genemap variable, so we don't rely on biomaRt

# ran this section of code once to generate the HSapiens_gene_ensembl.csv file
# library( "biomaRt" )
# listMarts(host="www.ensembl.org";;)
# mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org";;)
# genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                     filters = "ensembl_gene_id",
#                     values = rownames(genecounts),
#                     mart)
# write.csv(genemap,"HSapiens_gene_ensembl.csv",quote=FALSE,row.names=FALSE)

#genemap <- read.csv("HSapiens_gene_ensembl.csv", header = TRUE)
#devtools::use_data(genemap, overwrite = TRUE)
