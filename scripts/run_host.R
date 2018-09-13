
# run as Rscript run_host.R -c [counts matrix] -h [host model]

library("optparse")

option_list = list(
  make_option(c("-c", "--countsmatrix"), type="character", default=NULL, 
              help="Matrix of genecounts with samples in columns and gene names (ENSEMBL IDs) in row", metavar="character"),
  make_option(c("-m", "--hostmodel"), type="character", default = NULL,
              help="Tab delimited file containing the host model with ENSEMBL Gene ID, Weight", metavar="character"),
  make_option(c("-o", "--output"), type="character", default = "output.txt",
              help="Tab delimited file containing per-sample host metric score ( Per mBAL study, score > -4 indicates \"Infection\", score <=-4 indicates \" No Infection \")", 
              metavar="character"),
  make_option(c("-t", "--testsamples"), type="character", 
              help="File containing the sample IDs to be tested, if only a subset of the columns of the genecounts matrix", 
              metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$countsmatrix)){
  print_help(opt_parser)
  stop("Input counts matrix must be supplied.\n", call.=FALSE)
}
if (is.null(opt$hostmodel)){
  print_help(opt_parser)
  stop("Input host model must be supplied.\n", call.=FALSE)
}
if (is.null(opt$output)){
  print("No output file specified; find results at output.txt")
}


# Load dependencies
print("Loading dependencies: mBALPkg, DESeq2")
suppressMessages(library(mBALPkg))
suppressMessages(library(DESeq2))

########### Main program... ###########

# read in the raw genecounts file:
genecounts = read.csv(opt$countsmatrix, header=TRUE, row.names = 1)

print(dim(genecounts))

if(!is.null(opt$testsamples)){
  ids = read.table(opt$testsamples)
  print(dim(genecounts))
  genecounts = genecounts[,colnames(genecounts) %in% ids$V1]
  print(dim(genecounts))
}


# create dummy metadata file to enable:
#   1. generation of ExpressionSet Object underlyng all mBALPkg functions
#   2. normalization via DESeq
#metadata = as.data.frame(cbind("studyID" = colnames(raw_genecounts), 
#                 "effective_group" = as.factor(c(rep(1, ceiling(ncol(raw_genecounts)/2)),rep(2, floor(ncol(raw_genecounts)/2))))))
metadata = as.data.frame(cbind("studyID" = colnames(genecounts), 
                 "effective_group" = as.factor(c(rep(1, ceiling(ncol(genecounts)/2)),rep(2, floor(ncol(genecounts)/2))))))
rownames(metadata) <- metadata$studyID

print(dim(metadata))
print(rownames(metadata))

# remove non-protein-coding genes and apply filters...
# ...hard-coded values match those used in the mBAL analysis
# NOTE, in the core analysis the TEST_eset samples were isolated for the full dataset. 
# Differences in which samples are included prior to normalization may impact final output results
filtered_eset <-
  basic.setup(
    genecounts,
    metadata,
    genemap,
    removePC = TRUE,
    read_filter = 10000,
    nonzero_filter = .5,
    equalize = FALSE
  )

print(dim(filtered_eset))

print("created filtered eset")

# Run DESeq normalization on the input genecounts matrix using the same 
# parameters that were used for mBAL study TEST_eset
test_expressionset <- format_host_features(
  traditionalclassification.setup(
    filtered_eset,
    normalization_method = "DESeq",
    rld = FALSE
  ),
  zscore = FALSE
)


# Read in the model file (containing ENSEMBL gene names and the associated weights)
model <- read.csv(opt$hostmodel, header=TRUE, row.names=1, sep='\t')
print("Using input host model with weights: ")
print(model)
host_model <- as.numeric(model$weight)
names(host_model) <- rownames(model)


# Print out a matrix indicating which genes are present in your dataset...
# ...just in case some genes were filtered out due to low prevalence in input dataset

present <- c(names(host_model) %in% rownames(exprs(test_expressionset)))
if(sum(present) == length(names(host_model))){
	print("All Host Model genes are present in this dataset.")	
}else{
	print("WARNING: Not all host model genes are present in this dataset.\n The following genes are missing:")
	names(present) <- names(host_model)
	present <- as.matrix(present)
	colnames(present) <- c("gene_present_in_dataset")
	print(present[present$gene_present_in_dataset == FALSE,])
}


# This function uses the default center and scale values from the mBAL Training dataset. 
# If you want to use your own values...
#   1. Create a new input argument for a tab-delimited file with columns A) Center values and B) Scale values
#   2. Read in the file and parse the center and scale values as lists
#   3. Provide the lists to the function as center_params = A, and scale_params = B
classification_result <- classify_using_host_transcriptome.withmodel(test_expressionset, host_model)


# Print the result to console
print("Classification_result: ")
print(as.matrix(classification_result$test_classification))


# Write the results to output
write.table(as.matrix(classification_result$test_classification), file=opt$output, quote = FALSE)


