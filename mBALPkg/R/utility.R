#
# FUNCTIONS TO BE EXPORTED FOR USE IN PACKAGE
#

#' filterZeroRows
#'
#' Remove rows (genes) that are present in less than X% of samples
#' @param input_matrix matrix; full genes x samples matrix
#' @param is_nonzero_in_threshold float; the proportion of samples that must contain any given gene
#' @keywords filter
#' @export
#' @examples
#' filterZeroRows(input_matrix, is_nonzero_in_threshold)

filterZeroRows <- function(input_matrix, is_nonzero_in_threshold) {
  mabsent <- c()
  for (i in 1:ncol(input_matrix)) {
    absent <- as.numeric(input_matrix[, i]) == 0  # Zero spot = absent
    mabsent <-
      cbind(mabsent, absent)  # Add the sample to an "absent" matrix with all other samples
  }
  mpresent <- mabsent == 0   # Determine which genes are present
  
  # Find which genes are present in at least X% of the samples (is_nonzero_in_threshold param)
  mfilter <-
    apply(mpresent, 1, sum) > ncol(input_matrix) * is_nonzero_in_threshold
  mpresent30 <- input_matrix[mfilter == "TRUE",]
  genecounts <- mpresent30
  return(genecounts)
}


#' as.numeric.factor
#'
#' Convert numeric to factor
#' @param x numeric value
#' @keywords factor
#' @export
#' @examples
#' as.numeric.factor(x)

as.numeric.factor <- function(x) {
  as.numeric(levels(x))[x]
}


#
# FUNCTIONS REQUIRED FOR DATA PROCESSING; not exported for use in mBALAnalysisPkg
#

# Function: Take in the microbe scores and format into a matrix that can be appended to existing metadata
format_microbe_features3 <- function(microbe_scores, metadata) {
  
  # keep only the microbe scores for which RNAfilename is in the metadata
  microbe_scores <-
    microbe_scores[rownames(microbe_scores) %in% metadata$RNAfilename, ]  
  
  # sort
  microbe_scores <- microbe_scores[intersect(metadata$RNAfilename,rownames(microbe_scores)),]
  
  # assign appropriate row names
  rownames(microbe_scores) <-
    rownames(metadata[metadata$RNAfilename %in% rownames(microbe_scores), ])  
  
  # add effective_group variable
  microbe_scores$effective_group <-
    metadata[rownames(metadata) %in% rownames(microbe_scores), ]$effective_group  
  
  # keep only the following variables of interest
  microbe_model <-
    microbe_scores[, colnames(microbe_scores) %in% c(
      "X0","effective_group")]
  colnames(microbe_model) <- c("microbe_score","effective_group")
  
  return(microbe_model) 
}



# Function: Create matrix of diversity features for classification analysis.
format_diversity_features <- function(simpsons_div, metadata) {
  
  # select a subset of the metadata
  diversity_matrix <-
    metadata[, c("DNAfilename", "RNAfilename", "effective_group")]
  
  # select only rows for which DNA is available (since some samples were missing DNA)
  diversity_matrix <-
    diversity_matrix[!(diversity_matrix$DNAfilename == ""), ]

  # add new columns to the matrix for DNA and RNA diversity
  diversity_matrix$DNAdiv <- NA
  diversity_matrix$RNAdiv <- NA
  
  # Loop through the matrix and append the RNA and DNA diversity values
  for (i in rownames(diversity_matrix)) {
    dna_file <- as.character(diversity_matrix[i, ]$DNAfilename)
    rna_file <- as.character(diversity_matrix[i, ]$RNAfilename)
    if (!is.na(dna_file)) {
      new_dna_file <- gsub("-", ".", dna_file)
      diversity_matrix[i, ]$DNAdiv <- simpsons_div[new_dna_file]
    }
    if (!is.na(rna_file)) {
      new_rna_file <- gsub("-", ".", rna_file)
      diversity_matrix[i, ]$RNAdiv <- simpsons_div[new_rna_file]
    }
  }
  
  return(diversity_matrix)
}

#' basic.setup
#'
#' Take matrix of gene counts and metadata and transform them into filtered, ExpressionSet object to be used in downstream analysis.
#' @param genecounts matrix; the matrix of gene counts (samples in columns, genes in rows)
#' @param metadata matrix; the metadata matrix (samples in rows, data in columns)
#' @param genemap matrix; the matrix of ensembl to hgnc genes 
#' @param removePC boolean; boolean variable indicating whether or not to remove non-protein-coding genes (default = TRUE)
#' @param read_filter integer; threshold for number of reads required to include a sample (default = 10,000)
#' @param nonzero_filter float; value indicating the proportion of samples that must have non-zero gene expression value for a gene to be included (detault = .5)
#' @param equalize boolean; a variable indicating whether to subsample the dataset to create equal classes (NO LONGER IN USE)
#' @export
#' @examples
#' basic.setup(genecounts, metadata, removePC = TRUE, read_Filter = 10000, nonzero_filter = .5, equalize=FALSE)
basic.setup <-
  function(genecounts,
           metadata,
           genemap,
           removePC = TRUE,
           read_filter = 10000,
           nonzero_filter = .5,
           equalize = FALSE) {
    
    # if you are subsampling the training data to have equal representation of group 1 and group 2
    if (equalize) {
      metadata_g1 <-
        metadata[metadata$effective_group == 1, ][sample(seq(1, dim(metadata[metadata$effective_group ==
                                                                               1, ])[1]), dim(metadata[metadata$effective_group == 4, ])[1]), ]
      metadata_notg1 <-
        metadata[metadata$effective_group %in% c(2, 3, 4), ]
      metadata <- rbind(metadata_g1, metadata_notg1)
    }

    # make gene counts' sample IDs match the metadata rownames
    # select only data / metadata that are present in both matrices
    genecounts <-
      genecounts[, colnames(genecounts) %in% rownames(metadata)]  
    
    # select only the metadata that has represenation in the genecounts matrix 
    metadata <-
      metadata[rownames(metadata) %in% colnames(genecounts), ]  
    
    # re-order the metadata to match the order of the data
    ord <-
      match(colnames(genecounts), rownames(metadata)) 
    metadata <- metadata[ord,]
    
    if (removePC) {  # if you are only interested in the protein-coding genes, remove all non-pc genes...
      
      #get ENSEMBL gene names 
      genenames <-
        rownames(genecounts)
      
      # create an index in which you are rearranging ENSEMBL gene ids from genemap 
      # into same order as they are in genenames variable
      idx <- match(genenames, genemap$ensembl_gene_id)
      entrez <- genemap$entrezgene[idx]
      hgnc_symbol <- genemap$hgnc_symbol[idx]
      description <- genemap$description[idx]
      gene_biotype <- genemap$gene_biotype[idx]
      ensembl <- genemap$ensembl_gene_id[idx]
      
      #ga = gene annotations
      ga <- as.data.frame(cbind.data.frame(hgnc_symbol,gene_biotype, ensembl, entrez))     #description, 
      
      ga$gene_biotype <- as.character(ga$gene_biotype)
      ga <- ga[complete.cases(ga), ]
      
      #pc = protein coding
      pc <- ga[ga$gene_biotype %in% c("protein_coding"), ]    
      pc <- pc[!(is.na(pc$ensembl)), ]
      
      print("removing MT-*, RPL*, and RPS* genes")
      #index to remote MT-rna from protein coding list
      idxmt <- grep("MT-", pc$hgnc_symbol)    
      pc <- pc[-idxmt, ]
      idxrp <- grep("RPL", pc$hgnc_symbol)
      pc <- pc[-idxrp, ]
      idxrps <- grep("RPS", pc$hgnc_symbol)
      pc <- pc[-idxrps, ]

      # make an index for protein coding genes
      idxpc <- match(pc$ensembl, rownames(genecounts))   
      idxpc <- idxpc[!is.na(idxpc)]
      
      #subset read counts to only protein-coding based on above indexes
      genecounts.proteincoding <- genecounts[idxpc, ]    
      
      # filter samples for minimum number of total reads
      pass_quality_filter <-
        names(colSums(genecounts.proteincoding)[colSums(genecounts.proteincoding) > read_filter])  
      
      genecounts.proteincoding <-
        genecounts.proteincoding[, (colnames(genecounts.proteincoding) %in% pass_quality_filter)]
      
      # of the samples that pass filter, remove genes with many zeros
      genecounts.proteincoding <-
        genecounts.proteincoding[rownames(filterZeroRows(genecounts.proteincoding, nonzero_filter)), ]  
      
      # adjust the metadata to match the filtered data
      metadata <-
        metadata[rownames(metadata) %in% pass_quality_filter, ]  
      
      return(ExpressionSet(
        as.matrix(genecounts.proteincoding),
        phenoData = AnnotatedDataFrame(metadata)
      ))   #return ExpressionSet object
      
    }
    
    # if you didn't remove PC-genes and apply filters, then return the basic aligned metadata/genecounts data
    return(ExpressionSet(as.matrix(genecounts), phenoData = AnnotatedDataFrame(metadata)))   #return ExpressionSet object
    
  }

#' traditionalclassification.setup
#'
#' Take ExpressionSet and apply standard classification filtering
#' @param eset ExpressionSet; the ExpressionSet object (output from basic.setup, binding metadata with genecoutns)
#' @param normalization_method string; which normalization method to use, must be one of "DESeq" or "voom"; default = "DESeq"
#' @param rld boolean; boolean variable indicating whether or not to apply rld transformation; if not, then just return log(normalized counts)
#' @export
#' @examples
#' traditionalclassification.setup(eset, normalization_method = "DESeq", rld = FALSE)
traditionalclassification.setup <-
  function(eset,
           normalization_method = "DESeq",
           rld = FALSE) {
    
    # if you will be using DESeq normalization...
    if (normalization_method == "DESeq") {
      eset$effective_group <- as.factor(eset$effective_group)
      
      # create DESeq2 object
      dds <- DESeqDataSetFromMatrix(exprs(eset), pData(eset), design = ~ effective_group)  
      
      # if using built-in rld transformation...
      if (rld) {  
        rld <- rlog(dds)
        return(ExpressionSet(assay(rld), phenoData = AnnotatedDataFrame(pData(eset))))
      } else{  # else, take the normalized counts by estimateSizeFactors and log-transform
        dds <- estimateSizeFactors(dds)
        normalized <- counts(dds, normalized = TRUE)
        return(ExpressionSet(log(normalized + 1), phenoData = AnnotatedDataFrame(pData(eset))))
      }
      
    } else if (normalization_method == "voom") {  # if not DESeq norm, use VOOM
      v <- voom(eset)
      return(ExpressionSet(v$E, phenoData = AnnotatedDataFrame(pData(eset))))
    }
  }


#' format_host_features
#'
#' Take ExpressionSet and apply z-score scaling and centering of each gene
#' @param eset ExpressionSet; the ExpressionSet object (output from basic.setup, binding metadata with genecoutns)
#' @param zscore boolean; boolean variable indicating whether or not to apply gene-wise z-score
#' @export
#' @examples
#' format_host_features(eset, zscore = FALSE)
format_host_features <- function(eset, zscore = FALSE) {
  eset.orig <- eset
  if (zscore) {
    eset <-
      ExpressionSet(t(scale(
        t(exprs(eset)), center = TRUE, scale = TRUE
      )), phenoData = AnnotatedDataFrame(pData(eset))) #makes each column a z-score
  }
  return(eset)
}

format_host_features_get_z_values <- function(eset, zscore = FALSE) {
  eset.orig <- eset
  s <- c()
  if (zscore) {
    s <- scale(
      t(exprs(eset)), center = TRUE, scale = TRUE
    )
    eset <-
      ExpressionSet(t(s), phenoData = AnnotatedDataFrame(pData(eset))) #makes each column a z-score
  }
  out <- list(eset,attributes(s)$`scaled:center`,attributes(s)$`scaled:scale`)
  names(out) <- c("eset","centering_values","scaling_values")
  return(out)
}


# Function: This function allows you to scale a dataframe (ie metadata), 
#           while leaving known categorical variables unscaled
smart_scale <- function(df) {
  
  # subset original dataframe - these variables should be scaled
  to_scale <-
    df[, !(
      colnames(df) %in% c(
        "effective_group",
        "original_group",
        "is_bacteria",
        colnames(metadata_with_microbe_with_diversity)[sapply(metadata_with_microbe_with_diversity,
                                                                                  is.factor)]
      )
    )]
  
  # subset original dataframe - these variables should not be scaled
  not_to_scale <-
    df[, (
      colnames(df) %in% c(
        "effective_group",
        "original_group",
        "is_bacteria",
        colnames(metadata_with_microbe_with_diversity)[sapply(metadata_with_microbe_with_diversity,
                                                                                  is.factor)]
      )
    )]
  
  scaled <- scale(to_scale, center = FALSE, scale = TRUE)  # scale the subset of original dataframe intended for scaling

  # merge the scaled and un-scaled dataframes
  final <-
    cbind(scaled, df[, colnames(df) %in% c(
      "effective_group",
      "original_group",
      "is_bacteria",
      colnames(metadata_with_microbe_with_diversity)[sapply(metadata_with_microbe_with_diversity,
                                                                                is.factor)]
    )])
  
  # if there were multiple variables unscaled, it will be type = DataFrame,
  # so not_to_scale will have a colname() attribute that can be used to name the merged dataframe
  if (is.data.frame(not_to_scale)) {
    colnames(final) <- c(colnames(to_scale), colnames(not_to_scale))
  } else{
    colnames(final) <-
      c(colnames(to_scale), "effective_group") #if only one group was not scaled, this was the effective_group
  }

  final[is.na(final)] <- 0  # set NA to zero
  return(final)
}
