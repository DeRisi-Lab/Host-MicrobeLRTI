#' Raw gene counts
#'
#' This data contains RNA-seq gene counts for all patients sequencing in mBAL study.
#'
#' @format A data frame with [number of genes] rows and [number of patients] columns:
#' \describe{
#'   \item{genes}{All genes}
#'   \item{patients}{All patients}
#'   ...
#' }
#'
"raw_genecounts"

#' Raw metadata
#'
#' This data contains clinial metadata for all patients sequencing in mBAL study.
#'
#' @format A data frame with [number of patients] rows and [number of clinical features] columns:
#' \describe{
#'   \item{clinical features}{clinical and NGS-based metadata features collected for the study}
#'   \item{patients}{All patients}
#'   ...
#' }
#'
"raw_metadata"

#' mBAL filtered ExpressionSet
#'
#' This data contains the quality-filtered genecounts aligned to the mBAL metadata patient IDs;
#' The mBAL metadata contains diversity, microbe, and microbe_vector features at this point.
#'
#' @format An ExpressionSet object with filtered gene counts data in the exprs() field and metadata in the pData() field
#'
"filtered_eset"

#' hostCD
#'
#' This data contains the quality-filtered genecounts aligned to the mBAL metadata patient IDs;
#' The mBAL metadata contains diversity, microbe, and microbe_vector features at this point. NO scaling has been done to metadata features.
#'
#' @format An ExpressionSet object with filtered AND NORMALIZED gene counts data in the exprs() field and metadata in the pData() field
#'
"fulldata_eset"


#' classification_eset
#'
#' This data contains the quality-filtered genecounts aligned to the mBAL metadata patient IDs for only group 1 and group 4 patients;
#' This includes all features for training on. SCALING has been applied to all metadata features
#'
#' @format An ExpressionSet object with filtered AND NORMALIZED gene counts data in the exprs() field and metadata in the pData() field
#'
"classification_eset"

#' unknown_eset
#'
#' This data contains the quality-filtered genecounts aligned to the mBAL metadata patient IDs for only group 2 and group 3 patients;
#' This includes all features for training on. SCALING has been applied to all metadata features
#'
#' @format An ExpressionSet object with filtered AND NORMALIZED gene counts data for group 2 and group 3 in the exprs() field and metadata in the pData() field
#'
"unknown_eset"

#' genemap
#'
#' This matrix contains a map between ensembl and HGNC gene identifiers
#'
#' @format A matrix containing ensembl gene ids and hgnc gene ids
#'
"genemap"


#' tsalik_ensembl
#'
#' A list of ensembl identifiers for genes used in the Tsalik et al paper
#'
#' @format A list of ensembl identifiers for genes used in the Tsalik et al paper
#'
"tsalik_ensembl"


#' tsalik_ensembl
#'
#' A list of ensembl identifiers for genes used in the Suarez et al paper
#'
#' @format A list of ensembl identifiers for genes used in the Suarez et al paper
#'
"suarez_ensembl"


#' cell_proportions
#'
#' A matrix of estimated cell type proportions estimated from EpiDISH
#'
#' @format A matrix with rows corresponding to sample IDs and columns corresponding to cell types
#'
"cell_proportions"
