###
### Script developed by KKalantar
### Date September 2, 2016
### Purpose: Maintain functions that are used across many scripts
### source("C:\\cygwin64\\home\\kkalantar\\ResearchDL\\corefns_kk.R") inside files using these functions
###

## Determine which probes are present for each sample
filterZeroRows <- function(input_matrix, is_nonzero_in_threshold){
  mabsent <- c()
  for (i in 1:ncol(input_matrix)) {
    absent <- as.numeric(input_matrix[,i]) == 0 ## Zero spot =absent
    mabsent <- cbind(mabsent, absent) ## Add the sample to an "absent" matrix with all other samples
  }
  mpresent <- mabsent==0  ## Determine which genes are present
  ## Find which genes are present in at least 30% of the samples
  mfilter<- apply(mpresent,1,sum) > ncol(input_matrix)*is_nonzero_in_threshold
  mpresent30<-input_matrix[mfilter=="TRUE",]
  
  genecounts <- mpresent30 #8/17 
  return( genecounts )
}


filterOutliers <- function(input_matrix, remove_outliers_threshold){
  mremove <- c()
  mpresent30<-data.matrix(input_matrix) #convert to numeric so that the boxplot() method can handle data
  if(remove_outliers_threshold > 0){
    print("removing outliers")
    mpresentoutlier <- c()
    rownames_mpresentoutlier <- c()
    for (i in 1:nrow(mpresent30)) {
      # number of outliers less than threshold, but greater than 0 = 
      bpo<-length(boxplot.stats(mpresent30[i,])$out)
      print(bpo)
      remove <- as.numeric(0 < bpo && bpo < ncol(mpresent30)*remove_outliers_threshold)
      # Add the sample to a "remove" matrix with all other samples
      mremove <- cbind(mremove, remove)
      if(remove==0){
        mpresentoutlier<- rbind(mpresentoutlier, mpresent30[i,])
        rownames_mpresentoutlier <- cbind(rownames_mpresentoutlier, row.names(mpresent30)[i])
      }
    }
    
    row.names(mpresentoutlier)<-rownames_mpresentoutlier
    #dim(mpresentoutlier) 
    return(mpresentoutlier)
  }
}

filter_mitochondrial <- function(input_matrix){
  idxmt <- grep("MT-",rownames(input_matrix))
  output_matrix<-input_matrix[-idxmt,]
  return(output_matrix)
}

#function assumes that data is in the form of ENSEMBL IDs in the rows
#filters for only genes based on biotype filter from ENTREZ annotation
filter_biotype <- function(read_counts, filter){
  library( "biomaRt" )
  genenames<- sapply( strsplit( rownames(read_counts), split="\\." ), "[", 1 )
  listMarts(host="www.ensembl.org")
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "gene_biotype", "description"),
                    filters = "ensembl_gene_id",
                    values = genenames,
                    mart)
  
  #create an index in which you are rearranging ensmbl gene ids from genemap into same order as they are in genenames
  idx <- match(genenames, genemap$ensembl_gene_id )
  
  entrez <- genemap$entrezgene[ idx ]
  hgnc_symbol <- genemap$hgnc_symbol[ idx ]
  description <- genemap$description[ idx ]
  gene_biotype <- genemap$gene_biotype[ idx ]
  ensembl <- genemap$ensembl_gene_id[ idx ]
  
  #ga = gene annotations
  ga<-cbind(hgnc_symbol,description,gene_biotype,ensembl,entrez)
  #make ga into a data frame from a matrix
  ga<-as.data.frame(ga)
  ga$gene_biotype<-as.character(ga$gene_biotype)
  
  pc<-ga[ga$gene_biotype%in%filter,]
  length(pc)
  pc<-pc[!(is.na(pc$ensembl)),]
  
  if ("MT" %in% filter){
    pc<-pc[!("MT-"%in%pc$hgnc_symbol),]
  }
  
  
  return(read_counts[pc$ensembl,])
}


runLASSO <- function(input_matrix, results, metadata){
  genes_with_05padj<-input_matrix[row.names(input_matrix)%in%row.names(results),]
  x <- t(as.matrix(genes_with_05padj))
  y <- as.numeric(metadata$infxn)
  fit <- glmnet(x, y, family="gaussian", alpha=0.5, lambda=0.001)  # fit model
  summary(fit) # summarize the fit
  c <- as.matrix(coef(fit))
  inds<-which(c!=0)
  variables<-row.names(c)[inds]
  variables<-variables[!variables %in% '(Intercept)']
  data_in_features = input_matrix[row.names(input_matrix)%in%variables,]
  dim(data_in_features)
  output = data_in_features
  print(length(rownames(output)))
  return(output)
}

runLASSOonly <- function(input_matrix, metadata){
  genes_with_05padj<-input_matrix
  x <- t(as.matrix(genes_with_05padj))
  y <- as.numeric(metadata$infxn)
  fit <- glmnet(x, y, family="gaussian", alpha=0.5, lambda=0.001)  # fit model
  summary(fit) # summarize the fit
  c <- as.matrix(coef(fit))
  inds<-which(c!=0)
  variables<-row.names(c)[inds]
  variables<-variables[!variables %in% '(Intercept)']
  data_in_features = input_matrix[row.names(input_matrix)%in%variables,]
  dim(data_in_features)
  output = data_in_features
  print(length(rownames(output)))
  return(output)
}