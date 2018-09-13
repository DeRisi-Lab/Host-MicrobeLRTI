

#' classify_using_host_transcriptome.loocv
#'
#' Use this method to run host transcriptome classification. 
#' It takes in a training and test dataset consisting of gene expression values for each individual and 
#' returns the performance of the host expression classifier on the test samples. 
#' The basis of this classificaiton function is regularized linear regression, although use of a random forest can be specified.
#' @param eset ExpressionSet output from traditionalclassification.setup
#' @param train list of training samples
#' @param test List of test samples
#' @param control list of covariats to control
#' @param master_genemap matrix of HGNC to ENSEMBL ID mappings
#' @param full_model Boolean, whether to use standard logistic regression probabilities
#' @param save_pdf Boolean, whether to generate a saved PDF
#' @param plot Boolean, whether to plot the training features in a heatmap
#' @param alpha elastic net alpha parameter
#' @return List of outputs including 
#' 1. test_classification - numeric vector of classification probabilities for each test sample, 
#' 2. subset - matrix of samples x genes selected as predictive, 
#' 3. train_classification - numeric vector of classification probabilities for each training sample, 
#' @keywords classify host
#' @export
#' @examples
#' classify_using_host_transcriptome.loocv(eset, train, test, master_genemap, full_model=FALSE, save_PDF=FALSE, plot=FALSE)

classify_using_host_transcriptome.loocv <-
  function(eset,
           train,
           test,
           control,
           master_genemap,
           full_model = FALSE,
           save_PDF = FALSE,
           plot = FALSE,
           alpha = .5) {
    

    eset.train <-
        eset[, colnames(eset) %in% train]  # create training data based on input names
    
    scaled_training <- scale(t(exprs(eset.train)))
    center_vals <- attributes(scaled_training)[3]
    scale_vals <- attributes(scaled_training)[4]
    exprs(eset.train) <- t(as.matrix(scaled_training))
    
    eset.test <-
        eset[, colnames(eset) %in% test]   # create test data based on input names
    
    centered_test <- exprs(eset.test) - center_vals[[1]]
    scaled_test<-centered_test/scale_vals[[1]]
    exprs(eset.test) <- as.matrix(scaled_test)  # try to put the scaled data back into the original exprs object...NOT SURE THIS WILL WORK
    
    #print(dim(eset.train))
    #print(dim(eset.test))
    
    
    g_fit_cor <- c()  # declare variables
    m <- c()
    m_test <- c()
    
    # if there are covariates to control for, create model to correct for batch effect
    if (length(control) > 0) {
      m <-
        rbind(exprs(eset.train),
              t(sapply(pData(eset.train)[, control], as.numeric)),
              as.integer(pData(eset.train)[, 'effective_group'] == 1))   # full matrix of genes and features to control for (training)
      
      #print(dim(m))
      
      rownames(m) <-
        c(rownames(exprs(eset.train)),
          control,
          'effective_group')  # set rownames for matrix
      
      #print(tail(rownames(m)))
      
      eset_m <- ExpressionSet(as.matrix(m),phenoData = AnnotatedDataFrame(pData(eset.train)))  # create modified ExpressionSet eset_m for training
      
      if (ncol(exprs(eset.test)) > 1) {
        m_test <-
          as.matrix(rbind(exprs(eset.test),
                          t(sapply(
                            pData(eset.test)[, control], as.numeric
                          )),
                          as.integer(pData(eset.test)[, 'effective_group'] == 1))) #as.matrix()    # full matrix of genes and features to control for (test)
      
        rownames(m_test) <-
          c(rownames(exprs(eset.train)),
            control,
            'effective_group')  
      } else{
        
        #print("hello")
        m_test <-
          as.matrix(c(exprs(eset.test),
                      t(sapply(
                        pData(eset.test)[, control], as.numeric
                      )),
                      as.integer(pData(eset.test)[, 'effective_group'] == 1))) #as.matrix()    # full matrix of genes and features to control for (test)
        m_test <- as.matrix(m_test)
        rownames(m_test) <-
          c(rownames(exprs(eset.train)),
            control,
            'effective_group')    
        colnames(m_test) <- colnames(eset.test)
      }
      
      eset_mtest <- ExpressionSet(as.matrix(m_test),phenoData = AnnotatedDataFrame(pData(eset.test)))  # create modified ExpressionSet eset_m for training
      
      rownames(m_test) <-
        c(rownames(exprs(eset.test)),
          control,
          'effective_group')  # set rownames for matrix
      colnames(m_test) <-
        colnames(test)  # set colnames for matrix
    }
    #print(dim(m))
    
    
      #print("control GREATER than 1")
      eset.withcovs <-
        ExpressionSet(m, phenoData = AnnotatedDataFrame(pData(eset.train)))
      
      #print(dim(exprs(eset.withcovs)))
      
      
      # cv <-
      #   cv.glmnet(    # used in case we want to determine lambda by cross-validation
      #     t(m[!(rownames(m) %in% c("effective_group")), ]),
      #     as.integer(eset.train$effective_group == 1),
      #     family = c("binomial"),
      #     alpha = alpha
      #   )
      # print("LAMBDA")
      # print(cv$lambda.1se)
      
      
      #APPLY A CORRECTION
      #fit <- glm(eset.train$batch)
      #x <- as.data.frame(cbind(classification_eset$effective_group, as.factor(classification_eset$batch), t(exprs(classification_eset))))
      #colnames(x) <- c("effective_group","batch",colnames(t(exprs(classification_eset))))
      #g <- glm(effective_group ~ batch, family="gaussian", data = x)
      
      
      
      #print(dim(t(m[!(rownames(m) %in% c("effective_group")),])))
      
      g_fit_cor <-
        glmnet(
          t(m[!(rownames(m) %in% c("effective_group")),]),
          as.integer(eset.train$effective_group == 1),
          family = c("binomial"),
          alpha = alpha,
          lambda = .6 #cv$lambda.1se #.6
        )
    
    #print(g_fit_cor)
    b <- as.matrix(g_fit_cor$beta)
    b <- b[b != 0,]
    
    #g_fit_cor <- glmnet(t(exprs(eset.train)), model$residuals, family=c("gaussian"),alpha=.2,lambda = .6)
    cmat_cor <-
      as.matrix(coef(g_fit_cor))  # get coefficients for all variables at that lambda
    dmat_cor <-
      cmat_cor[cmat_cor != 0,]  # this gives the genes which were used in prediction with that lambda value that are not set to zero - many IFN genes
    
    #print(head(dmat_cor))
    genes_of_interest <-names(dmat_cor)[2:length(dmat_cor)] # names(b)
    
    #ADDED 1/5/18
    #genes_of_interest <- c(genes_of_interest, c("ENSG00000119922","ENSG00000169245","ENSG00000135114"))
    
    #print(length(genes_of_interest))
    
    #REMOVE UNWANTED GENES OF INTEREST - ie if a covariate is in the predictive genes - NOTE: don't need to use this anymore, fixed elsewhere
    #genes_of_interest <- intersect(genes_of_interest,rownames(exprs(eset.train)))
            
    number_of_features_cor <- length(genes_of_interest)
    eset.semi.mvp <- eset_m[genes_of_interest,]
    eset.semi.Y <-
      eset.semi.mvp[, eset.semi.mvp$effective_group %in% c(1)]
    eset.semi.N <-
      eset.semi.mvp[, eset.semi.mvp$effective_group %in% c(4)]
    weights <-
      rowMeans(exprs(eset.semi.Y)) - rowMeans(exprs(eset.semi.N))
    weights[weights > 0] <- 1
    weights[weights < 0] <- -1
    
    
    df <- exprs(eset.semi.mvp) * weights
    df_ordered <- exprs(eset.semi.mvp)[, order(colSums(df))]
    metadata_ordered <- pData(eset.semi.mvp)[colnames(df_ordered),]
    eset_ordered2 <-
      ExpressionSet(as.matrix(df_ordered), phenoData = AnnotatedDataFrame(metadata_ordered))
    eset.semi.mvp.test <- eset_mtest[genes_of_interest,] #eset.test[genes_of_interest,]
    df3 <- exprs(eset.semi.mvp.test)
    colnames(df3) <- colnames(exprs(eset.semi.mvp.test))
    rownames(df3) <- rownames(exprs(eset.semi.mvp.test))
    df3 <- df3 * weights
    
    df_ordered <- exprs(eset.semi.mvp.test)[, order(colSums(df3))]
    
    final_projection <-
      cbind(exprs(eset_ordered2), df_ordered)#exprs(eset_ordered3))
    colnames(final_projection) <-
      c(colnames(exprs(eset_ordered2)), colnames(exprs(eset.test)))
    fpsums <- cbind(df, df3)
    fpsums <- fpsums[, order(colSums(fpsums))]
    fp <- final_projection[, colnames(fpsums)]#final_projection))]
    metadata_ordered <- pData(eset)[colnames(fp),]
    n <-
      ExpressionSet(as.matrix(fp), phenoData = AnnotatedDataFrame(metadata_ordered))
    
    full <-
      ExpressionSet(cbind(exprs(eset.semi.mvp), exprs(eset.semi.mvp.test)), phenoData = AnnotatedDataFrame(rbind(
        pData(eset.semi.mvp), pData(eset.semi.mvp.test)
      )))
    
    if (save_PDF) {
      mat1 <- exprs(full)
      # b <-
      #   master_genemap[master_genemap$ensembl_gene_id %in% rownames(mat1),]
      # rownames(mat1) <-
      #   b[order(match(b$ensembl_gene_id, rownames(mat1))),]$hgnc_symbol
      mat2 <- exprs(n)
      # b <-
      #   master_genemap[master_genemap$ensembl_gene_id %in% rownames(mat2),]
      # rownames(mat2) <-
      #   b[order(match(b$ensembl_gene_id, rownames(mat2))),]$hgnc_symbol
      
      heatmap.2(
        mat1,
        ColSideColors = infxn_colors[full$effective_group],
        trace = "none",
        key = TRUE,
        col = myPalette,
        margin = c(10, 10),
        dendrogram = "both",
        main = test
      )
      
      heatmap.2(
        mat2,
        ColSideColors = infxn_colors[n$effective_group],
        labCol = colnames(exprs(n)),
        trace = "none",
        key = TRUE,
        col = myPalette,
        margin = c(10, 10),
        Colv = FALSE,
        dendrogram = "row",
        main = test
      )
    }
    
    if (length(control) < 1) {
      mat2use.train <- as.matrix(t(exprs(eset.train)))
      mat2use.test <- as.matrix(t(exprs(eset.test)))
    }
    else{
      mat2use.train <- t(m[!(rownames(m) %in% c("effective_group")), ])
      mat2use.test <-
        t(m_test[!(rownames(m_test) %in% c("effective_group")), ])
      rownames(mat2use.test) <- colnames(test)
    }
    
    if (full_model) {
      pred.obj <-
        prediction(
          predict(g_fit_cor, mat2use.train, type = "response")[, c("s0")],
          as.numeric(eset.train$effective_group == 1)  # eset.train
        )
      
      if (plot) {
        barplot(
          predict(g_fit_cor, mat2use.test, type = "response")[, c("s0")],
          col = barplot_colors[(eset.test$effective_group)],
          main = 'host-only full predictions',
          las = 2,
          cex.main = .8
        )
      }
      
      test_output <-
        as.double(predict(g_fit_cor, mat2use.test)[, c("s0")])
      names(test_output) <- rownames(mat2use.test)
      output <- list(
        test_output,
        as.data.frame(t(exprs(eset.semi.mvp))),
        predict(g_fit_cor, mat2use.train, type = "response")[, c("s0")],
        b,
        g_fit_cor,
        weights
      )
      names(output) <-
        c("test_classification",
          "subset",
          "train_classification",
          "betas",
          "model",
          "weights")
      return(output)
    } else{
      if (plot) {
        barplot(
          colSums(fpsums)[test],
          col = barplot_colors[eset.test$effective_group + 1],
          main = 'host-only score predictions',
          las = 2,
          cex.main = .8
        )
      }
      
      output <- list(colSums(fpsums)[test],
                     exprs(full),#as.data.frame(t(exprs(eset.semi.mvp))),
                     colSums(fpsums)[train],
                     b,
                     g_fit_cor,
                     weights,
                     scale_vals,
                     center_vals)
      names(output) <-
        c("test_classification",
          "subset",
          "train_classification",
          "betas",
          "model",
          "weights",
          "scale_vals",
          "center_vals")
      return(output)
      #return(colSums(fpsums)[test])
    }
    
    
  }



#' classify_using_tsalik_et_al.loocv
#'
#' Take pre-specified gense of interest, along with a training and test set consisting of host transcriptome gene counts and 
#' return the performance of the host expression classifier trained on that subset of genes. This function was developed to evaluate the performance of existing 
#' host gene expression classifiers from the literature.
#' @param genes_of_interest list of genes from the original paper
#' @param eset ExpressionSet output from traditionalclassification.setup
#' @param train list of training samples
#' @param test List of test samples
#' @param full_model Boolean, whether to use standard logistic regression probabilities
#' @param save_pdf Boolean, whether to generate a saved PDF
#' @param plot Boolean, whether to plot the training features in a heatmap
#' @keywords classify host tsalik
#' @return List of outputs including 
#' 1. test_classification - numeric vector of classification probabilities for each test sample, 
#' 2. subset - matrix of samples x genes selected as predictive, 
#' 3. train_classification - numeric vector of classification probabilities for each training sample, 
#' @export
#' @examples
#' classify_using_tsalik_et_al.loocv(genes_of_interest, eset, train, test, full_model = FALSE, save_PDF=FALSE, plot=FALSE)

classify_using_tsalik_et_al.loocv <-
  function(genes_of_interest,
           eset,
           train,
           test,
           full_model = FALSE,
           save_PDF = FALSE,
           plot = FALSE) {
    #genes taken from Supplemental table S5  http://stm.sciencemag.org/content/scitransmed/suppl/2016/01/15/8.322.322ra11.DC1/8-322ra11_SM.pdf
    eset.train <- eset[, colnames(eset) %in% train]
    eset.test <- eset[, colnames(eset) %in% test]
    
    eset.semi.mvp <-
      eset.train[(rownames(eset.train) %in% genes_of_interest),]
    eset.semi.mvp.test <-
      eset.test[rownames(eset.test) %in% genes_of_interest,]
    
    new_mat <-
      rbind(as.integer(eset.train$effective_group == 1), scale(exprs(eset.semi.mvp)))
    rownames(new_mat) <-
      c("effective_group", rownames(exprs(eset.semi.mvp)))
    g_fit_cor <-
      glm(
        formula = effective_group ~ .,
        family = binomial(link = 'logit'),
        data = as.data.frame(t(new_mat))
      )
    
    if (full_model) {
      if (plot) {
        barplot(
          predict(g_fit_cor, as.data.frame(t(
            exprs(eset.semi.mvp.test)
          )), type = "response"),
          col = barplot_colors[(eset.test$effective_group == 1) + 1],
          main = 'host-only full predictions',
          las = 2,
          cex.main = .8
        )
      }
      
      output <- list(
        predict(g_fit_cor, as.data.frame(t(
          exprs(eset.semi.mvp.test)
        )), type = "response"),
        as.data.frame(t(exprs(eset.semi.mvp))),
        predict(g_fit_cor, as.data.frame(t(
          exprs(eset.semi.mvp)
        )), type = "response")
      )
      names(output) <-
        c("test_classification",
          "subset",
          "train_classification")
      return(output)
      
      
      
    } else{
      #range01 <- function(x){(x-min(x))/(max(x)-min(x))}
      if (plot) {
        barplot(
          colSums(fpsums)[test],
          col = barplot_colors[eset.test$effective_group + 1],
          main = 'host-only score predictions',
          las = 2,
          cex.main = .8
        )
      }
      return(colSums(fpsums)[test])
    }
    
  }













#' classify_using_host_transcriptome.withmodel
#'
#' Use this method to run host transcriptome classification. 
#' It takes in a training and test dataset consisting of gene expression values for each individual and 
#' returns the performance of the host expression classifier on the test samples. 
#' The basis of this classificaiton function is regularized linear regression, although use of a random forest can be specified.
#' @param eset.test ExpressionSet output from traditionalclassification.setup
#' @param weights List variable containing the host weights
#' @param center_params List of genenames with parameters for centering data (obtained by centering training data); default is the centering parameters from the mBAL Study training data
#' @param scale_params List of genenames with parameters for scaling data (obtained by scaling training data); default is the scaling parameters from the mBAL Study training data
#' @return List of outputs including 
#' 1. test_classification - numeric vector of classification probabilities for each test sample, 
#' @keywords classify host
#' @export
#' @examples
#' classify_using_host_transcriptome.withmodel(test_eset, weights, center_params, scale_params)

classify_using_host_transcriptome.withmodel <-
  function(eset.test,
           weights,
           center_params = scaling_values_from_training_data$centering_values,
           scale_params = scaling_values_from_training_data$scaling_values) {
    
    common_genes <- intersect(names(scaling_values_from_training_data$centering_values),rownames(exprs(eset.test)))
    eset.test <- eset.test[common_genes,]
    center_params <- center_params[common_genes]
    scale_params <- scale_params[common_genes]
    
    
    centered_test <- exprs(eset.test) - center_params
    scaled_test<-centered_test/scale_params
    exprs(eset.test) <- as.matrix(scaled_test)  # try to put the scaled data back into the original exprs object...NOT SURE THIS WILL WORK
    
    genes_of_interest <-names(weights)
    
    out <- colSums(exprs(eset.test[names(weights),]) * weights)
    
    output <- list(out)
    names(output) <- c("test_classification")
    return(output)
  }



