

#
# SECTION #1: Utility functions for classification analyses
#


#' Leave One Out Split
#'
#' This function splits input data for leave-one-out cross-validation
#' @param df full data matrix
#' @keywords split
#' @export
#' @examples
#' tenfold.split(myData)

leave_one_out.split <- function(df) {
  training_data <- list()
  test_data <- list()
  for (j in seq(1:dim(df)[1])) {
    samp <- j
    df.train <- rownames(df[-samp,])
    df.test <- rownames(df[samp,])
    training_data[[j]] <- df.train
    test_data[[j]] <- df.test
  }
  return(list(training_data, test_data))
}



#' #' Leave One Out Equalize Split
#' #'
#' #' This function splits input data for train/test set with equalized training and test sets with some samples left out
#' #' @param eset the ExpressionSet object containing all of the individuals in group 1 and group 4
#' #' @keywords split
#' #' @export
#' #' @examples
#' #' leave_one_out.equalize.split(ExpressionSet)
#' #'
#' leave_one_out.equalize.split <- function(df) {
#'   g1 <- classification_eset[, classification_eset$effective_group == 1]
#'   g4 <-
#'     classification_eset[, classification_eset$effective_group == 4]
#'   
#'   n <- floor(.9 * (min(length(colnames(
#'     g1
#'   )), length(colnames(
#'     g4
#'   )))))
#'   
#'   training_data <-
#'     c(sample(colnames(g1), n, replace = FALSE),
#'       sample(colnames(g4), n, replace = FALSE))
#'   
#'   eligible <- c()
#'   used <- c()
#'   test_data <- c()
#'   
#'   if (length(colnames(g4)) < length(colnames(g1))) {
#'     eligible <- colnames(g1)[!(colnames(g1) %in% training_data)]
#'     used <-
#'       sample(eligible, length(colnames(g4)) - n, replace = FALSE)
#'     test_data <-
#'       c(used, colnames(g4)[!(colnames(g4) %in% training_data)])
#'   } else {
#'     eligible <- colnames(g4)[!(colnames(g4) %in% training_data)]
#'     used <-
#'       sample(eligible, length(colnames(g1)) - n, replace = FALSE)
#'     test_data <-
#'       c(colnames(g1)[!(colnames(g1) %in% training_data)], used)
#'   }
#'   
#'   #test_data <- c(colnames(g1)[!(colnames(g1) %in% training_data)],colnames(g4)[!(colnames(g4) %in% training_data)])
#'   
#'   return(list(training_data, test_data))
#' }






#' #' Specific Split - TODO: WARNING MAY NOT WORK AS EXPECTED, FIX
#' #'
#' #' Splits the data in a particular way - currently specific to current version with 32 effective_group 1s and 13 effective_group 4s
#' #' @param df full data matrix
#' #' @param folds integer number of folds to return
#' #' @keywords split
#' #' @export
#' #' @examples
#' #' tenfold.split(myData, 5)
#' 
#' tenfold.split <- function(df, folds = 10) {
#'   #folds = 10
#'   folds.g1 <-
#'     split(sample(which(df$effective_group == 1)),
#'           c(
#'             rep(1, 3),
#'             rep(2, 3),
#'             rep(3, 3),
#'             rep(4, 3),
#'             rep(5, 3),
#'             rep(6, 3),
#'             rep(7, 3),
#'             rep(8, 3),
#'             rep(9, 2),
#'             rep(10, 2)
#'           ))
#'   folds.g4 <-
#'     split(sample(which(df$effective_group == 4)), 
#'       c(
#'         rep(1, 2),
#'         rep(2, 2),
#'         rep(3, 2),
#'         rep(4, 2),
#'         rep(5, 2),
#'         rep(6, 1),
#'         rep(7, 1),
#'         rep(8, 1),
#'         rep(9, 1),
#'         rep(10, 1)
#'     ))
#'   
#'   training_data <- list()
#'   test_data <- list()
#'   
#'   for (j in seq(1:folds)) {
#'     samp <- c(folds.g1[[j]], folds.g4[[j]])
#'     #print(samp)
#'     df.train <- colnames(df)[-samp]
#'     df.test <- colnames(df)[samp]
#'     training_data[[j]] <- df.train
#'     test_data[[j]] <- df.test
#'   }
#'   return(list(training_data, test_data))
#' }






#' #' classify_by_subset.loocv
#' #'
#' #' Take matrix of gene counts and metadata and transform them into filtered, ExpressionSet object to be used in downstream analysis.
#' #' @param train ExpressionSet of training data
#' #' @param test ExpressionSet of test data
#' #' @param variables_of_interest list of metadata variables to use in prediction
#' #' @param plot Boolean, indicating whether to plot data
#' #' @param method "glm" or "rf" or "glmnet"
#' #' @param optimize "acc" or "npv"
#' #' @keywords classify clinical
#' #' @export
#' #' @examples
#' #' classify_using_clinical_data.loocv(train, test, plot=FALSE, method="glm", optimize="acc"))
#' 
#' 
#' classify_by_subset.loocv <-
#'   function(train,
#'            test,
#'            variables_of_interest,
#'            plot = FALSE,
#'            method = "glm",
#'            optimize = "acc") {
#'     current_data <-
#'       pData(train)[, c("effective_group", variables_of_interest)]
#'     current_data[, c("effective_group")] <-
#'       as.numeric(current_data[, c("effective_group")] == 1)
#'     
#'     if (method == "glm") {
#'       model <-
#'         glm(
#'           formula = as.formula(paste(
#'             "effective_group ~ ",
#'             paste(variables_of_interest, collapse = "+")
#'           )),
#'           family = binomial(link = 'logit'),
#'           data = as.data.frame(train)
#'         )
#'       threshold <- .5
#'       
#'       # if (plot) {
#'       #   barplot(
#'       #     predict(model, test, type = "response") - .5,
#'       #     col = barplot_colors[test$effective_group + 1],
#'       #     main = 'clinical-only predictions',
#'       #     las = 2,
#'       #     cex.main = .8
#'       #   )
#'       #   abline(h = threshold)
#'       # }
#'       
#'       output <- list(
#'         predict(model, test, type = "response"),
#'         predict(model, train, type = "response"),
#'         threshold
#'       )
#'       names(output) <-
#'         c("test_classification",
#'           "train_classification",
#'           "threshold")
#'       return(output)
#'       
#'       # data_matrix = as.matrix(pData(train)[,variables_of_interest])
#'       # data_matrix[is.na(data_matrix)] <- 0
#'       # print(length(as.integer(train$effective_group==1)))
#'       # model <- glmnet(data.matrix(data_matrix), as.integer(train$effective_group==1), family=c("binomial"),alpha=0,lambda = .6)
#'       # threshold <- .5
#'       #
#'       # coef <- as.matrix(coef(model))      #get coefficients for all variables at that lambda
#'       # coef_used <-
#'       #   coef[coef != 0, ]       #this gives the genes which were used in prediction with that lambda value that are not set to zero - many IFN genes
#'       # features_of_interest <- names(coef_used)
#'       #
#'       # output <- list(
#'       #   predict(model, as(as.matrix(pData(test)[,variables_of_interest]),"dgCMatrix"))[, c("s0")],
#'       #   features_of_interest,
#'       #   predict(model, as(as.matrix(pData(train)[,variables_of_interest]),"dgCMatrix"))[, c("s0")],
#'       #   threshold
#'       # )
#'       #
#'       # names(output) <- c("test_classification","features","train_classification","threshold")
#'       # return(output)
#'       
#'     } else if (method == "glmnet") {
#'       data_matrix = as.matrix(pData(train)[, variables_of_interest])
#'       data_matrix[is.na(data_matrix)] <- 0
#'       model <-
#'         glmnet(
#'           data.matrix(data_matrix),
#'           as.integer(train$effective_group == 1),
#'           family = c("binomial"),
#'           alpha = .2,
#'           lambda = .6
#'         )
#'       threshold <- .5
#'       
#'       coef <-
#'         as.matrix(coef(model))      #get coefficients for all variables at that lambda
#'       coef_used <-
#'         coef[coef != 0,]       #this gives the genes which were used in prediction with that lambda value that are not set to zero - many IFN genes
#'       features_of_interest <- names(coef_used)
#'       
#'       output <- list(predict(model, as(as.matrix(pData(
#'         test
#'       )[, variables_of_interest]), "dgCMatrix"))[, c("s0")],
#'       features_of_interest,
#'       predict(model, as(as.matrix(pData(
#'         train
#'       )[, variables_of_interest]), "dgCMatrix"))[, c("s0")],
#'       threshold)
#'       
#'       names(output) <-
#'         c("test_classification",
#'           "features",
#'           "train_classification",
#'           "threshold")
#'       return(output)
#'       
#'     } else if (method == "rf") {
#'       train$effective_group <- factor(train$effective_group)
#'       test$effective_group <- factor(test$effective_group)
#'       levels(test) <- levels(train)
#'       train <- as.data.frame(train)
#'       train[is.na(train)] <- 0
#'       model <-
#'         randomForest(#as.factor(effective_group) ~  TempMinSAPS + TempMaxSAPS + WBCMinSAPS + WBCMaxSAPS + APACHEIII,
#'           as.formula(paste(
#'             "effective_group ~ ",
#'             paste(variables_of_interest, collapse = "+")
#'           )),
#'           data = as.data.frame(train),
#'           classwt = class_weight)
#'       pred.obj <-
#'         prediction(predict(model, train, type = "prob")[, 2],
#'                    train$effective_group)
#'       
#'       if (optimize == "acc") {
#'         #get the best threshold
#'         ss <- performance(pred.obj, "sens", "spec")
#'         threshold <-
#'           ss@alpha.values[[1]][which.max(ss@x.values[[1]] + ss@y.values[[1]])]
#'       } else if (optimize == "npv") {
#'         #get the threshold which maximizes negative predictive value
#'         npvss <- performance(pred.obj, "npv")
#'         threshold <-
#'           npvss@x.values[[1]][which.max(npvss@y.values[[1]])]
#'       }
#'       
#'       # if (plot) {
#'       #   barplot(
#'       #     predict(model, test, type = "prob")[, c("1")],
#'       #     col = barplot_colors[test$effective_group + 1],
#'       #     main = 'microbe-only predictions',
#'       #     las = 2,
#'       #     cex.main = .8
#'       #   )
#'       #   abline(h = threshold)
#'       # }
#'       
#'       output <- list(
#'         predict(model, test, type = "prob")[, c("1")],
#'         predict(model, train, type = "prob")[, c("1")],
#'         threshold
#'       )
#'       names(output) <-
#'         c("test_classification",
#'           "train_classification",
#'           "threshold")
#'       return(output)
#'       
#'     }
#'     
#'   }
