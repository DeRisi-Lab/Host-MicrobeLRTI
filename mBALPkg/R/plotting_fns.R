#
# FUNCTIONS THAT PLOT ANY OUTPUT USED MULTPLE TIMES
#

#' cohort summary
#'
#' Plot a variety of summary metrics about the cohort including...
#' effective group ID (Grouped by clinical diagnosis into Groups 1 - 4), 
#' pathogen type, 
#' gram stain,
#' gender,
#' age
#' Note: the metadata plots are hard-coded to work with the metadata of this particular study.
#' @param eset ExpressionSet object containing the genecounts and metadata for the cohort
#' @keywords cohort
#' @export
#' @examples
#' cohort_summary(eset)

cohort_summary <- function(eset) {
  
  # generate barpltos of basic summary stats
  barplot(colSums(exprs(eset)), col = infxn_colors[eset$effective_group], las =
            2)  # barplot of total aligned reads
  barplot(colSums(exprs(eset) > 0), col = infxn_colors[eset$effective_group], las =
            2)  # barplot of total number of unique genes covered
  table(eset$effective_group)
  
  # metadata plots:
  
  # effective group ID 
  # Group 1 = CDC Pneumonia + Microbiology confirmation; 
  # Group 2 = CDC Pneumo, no pathogen detected; 
  # Group 3 = unknown; 
  # Group 4 = Non-infectious LRT Disease
  barplot(
    as.matrix(table(eset$effective_group)),
    col = c("red3", "orangered", "royalblue1", "blue3"),
    legend = TRUE
  )
  print(table(eset$effective_group))
  
  # pathogen type - Virus, Bacteria, Fungus, Co-infxn
  barplot(
    as.matrix(table(eset[, eset$effective_group == 1]$pathogenType)),
    col = c(
      "grey",
      "green3",
      "turquoise2",
      "orange2",
      "turquoise4",
      "grey"
    ),
    legend = TRUE
  )
  print(table(eset[, eset$effective_group == 1]$pathogenType))
  
  # gender - M/F
  barplot(
    as.matrix(table(eset$Gender)),
    col = c("grey", "slateblue2", "tomato2", "tomato2"),
    legend = TRUE
  )
  
  # age - numeric value
  hist(
    eset$age,
    breaks = 20,
    xlab = "Age",
    main = "Distribution of Participant Age",
    col = 'slategray4'
  )
}


#' #' plot auc
#' #'
#' #' Plot the ROC curves for input model; will plot one curve per round of cross-validation if multiple rounds were done.
#' #' Note: plots with not be smoothed unless smooth=TRUE parameter is passed when generating the pROC::ROC ojects initially.
#' #' @param predictionObj list containing ROCR auc curve objects for each run of cross-validation
#' #' @param description Title for plot
#' #' @param c color for lines
#' #' @keywords roc curve
#' #' @export
#' #' @examples
#' #' plot_auc(predictionObj, description, c = 'darkblue')
#' 
#' plot_auc <- function(predictionObj, description, c = 'darkblue') {
#'   plot(predictionObj[[1]], col = alpha(c,.2), main = description)  # plot the first round in full line color
#'   text(.5, .1, paste("mean AUC:", round(mean(unlist(  # add text with mean AUC value for all rounds
#'     lapply(predictionObj, function(x) {
#'       x$auc[[1]]
#'     })
#'   )), 2)))
#'   text(.5, 0, paste("SD:", round(sd(unlist(  # add text with std dev (if available)
#'     lapply(predictionObj, function(x) {
#'       x$auc[[1]]
#'     })
#'   )), 2)))
#'   if (length(predictionObj) > 1) {  # if there were multiple rounds, plot all AUC values in color with alpha = .2
#'     for (i in seq(2:length(predictionObj))) {
#'       lines(predictionObj[[i]], col = alpha(c, .2))
#'     }
#'   }
#' }
