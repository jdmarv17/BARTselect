#' Function to test posterior CIPs, CMDs, MIPs against supplied thresholds
#'
#' @param global_thresh A global CIP threshold value
#' @param local_thresholds A vector of local CIP thresholds - must be same length as marginals
#' @param diff_thresholds A vector of CMD thresholds - must be same length as marginals
#' @param conditional_df A data frame or matrix of posterior CIPs
#' @param diff_df A data frame or matrix of posterior CMDs
#' @param marginals A vector of posterior MIPs
#' @param method A value %in% c("global", "local", "both") to specify which selections are returned
#'
#' @return A list of selected interactions. If method %in% c("global", "local") the list will contain one data frame 
#' of selected interactions and associated posterior info. If method = "both" a data frame will be returned for each method
#' 
#'
#' @examples
#' #' library(MASS)
#' library(dbarts)
#' library(dplyr)
#' library(janitor)
#' data(birthwt)
#' data = birthwt %>% dplyr::select(., -low)
#' formula = bwt ~ .
#' dbarts_data = dbartsData(formula, data)
#' dbarts_x = clean_names(dbarts_data@x)
#' dbarts_y = dbarts_data@y
#' vars = colnames(dbarts_x)
#' 
#' thresholds = get_thresholds(formula, data, vars = vars)
#' bart_run = real_bart_run(formula, data, vars = vars)
#' 
#' select = test_thresholds(thresholds[[1]], thresholds[[2]], thresholds[[3]],
#'                          bart_run[[1]], bart_run[[2]], bart_run[[3]],
#'                          method = "global")
#'                            
test_thresholds = function(global_thresh, local_thresholds, diff_thresholds,
                           conditional_df, diff_df, marginals, method = "global") {
  
  # initialize records
  global_record = tibble()
  local_record = tibble()
  
  # get matrix form for easier indexing
  conditional_mat = as.matrix(conditional_df)
  diff_mat = as.matrix(diff_df)
  
  # go over real run posterior probabilities to test thresholds
  for (i in 1:ncol(conditional_mat)) { # go across cols
    for (j in 1:nrow(conditional_mat)) { # go down rows
      
      # conditional threshold check and conditional - marginal check
      if (!is.na(conditional_mat[i, j]) & # make sure we aren't on diagonal 
          !is.na(conditional_mat[j, i])) {
        
        if (conditional_mat[i, j] > global_thresh & # P(x_i in x_j tree) > conditional_thresh
            conditional_mat[j, i] > global_thresh & # P(x_j in x_i tree) > conditional_thresh
            diff_mat[i, j] > diff_thresholds[i] & # make sure both vars have conditional - marginal > threshold
            diff_mat[j, i] > diff_thresholds[j] ) {
          # get some summary stats
          new_row = tibble(var1 = colnames(conditional_mat)[i],
                           var2 = colnames(conditional_mat)[j],
                           var1_in_var2_trees = conditional_mat[i, j],
                           var2_in_var1_trees = conditional_mat[j, i],
                           var1_diff = diff_mat[i, j],
                           var2_diff = diff_mat[j, i],
                           marg_var1 = marginals[i],
                           marg_var2 = marginals[j],
                           var1_thresh = local_thresholds[i], #local_thresholds$var_maxes[i],
                           var2_thresh = local_thresholds[j], #local_thresholds$var_maxes[j],
                           avg_inclusion_prob = ((conditional_mat[i, j] + conditional_mat[j, i]) / 2),
                           conditional_threshold = global_thresh)
          
          # bind back
          global_record =
            bind_rows(global_record, new_row)
          
          
        } # end of global threshold check
      }
      
      
      
      # if both pass each others variable specific threshold and neither is NA
      if (!is.na(conditional_mat[i, j]) & # make sure we aren't on diagonal 
          !is.na(conditional_mat[j, i])) {
        if (conditional_mat[i,j] > local_thresholds[i] & #local_thresholds$var_maxes[i] & # P(x_i in x_j tree) > max(P(x_i in x_j tree) i != j under null)
            conditional_mat[j,i] > local_thresholds[j] & #local_thresholds$var_maxes[j] & # P(x_j in x_i tree) > max(P(x_j in x_i tree) j != i under null)
            diff_mat[i, j] > diff_thresholds[i] & # make sure both vars have conditional - marginal > 0
            diff_mat[j, i] > diff_thresholds[j]) {
          # get average co-inclusion prob
          new_row = tibble(var1 = colnames(conditional_mat)[i],
                           var2 = colnames(conditional_mat)[j],
                           var1_in_var2_trees = conditional_mat[i, j],
                           var2_in_var1_trees = conditional_mat[j, i],
                           var1_diff = diff_mat[i, j],
                           var2_diff = diff_mat[j, i],
                           marg_var1 = marginals[i],
                           marg_var2 = marginals[j],
                           var1_thresh = local_thresholds[i], #local_thresholds$var_maxes[i],
                           var2_thresh = local_thresholds[j], #local_thresholds$var_maxes[j],
                           avg_inclusion_prob = ((conditional_mat[i, j] + conditional_mat[j, i]) / 2),
                           conditional_threshold = global_thresh)
          
          # bind back
          local_record =
            bind_rows(local_record, new_row)
        }
        
      }
      
      
    } # end of j index
  } # end of i index
  
  if (method == "global") {
    to_return = list(global_record)
  } else if (method == "local") {
    to_return = list(local_record)
  } else if (method == "both") {
    to_return = list(global_record, local_record)
  } else {
    cat("Invalid thresholding method supplied. \n")
  }
  
  return(to_return)
}