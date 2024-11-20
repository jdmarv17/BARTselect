#' Function to generate global SE threshold
#'
#' @param record Matrix of CIPs from get_posteriors()
#' @param alpha Determines 1 - alpha quantile within threshold
#' @param use_means TRUE/FALSE whether to use variable level CMD means. If FALSE means are set to 0.
#'
#' @return A data frame with 1 - alpha quantiles, means, sd's, and global SE threshold for each variable
#' @importFrom stats quantile sd var
#' @importFrom dplyr mutate case_when
#' @keywords internal
#'
#' @examples
#' library(MASS)
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
#' bart_run = run_bart(formula, data, vars = vars)
#' 
#' posteriors = get_posteriors(bart_run[[1]], vars, parallel::detectCores())
#' thresholds = get_multiplier_thresh(record = posteriors[[1]], alpha = 0.01, use_means = TRUE)
#' 
get_multiplier_thresh = function(record, alpha = 0.01, use_means = TRUE) {

  # get means and sd's for each 
  if (use_means == TRUE) {
    means = colMeans(record, na.rm = TRUE)
  }else {
    means = rep(0, times = ncol(record))
  }
  std_dev = apply(record, MARGIN = 2, FUN = sd, na.rm = TRUE)
  
  
  # get list of variable specific quantiles
  quant_record = apply(record, MARGIN = 2, FUN = quantile, probs = (1 - alpha), na.rm = TRUE)
  
  # get multipliers
  stats = 
    data.frame(quant_record, means, std_dev) %>%
    mutate(C = (quant_record - means) / std_dev,
           C_thresh = (means + max(.data$C) * std_dev)) %>%
    mutate(C_thresh = case_when(
      C_thresh < 0 ~ 0, # if less than 0 set to 0
      C_thresh > 1 ~ 1, # if greater than 1 set to 1
      TRUE ~ C_thresh # otherwise leave the same
    ))
  
  # return C's and max{C_i}
  to_return = list(stats)
  return(to_return)
}
