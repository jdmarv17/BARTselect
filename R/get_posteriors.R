#' Function to get posterior MIPs, CIPs, and CMDs
#'
#' @param indicators A data frame of indicators. This data frame is output by run_bart and has one column
#' for each BART covariate and one row for each posterior tree.
#' @param vars Column names resulting from dbartsData(data).
#' @param num_threads_wrangle Number of threads for processing - can be larger than 'num_chains'
#'
#' @return A list with 3 objects. First is the matrix of CIPs, second is matrix of CMDs, third is vector of MIPs
#'
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
#' 
get_posteriors = function(indicators, vars, num_threads_wrangle) {

  # functions to help replace loops in get_posteriors()
  co.include = function(index, df) { colMeans(subset(df, df[, index] == 1)) }
  diff.include = function(index, record.df, margin.df) { (record.df[index, ] - margin.df[index]) }
  
  # grab marginal inclusions for difference matrix
  marginals = colMeans(indicators)
  
  # going over each predictor, get co-inclusion probability of that predictor in all other predictors' trees
  record = matrix(unlist(
    parallel::mclapply(X = 1:ncol(indicators), FUN = co.include, df = indicators, mc.cores = num_threads_wrangle) 
  ), ncol = ncol(indicators), nrow = ncol(indicators))
  
  colnames(record) = vars
  diag(record) = NA
  
  # conditional - marginal matrix
  diff_matrix = matrix(unlist(
    parallel::mclapply(X = 1:ncol(indicators), FUN = diff.include,
                       record.df = record, margin.df = marginals, mc.cores = num_threads_wrangle)
  ), ncol = ncol(indicators), nrow = ncol(indicators)) # this puts d_{i|j} = p_{i|j} - p_{i} in i-th col
  
  # fix names
  colnames(diff_matrix) = vars
  
  # create return list
  to_return = list(record, diff_matrix, marginals)
  
  return(to_return)
}