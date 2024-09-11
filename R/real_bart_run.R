#' Function to run BART on the original data
#'
#' @param formula Object of class formula.
#' @param data A data frame of original data.
#' @param num_trees Number of BART trees per sample.
#' @param num_samps Number of posterior samples to take after burn-in.
#' @param num_burn Number of burn-in samples.
#' @param num_chains Number of independent MCMC chains to run.
#' @param num_thin Keep every 'num_thin' draw.
#' @param prior_power Power parameter for tree prior.
#' @param prior_base Base parameter for tree prior.
#' @param num_threads_bart Number of threads used to run BART - not recommended to be greater than 'num_chains'.
#' @param num_threads_wrangle Number of threads for post processing - can be larger than 'num_chains'.
#' @param vars Column names resulting from dbartsData(data).
#'
#' @return A list of four objects. The first three are the return objects from get_posteriors(), 
#' the fourth is a dataframe of filtered trees.
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
#' 
#' bart_run = real_bart_run(formula, data, vars = vars)
real_bart_run = function(formula, data,
                         num_trees = 10, num_samps = 5000,
                         num_burn = 2500, num_chains = parallel::detectCores(),
                         num_thin = 5, prior_power = 8, prior_base = .99,
                         num_threads_bart = min(num_chains, parallel::detectCores()),
                         num_threads_wrangle = parallel::detectCores(), vars) {
  
  # actual BART run
  bart_run = run_bart(formula, data = data, 
                      num_burn = num_burn, num_samps = num_samps,
                      num_chains = num_chains, num_trees = num_trees,
                      num_thin = num_thin, prior_power = prior_power,
                      prior_base = prior_base, num_threads_bart = num_threads_bart,
                      num_threads_wrangle = num_threads_wrangle, vars = vars)
  
  
  # process real bart inclusions to get posterior probabilities 
  posteriors = get_posteriors(bart_run[[1]], vars, num_threads_wrangle) # bart_run[[1]] is indicators df
  
  # list of return objects
  to_return = list(posteriors[[1]], posteriors[[2]], posteriors[[3]], bart_run[[2]])
  
  
  return(to_return)
  
}
