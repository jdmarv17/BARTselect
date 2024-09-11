#' Function to get interaction and variable selection thresholds
#' 
#' @param formula Object of class formula.
#' @param data A data frame of original data.
#' @param num_trees Number of BART trees per sample.
#' @param num_samps Number of posterior samples to take after burn-in.
#' @param num_burn Number of burn-in samples.
#' @param num_thin Keep every 'num_thin' draw.
#' @param prior_power Power parameter for tree prior.
#' @param prior_base Base parameter for tree prior.
#' @param num_chains Number of independent MCMC chains to run.
#' @param num_threads_bart Number of threads used to run BART - not recommended to be greater than 'num_chains'.
#' @param num_threads_wrangle Number of threads for post processing - can be larger than 'num_chains'.
#' @param alpha_g Alpha level for global threshold
#' @param alpha_l Alpha level for local thresholds
#' @param alpha_d Alpha level for global SE thresholds
#' @param alpha_g_vip Alpha level for variable selection threshold
#' @param num_run Number of null runs to do for threshold generations
#' @param use_differences TRUE/FALSE whether global SE threshold should be used for CMDs. If FALSE user must supply diff_thresh 
#' @param use_means TRUE/FALSE whether to use means for CMD threshold. If false the variable level CMD means are set to 0.
#' @param diff_thresh A value for the difference threshold if use_differences set to FALSE. 
#' use_differences set to TRUE overrides a value of diff_thresh
#' @param vars Column names resulting from dbartsData(data).
#'
#' @return A list of 4 or 5 objects depending on whether use_differences = TRUE. First object is a global CIP threshold.
#' Second object is a vector of local CIP thresholds. Third object is vector of CMD thresholds. Fourth object is a global
#' variable selection threshold. Fifth object (if use_differences = TRUE) is the return object from get_multiplier_thresh
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
#' thresholds = get_thresholds(formula, data, vars = vars)
#' 
get_thresholds = function(formula, data,
                          num_trees = 10, num_samps = 5000,
                          num_burn = 5000, num_thin = 5,
                          prior_power = 4, prior_base = 0.99, 
                          num_chains = 4, num_threads_bart = min(num_chains, parallel::detectCores()),
                          num_threads_wrangle = parallel::detectCores(),
                          alpha_g = 0.01, alpha_l = 0.01, alpha_d = 0.01, alpha_g_vip = 0.01,
                          num_run = 6, use_differences = TRUE, use_means = TRUE, 
                          diff_thresh, vars) {
  
  
  if (use_differences == TRUE & use_means %ni% c(TRUE, FALSE)) {
    cat("Warning: 'use_means' missing or incorrectly specificied. Defaulting to use_means = TRUE.")
    use_means = TRUE
  }
  
  response = as.character(formula[2])
  
  # record to record vip's
  vip.record = matrix(ncol = length(vars), nrow = num_run)
  colnames(vip.record) = vars
  
  # initialize empty vectors to store results
  conditional_record = tibble()
  difference_record = c()
  
  # do the null runs and bind results after each run
  for (j in 1:num_run) { # start of j loop
    
    if (j != num_run) {
      cat("Starting null run", j, "of", num_run, "... \n")
    }
    if (j == num_run) {
      cat("Starting final null run. \n")
    }
    
    # create null data frame where y vector is randomized 
    null_response = sample(unlist(data[,response]), size = nrow(data), replace = FALSE)
    data[,response] = null_response
    
    # null BART run
    null_bart_run = run_bart(formula = formula, data = data, 
                             num_burn = num_burn, num_samps = num_samps,
                             num_chains = num_chains, num_trees = num_trees, 
                             num_thin = num_thin, prior_power = prior_power,
                             prior_base = prior_base, num_threads_bart = num_threads_bart,
                             num_threads_wrangle = num_threads_wrangle, vars = vars)
    
    # get VIPs in case we do var selection
    vip = get_vip(filtered_trees = null_bart_run[[2]])
    
    # This step accounts for potentially not seeing a variable in any tree 
    if (nrow(vip) == length(vars)) {
      vip.record[j,] = vip$prop.split
    } else {
      ind.excluded = which(1:length(vars) %ni% vip$var)
      ind.included = which(1:length(vars) %in% vip$var)
      vip.record[j,ind.excluded] = NA
      vip.record[j,ind.included] = vip$prop.split[ind.included]
    }
    
    # save indicators
    indicators = null_bart_run[[1]]
    
    # get null co-inclusion probabilities using get_posteriors
    posteriors = get_posteriors(indicators, vars, num_threads_wrangle) 
    
    conditionals_df = posteriors[[1]]
    diff_df = posteriors[[2]]
    marginals = posteriors[[3]]
    
    conditional_record = rbind(conditional_record, conditionals_df)
    difference_record = rbind(difference_record, diff_df)
    
    
  } # end of j loop
  
  # check for variables never selected and remove from threshold consideration
  NA.col.ind = unname(which(colSums(is.na(conditional_record)) == nrow(conditional_record)))
  non.NA.col.ind = which(1:length(vars) %ni% NA.col.ind)
  conditionals = conditional_record[,non.NA.col.ind]
  differences = difference_record[,non.NA.col.ind]
  
  # get difference threshold accounting for possibility of imputation
  if (use_differences == TRUE) {
    multipliers = get_multiplier_thresh(record = differences, 
                                        alpha = alpha_d, 
                                        use_means = TRUE)
    diff_thresholds_short = multipliers[[1]]$C_thresh
    
    if (length(NA.col.ind) > 0) { # if imputation needed 
      impute_val = mean(diff_thresholds_short) # get impute value
      diff_thresholds = rep(0, times = length(vars)) # create new vector to store values
      diff_thresholds[non.NA.col.ind] = diff_thresholds_short # place calculated values in normal positions
      diff_thresholds[NA.col.ind] = impute_val # put impute value in correct place
    } else {
      diff_thresholds = diff_thresholds_short
    }
    
  } else { # if not using global SE method then set all diff_thresholds to the user supplied value
    diff_thresholds = rep(diff_thresh, times = ncol(conditionals))
  }
  
  
  # get local threshold accounting for possibility of imputation
  if (length(NA.col.ind) > 0) {
    local_thresholds_short = apply(conditionals, MARGIN = 2,
                                   FUN = quantile, probs = (1 - alpha_l), na.rm = TRUE)
    
    impute_val = mean(local_thresholds_short) # get impute value
    local_thresholds = rep(0, times = length(vars)) # create new vector to store values
    local_thresholds[non.NA.col.ind] = local_thresholds_short # place calculated values in normal positions
    local_thresholds[NA.col.ind] = impute_val # put impute value in correct place
  } else {
    local_thresholds = apply(conditionals, MARGIN = 2,
                             FUN = quantile, probs = (1 - alpha_l), na.rm = TRUE)
  }
  
  
  # get global threshold
  global_threshold =
    unname(quantile(conditionals, probs = (1 - alpha_g), na.rm = TRUE))
  
  # variable selection threshold
  vip.global =
    unname(quantile(vip.record, probs = (1 - alpha_g_vip), na.rm = TRUE))
  
  
  # list things to return
  if (use_differences == TRUE) {
    to_return = list(global_threshold, local_thresholds, diff_thresholds, vip.global, multipliers)
    
  } else { # if not using global SE method then don't return multipliers either
    to_return = list(global_threshold, local_thresholds, diff_thresholds, vip.global)
  }  
  
  return(to_return)
}