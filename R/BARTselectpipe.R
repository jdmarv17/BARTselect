#' Function to run automatic BART-based interaction and variable selection. Includes automatic GLM modeling
#' with the selected variables and interactions. If the dataset at hand includes factor covariates the modeling
#' may only include specific factor levels as covariates instead of the original factor covariate.
#'
#' @param formula Object of class formula.
#' @param data A data frame of original data.
#' @param num_trees Number of BART trees per sample.
#' @param num_samps Number of posterior samples to take after burn-in.
#' @param num_burn Number of burn-in samples.
#' @param num_thin Keep every 'num_thin' draw.
#' @param num_chains Number of independent MCMC chains to run.
#' @param alpha_g Alpha level for global threshold
#' @param alpha_l Alpha level for local thresholds
#' @param alpha_d Alpha level for global SE thresholds
#' @param alpha_g_vip Alpha level for variable selection threshold
#' @param num_null_run Number of null runs to do for threshold generations
#' @param method A value %in% c("global", "local", "both") to specify which selections are returned
#' @param prior_power Power parameter for tree prior.
#' @param prior_base Base parameter for tree prior.
#' @param num_threads_bart Number of threads used to run BART - not recommended to be greater than 'num_chains'.
#' @param num_threads_wrangle Number of threads for post processing - can be larger than 'num_chains'.
#' @param set_diff_thresh TRUE/FALSE whether user wants to manually set the CMD threshold
#' @param hierarchical TRUE/FALSE whether to force both main effects in model if corresponding pairwise interaction is included (only for interpretable model after selections).
#' @param response_type A value %in% c("continuous", "binary") to specify outcome type for interpretable model
#' @param diff_thresh A user supplied value for CMD threshold. Only used if set_diff_thresh = TRUE.
#'
#' @return A list of return objects. Object (1) is selected interactions and summary from global threshold. Object (2) is selected interactions
#' and summary from local thresholds. Objects (3)-(5) are the posterior CIPs, CMDs, and MIPs respectively. Objects (6)-(8)
#' are the global CIP threshold, local CIP thresholds, and CMD thresholds. Object (9) is a character specifying the threshold type used. Objects (10) and (11)
#' are selected variable indices and names from variable selection procedure. Objects (12) and (13) are the selected interactions from the global and local thresholds respectively.
#' Objects (14) and (15) are the fit interpretable models using global and local methods respectively.
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
#' 
#' select = BARTselectpipe(formula, data, alpha_g = 0.05, alpha_g_vip = 0.1)
#'  
BARTselectpipe = function(formula, data, 
                      num_trees = 10, num_samps = 5000,
                      num_burn = 5000,  num_thin = 5,
                      num_chains = min(4, dbarts::guessNumCores()),
                      alpha_g = 0.01, alpha_l = 0.01, alpha_d = 0.001,
                      alpha_g_vip = 0.01, num_null_run = 6, method = "global",
                      prior_power = 4, prior_base = .95,
                      num_threads_bart = min(dbarts::guessNumCores(), num_chains),
                      num_threads_wrangle = dbarts::guessNumCores(),
                      set_diff_thresh = FALSE, hierarchical = TRUE,
                      response_type = "continuous", diff_thresh) {
  
  `%ni%` = base::Negate(`%in%`)
  
  # grab the data frame bart used
  data_to_use = dbarts::dbartsData(formula, data)
  dbarts_x = janitor::clean_names(data_to_use@x)
  dbarts_y = data_to_use@y
  vars = colnames(dbarts_x)
  response = as.character(formula[2])
  
  # print progress
  if (set_diff_thresh == FALSE) {
    # get thresholds adaptively and print result
    thresh_list = get_thresholds(formula = formula, data = data, num_trees = num_trees,
                                 num_samps = num_samps, num_burn = num_burn, num_thin = num_thin,
                                 num_chains = num_chains, num_threads_bart = num_threads_bart,
                                 num_threads_wrangle = num_threads_wrangle,
                                 alpha_g = alpha_g, alpha_l = alpha_l, alpha_d = alpha_d,
                                 alpha_g_vip = alpha_g_vip, num_run = num_null_run,
                                 use_differences = !set_diff_thresh, prior_power = prior_power,
                                 prior_base = prior_base, vars = vars)
    cat("Difference thresholds found with global multiplier \n")
    cat("Average difference threshold =",  mean(thresh_list[[3]]), "\n")
  } else if (set_diff_thresh == TRUE) {
    # get thresholds but set differece threshold
    thresh_list = get_thresholds(formula = formula, data = data, num_trees = num_trees,
                                 num_samps = num_samps, num_burn = num_burn, num_thin = num_thin,
                                 num_chains = num_chains, num_threads_bart = num_threads_bart,
                                 num_threads_wrangle = num_threads_wrangle,
                                 alpha_g = alpha_g, alpha_l = alpha_l, alpha_d = alpha_d,
                                 alpha_g_vip = alpha_g_vip, num_run = num_null_run,
                                 use_differences = !set_diff_thresh, prior_power = prior_power,
                                 prior_base = prior_base, diff_thresh = diff_thresh,
                                 vars = vars)
    cat("Difference thresholds set:", diff_thresh, "\n")
  }
  
  cat("Variable selection VIP threshold found:", thresh_list[[4]], "\n")
  if (method == "global") {
    cat("Interaction selection CIP global threshold found:", thresh_list[[1]], "\n")
    cat("Starting real BART run \n")
  } else if (method == "local") {
    cat("Interaction selection CIP local thresholds found \n")
    cat("Starting real BART run \n")
  } else if (method == "both") {
    cat("Interaction selection CIP global threshold found:", thresh_list[[1]], "\n")
    cat("Interaction selection CIP local thresholds found \n")
    cat("Starting real BART run \n")
  }
  
  
  
  # get conditional co-inclusion df, conditional - marginal df, and marginals vector
  real_run = real_bart_run(formula = formula, data = data, num_trees = num_trees,
                           num_samps = num_samps, num_burn = num_burn, num_chains = num_chains,
                           num_thin = num_thin, prior_power = prior_power,
                           prior_base = prior_base, num_threads_bart = num_threads_bart,
                           num_threads_wrangle = num_threads_wrangle,
                           vars = vars)
  
  vip.real = get_vip(filtered_trees = real_run[[4]])
  
  
  # compare results to thresholds
  passed_thresh = test_thresholds(thresh_list[[1]], thresh_list[[2]], thresh_list[[3]], # these 3 are threshold args from null
                                  real_run[[1]],  real_run[[2]], real_run[[3]], # these 3 are posterior probs from real run
                                  method = method)
  # first return object is global selected, second is local selected
  
  
  
  if (method == "global" & nrow(passed_thresh[[1]]) != 0) {
    
    # get rid of duplicated rows if not empty
    selected_global = dplyr::distinct(passed_thresh[[1]], avg_inclusion_prob, .keep_all = TRUE)
    selected_local = NULL
    
    interact.global = tidyr::unite(selected_global, "pair", var1, var2, sep = ":")$pair
    interact.local = NULL
    
  } else if (method == "global" & nrow(passed_thresh[[1]]) == 0) {
    
    selected_global = passed_thresh[[1]]
    selected_local = NULL
    
    interact.global = c()
    interact.local = NULL
    
  } else if (method == "local" & nrow(passed_thresh[[1]]) != 0) {
    
    # get rid of duplicated rows if not empty
    selected_local = dplyr::distinct(passed_thresh[[1]], avg_inclusion_prob, .keep_all = TRUE)
    selected_global = NULL
    
    interact.global = NULL
    interact.local = tidyr::unite(selected_local, "pair", var1, var2, sep = ":")$pair
    
  } else if (method == "local" & nrow(passed_thresh[[1]]) == 0) {
    
    selected_local = passed_thresh[[1]]
    selected_global = NULL
    
    interact.global = NULL
    interact.local = c()
  } 
  
  
  # if using both methods do same thing but separately
  if (method == "both") {
    if (nrow(passed_thresh[[1]]) != 0) {
      
      # get rid of duplicated rows if not empty
      selected_global = dplyr::distinct(passed_thresh[[1]], avg_inclusion_prob, .keep_all = TRUE)
      interact.global = tidyr::unite(selected_global, "pair", var1, var2, sep = ":")$pair
      
    } else {
      
      selected_global = passed_thresh[[1]]
      interact.global = c()
      
    }
    
    if (nrow(passed_thresh[[2]]) != 0) {
      
      # get rid of duplicated rows if not empty
      selected_local = dplyr::distinct(passed_thresh[[2]], avg_inclusion_prob, .keep_all = TRUE)
      interact.local = tidyr::unite(selected_local, "pair", var1, var2, sep = ":")$pair
      
    } else {
      
      selected_local = passed_thresh[[2]]
      interact.local = c()
      
    }
    
  }
  
  
  # check which variables selected
  var.select = 
    vip.real %>%
    dplyr::mutate(select = ifelse(prop.split >= thresh_list[[4]], 1, 0)) %>%
    dplyr::filter(select == 1)
  
  var.names = vars[var.select$var]
  
  # now fit model using selected variables and interactions
  if (method == "global") {
    mod.global = fit_model(data = data, response_type = response_type, response = response,
                           main_effects = var.names, interactions = interact.global,
                           hierarchical = hierarchical,
                           dbarts_x = dbarts_x, dbarts_y = dbarts_y)
    mod.local = NULL
  } else if (method == "local") {
    mod.global = NULL
    mod.local = fit_model(data = data, response_type = response_type, response = response,
                          main_effects = var.names, interactions = interact.local,
                          hierarchical = hierarchical,
                          dbarts_x = dbarts_x, dbarts_y = dbarts_y)
  } else if (method == "both") {
    mod.global = fit_model(data = data, response_type = response_type, response = response,
                           main_effects = var.names, interactions = interact.global,
                           hierarchical = hierarchical,
                           dbarts_x = dbarts_x, dbarts_y = dbarts_y)
    mod.local = fit_model(data = data, response_type = response_type, response = response,
                          main_effects = var.names, interactions = interact.local,
                          hierarchical = hierarchical,
                          dbarts_x = dbarts_x, dbarts_y = dbarts_y)
  }
  
  
  
  # things to return
  to_return = list(selected_global, selected_local, 
                   real_run[[1]], real_run[[2]], real_run[[3]],
                   thresh_list[[1]], thresh_list[[2]], thresh_list[[3]],
                   method, var.select$var, var.names, interact.global, interact.local,
                   mod.global, mod.local)
  #,mod.global, mod.local)
  
  
  return(to_return)
  
}