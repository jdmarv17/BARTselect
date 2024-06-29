library(dbarts)
library(tidyverse)
library(parallel)
library(janitor)

##### VETTED FUNCTIONS AFTER 7/18/23 #####
`%ni%` = Negate(`%in%`)





fit.model = function(data, dbarts_x, dbarts_y, response.type = "continuous", response,
                     main.effects = c(), interactions = c(), hierarchical = TRUE) {
  
  # This function fits lm or glm with passed selected main effects and interactions
  # where it assumes interactions are passed in colon form
  ##### Return 1: a tibble of estimated coefficients
  ##### Return 2: lm/glm object using passed main effects and interactions
  # ARGS:
  # `data` = data used to select (same used to fit)
  # `response.type` = binary or continuous response
  # `response` = variable name of response (character valued)
  # `main.effects` = vector of predictors (character valued)
  # `interactions` = vector of interactions (character valued, w/ colon e.g. c("x1:x3", "x5:x9"))
  # `hierarchical` = TRUE/FALSE should both main effects be included if interaction is included 
  
  # do we have interactions and main effects?
  nonzero.int = ifelse(length(interactions) > 0, TRUE, FALSE)
  nonzero.var = ifelse(length(main.effects) > 0, TRUE, FALSE)
  
  # grab outcome to attach to dbarts_data
  outcome = unname(unlist(data[,response]))
  model.data = as.data.frame(cbind(dbarts_x, dbarts_y))
  colnames(model.data) = c(colnames(dbarts_x), response)
  
  
  if (hierarchical) {
    
    if (nonzero.int & nonzero.var) { # detected both main effects and interactions
      # replace ":" with "*" in interactions before formula
      interactions = str_replace_all(string = interactions, 
                                     pattern = ":",
                                     replacement = "*")
      
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main.effects, collapse = " + "), " + ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (nonzero.int & !nonzero.var) { # detected only interactions
      # replace ":" with "*" in interactions before formula
      interactions = str_replace_all(string = interactions, 
                                     pattern = ":",
                                     replacement = "*")
      
      formula = formula(paste(as.name(response), " ~ ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & nonzero.var) { # detected only main effects
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main.effects, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & !nonzero.var) { # intercept only model
      formula = formula(paste(as.name(response), " ~ ",
                              1,
                              sep = ""))
    }
    
  } else {
    if (nonzero.int & nonzero.var) { # detected both main effects and interactions
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main.effects, collapse = " + "), " + ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (nonzero.int & !nonzero.var) { # detected only interactions
      formula = formula(paste(as.name(response), " ~ ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & nonzero.var) { # detected only main effects
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main.effects, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & !nonzero.var) { # intercept only model
      formula = formula(paste(as.name(response), " ~ ",
                              1,
                              sep = ""))
    }
  }
  
  
  # depending on response type fit a linear regression or logistic regression
  if (response.type == "continuous") {
    mod = lm(formula = formula, data = model.data)
    coefs = coef(mod)
  } else if (response.type == "binary") {
    mod = glm(formula = formula, data = model.data, family = binomial(link = "logit"))
    coefs = coef(mod)
  } else {
    cat("Expected response.type = 'continuous' or response.type = 'binary'.\n")
    cat("Other response.type values are invalid.\n")
  }
  
  coef.tibble = rownames_to_column(data.frame(coefs), "var")
  to.return = list(coef.tibble, mod)
  
  return(to.return)
}

# This function fits lm or glm with passed selected main effects and interactions
fit.model.sim = function(data, response.type = "continuous", response,
                     main.effects, interactions, intercept = TRUE) {
  # ARGS:
  # `original.formula` = formula for original selection (this tells us candidate covariates)
  # `data` = data used to select (same used to fit)
  # `response.type` = binary or continuous response
  # `response` = variable name of response (character valued)
  # `selected.vars` = index of selected variables from candidate covariates
  # `selected.interact` = pairswise interaction vector e.g. c("x1:x3", "x5:x9")
  NA.int = ifelse(NA %in% interactions, TRUE, FALSE)
  NA.var = ifelse(NA %in% main.effects, TRUE, FALSE)
  
  nonzero.int = ifelse(length(interactions) > 0, TRUE, FALSE)
  nonzero.var = ifelse(length(main.effects) > 0, TRUE, FALSE)
  
  if (intercept == TRUE) {
    if (!NA.int & nonzero.int & !NA.var & nonzero.var) {
      # detected both main effects and interactions
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main.effects, collapse = " + "), " + ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (!NA.int & nonzero.int & (NA.var | !nonzero.var)) {
      # detected only interactions
      formula = formula(paste(as.name(response), " ~ ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if ((NA.int | !nonzero.int) & !NA.var & nonzero.var) {
      # detected only main effects
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main.effects, collapse = " + "),
                              sep = ""))
      
    } else if ((NA.int | !nonzero.int) & (NA.var | !nonzero.var)) {
      formula = formula(paste(as.name(response), " ~ ",
                              1,
                              sep = ""))
    }
  } else if (intercept == FALSE) {
    if (!NA.int & nonzero.int & !NA.var & nonzero.var) {
      # detected both main effects and interactions
      formula = formula(paste(as.name(response), " ~ ", "-1", "+",
                              paste(main.effects, collapse = " + "), " + ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (!NA.int & nonzero.int & (NA.var | !nonzero.var)) {
      # detected only interactions
      formula = formula(paste(as.name(response), " ~ ", "-1", "+",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if ((NA.int | !nonzero.int) & !NA.var & nonzero.var) {
      # detected only main effects
      formula = formula(paste(as.name(response), " ~ ", "-1", "+",
                              paste(main.effects, collapse = " + "),
                              sep = ""))
      
    } else if ((NA.int | !nonzero.int) & (NA.var | !nonzero.var)) {
      formula = formula(paste(as.name(response), " ~ ",
                              1,
                              sep = ""))
    }
    
  }
  
  
  
  
  # depending on response type fit a linear regression or logistic regression
  if (response.type == "continuous") {
    mod = lm(formula = formula, data = data)
    coefs = coef(mod)
  } else if (response.type == "binary") {
    mod = glm(formula = formula, data = data, family = binomial(link = "logit"))
    coefs = coef(mod)
  } else {
    cat("Expected response.type = 'continuous' or response.type = 'binary'.\n")
    cat("Other response.type values are invalid.\n")
  }
  
  coef.tibble = rownames_to_column(data.frame(coefs), "var")
  to.return = list(coef.tibble, mod)
  
  return(to.return)
}




# Mostly for simulations, requires known true pairs
process.select.int = function(selected, version, correct.pairs) {
  # ARGS:
  # `selected` = vector of selected interactions
  # `version` = what method was used to select interactions
  # `correct.pairs` = what true interactions were used to generate data
  if (length(selected) != 0) {
    detected = 
      tibble(pair = selected) %>%
      separate(., pair, sep = ":", into = c("var1", "var2")) %>%
      unite(., "pair", var1, var2, sep = ":", remove = FALSE) %>%
      unite(., "pair2", var2, var1, sep = ":", remove = TRUE) %>%
      mutate(true_interaction = case_when(
        pair %in% correct.pairs ~ 1,
        pair2 %in% correct.pairs ~ 1,
        TRUE ~ 0 ),
        version = version) %>%
      dplyr::select(., -pair2)
  } else {
    detected = tibble(pair = NA, 
                      true_interaction = 0,
                      version = version)
  }
  
  return(detected)
  
}

# mostly for simulations
process.select.vars = function(selected, version, correct.vars) {
  
  if (length(selected) != 0) {
    selected.vars = 
      data.frame(var = selected) %>%
      mutate(correct = ifelse(var %in% correct.vars, 1, 0)) %>%
      mutate(version = version) 
  } else {
    selected.vars = tibble(var = NA, 
                           correct = 0,
                           version = version)
  }
  
  return(selected.vars)
  
}



make.indicators = function(var.indices, p) {

  ##### This function takes dimension p and a vector of integers and
  ##### return a vector of 0's and 1's.
  ##### The 1's are in the position equal to the value of the integers in the vars argument
  ##### ex: p = 10, vars = c(1, 2, 4, 8) --> return = c(1, 1, 0, 1, 0, 0, 0, 1, 0, 0)
  ##### Return 1: ind_vec is the indicator vector
  # ARGS:
  # `p` = dimension of return vector, corresponds with the highest value vars could take (not observed)
  # `var.indices` = vector of integers, no longer than p, with  1 < integers < p


  # get unique, sorted list and initialize the vector that will be returned
  ordered = unique(sort(var.indices))
  ind_vec = c()
  for (i in 1:length(ordered)) {
    if (i == 1) { # if its the first in the list
      next_vec = c(rep(0, times = (ordered[i] - 1)), 1)
      # if var_1 = 5, this makes a vector with four 0's and one 1

    } else { # otherwise need to repeat based on difference between last var
      next_vec = c(rep(0, times = (ordered[i] - ordered[i - 1] - 1)), 1)
      # this makes a vector with a series of 0's (possibly none), followed by a 1
    }

    ind_vec = c(ind_vec, next_vec) # bind vectors together to keep running list
    # this will have all the 1's we need after the loop is done
  }

  # now check if we need to append 0's to the end
  # if last var = p then we are done, otherwise append 0's
  if (ordered[length(ordered)] == p) {
    return(ind_vec)

  } else {
    to_append = p - ordered[length(ordered)]
    ind_vec = c(ind_vec, rep(0, times = to_append))
    return(ind_vec)
  }

}




run_bart = function(formula, data, 
                     num_trees = 10, num_samps = 5000,
                     num_burn = 5000, num_thin = 5,
                     prior_power = 4, prior_base = 0.99,
                     num_chains = 4, num_threads_bart = min(num_chains, parallel::detectCores()),
                     num_threads_wrangle = parallel::detectCores(), vars) {
  
  ##### This function will take BART params and returns two objects:
  ##### Return 1: a data frame of variable inclusion indicators
  #####           i.e. one row for each tree with indicator = 1 if var in the tree
  ##### Return 2: filtered trees from BART run
  # ARGS:
  # `formula` = regression formula
  # `data` = data frame with x's and y's
  # `num_trees` = param for number of BART trees
  # `num_samps` = number of MCMC samples
  # `num_burn` = number of MCMC samples to burn
  # `num_thin` = keep every 'num_thin'th sample
  # `prior_power` = power value for prior governing tree depth
  # `prior_base` = base value for prior governing tree depth
  # `num_chains` = number of MCMC chains
  # `num_threads` = number of cores to use (should be the same as num_chains for efficiency)
  
  # calculate number of samples after thinning to 
  actual_samps = num_samps / num_thin
  
  # if n <= p get a sigma estimate from std dev of the response, otherwise let bart2() estimate it
  sigma = ifelse(nrow(data) <= (ncol(data) - 1), sd(data[[response]], na.rm = TRUE), NA_real_)
  # ^ THIS CHECK PROBABLY NEEDS TO BE ADJUSTED IF USING DUMMY MATRIX NOW
  
  # run bart2()
  bart_run = bart2(formula, data = data, n.chains = num_chains, n.burn = num_burn,
                   power = prior_power, base = prior_base, n.trees = num_trees,
                   n.samples = num_samps, keepTrees = TRUE, printCutoffs = 0,
                   verbose = FALSE, n.thin = num_thin, sigest = sigma,
                   n.threads = num_threads_bart)
  
  
  # extract the trees
  bart_trees = dbarts::extract(bart_run, type = c("trees"))
  
  # get rid of terminal nodes (THIS IS FIRST RETURN OBJECT)
  filtered_trees = bart_trees[bart_trees$var != -1,]
  
  # split up by chain, tree, sample to take advantage of mclapply
  var.list = 
    filtered_trees %>%
    group_by(., chain, tree, sample) %>%
    group_split() %>% # list of filtered_tree tibbles by chain, tree, sample
    mclapply(., function(x) x[,"var"], mc.cores = num_threads_wrangle) %>% # select only var column
    lapply(., unlist) %>% # get as vector
    lapply(., unname) # get rid of names
  
  # run through make.indicators() with mclapply
  indicator.list = mclapply(X = var.list, FUN = make.indicators,
                            p = length(vars), mc.cores = num_threads_wrangle)
  
  # list to return
  to.return = list(as.data.frame(data.table::transpose(indicator.list), col.names = vars), filtered_trees)
  
  return(to.return)
}






get_posteriors = function(indicators, vars, num_threads_wrangle) {
  #### for notation let x_j be col name and x_i row
  
  ##### This function will take formula, data, and indicators df as arguments
  ##### and return posterior probabilities for conditional co-inclusion, 
  ##### conditional - marginal, and marginal
  ##### Return 1: a matrix that has conditional co-inclusion probabilities
  #####           with form x_i in x_j trees
  ##### Return 2: a matrix that has conditional - marginal with 
  #####           form x_i in x_j trees - marginal inclusion of x_i
  ##### Return 3: a vector of marginal inclusion probabilities
  # ARGS:
  # `formula` = regression formula
  # `data` = data frame with x's and y's
  # `indicators` = data frame of tree inclusion indicators
  
  # functions to help replace loops in get_posteriors()
  co.include = function(index, df) { colMeans(subset(df, df[, index] == 1)) }
  diff.include = function(index, record.df, margin.df) { (record.df[index, ] - margin.df[index]) }
  
  # grab marginal inclusions for difference matrix
  marginals = colMeans(indicators)
  
  # going over each predictor, get co-inclusion probability of that predictor in all other predictors' trees
  record = matrix(unlist(
    mclapply(X = 1:ncol(indicators), FUN = co.include, df = indicators, mc.cores = num_threads_wrangle) 
  ), ncol = ncol(indicators), nrow = ncol(indicators))
  
  colnames(record) = vars
  diag(record) = NA
  
  # conditional - marginal matrix
  diff_matrix = matrix(unlist(
    mclapply(X = 1:ncol(indicators), FUN = diff.include,
             record.df = record, margin.df = marginals, mc.cores = num_threads_wrangle)
  ), ncol = ncol(indicators), nrow = ncol(indicators)) # this puts d_{i|j} = p_{i|j} - p_{i} in i-th col
  
  # fix names
  colnames(diff_matrix) = vars
  
  # create return list
  to_return = list(record, diff_matrix, marginals)
  
  return(to_return)
}







get_multiplier_thresh = function(record, alpha = 0.01, use_means = TRUE) {
  ##### This function will take a probability record and an alpha and return a global 
  ##### multiplier C* according to Bleich et al. (2014) global SE strategy
  ##### 
  ##### Return 1: df with 1 - alpha quantiles, mean, sd, and multiplier thresh for each variable
  # ARGS:
  # `record` = record of null run probabilities (conditional or conditional - marginal)
  #          = contains "p" columns and num_run * p  rows
  # `alpha` = determines quantile in multiplier method
  # `use_means` = a boolean to determine whether to use means or set to 0
  
  
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
           C_thresh = (means + max(C) * std_dev)) %>%
    mutate(C_thresh = case_when(
      C_thresh < 0 ~ 0, # if less than 0 set to 0
      C_thresh > 1 ~ 1, # if greater than 1 set to 1
      TRUE ~ C_thresh # otherwise leave the same
    ))
  
  
  # return C's and max{C_i}
  to_return = list(stats)
  return(to_return)
}


get_thresholds = function(formula, data,
                          num_trees = 10, num_samps = 5000,
                          num_burn = 5000, num_thin = 5,
                          prior_power = 4, prior_base = 0.99, 
                          num_chains = 4, num_threads_bart = min(num_chains, parallel::detectCores()),
                          num_threads_wrangle = parallel::detectCores(),
                          alpha_g = 0.01, alpha_l = 0.01, alpha_d = 0.01, alpha.g.vip = 0.01,
                          num_run = 6, use_differences = TRUE, use_means = TRUE, 
                          diff.thresh, vars) {
  
  #### for notation let x_j be col name and x_i row
  
  ##### This function will take formula, data, BART params, alpha, and number of runs and 
  ##### returns several thresholds that will be used to evaluate posterior probabilities 
  ##### from a real BART run
  ##### Return 1: 1 - alpha_g conditional co-inclusion global threshold
  ##### Return 2: 1 - alpha_l conditional co-inclusion local threshold
  ##### Return 3: vector of difference thresholds using 1 - alpha_d in multiplier method
  # ARGS:
  # `formula` = regression formula
  # `data` = data frame with x's and y's
  # `num_trees` = param for number of BART trees
  # `num_samps` = number of MCMC samples
  # `num_burn` = number of MCMC samples to burn
  # `num_thin` = keep every 'num_thin'th sample
  # `prior_power` = power value for prior governing tree depth
  # `prior_base` = base value for prior governing tree depth
  # `num_chains` = number of MCMC chains
  # `num_threads` = number of cores to use (should be the same as num_chains for efficiency)
  # `alpha_g` = determines 1 - alpha_g quantile for global threshold
  # `alpha_l` = determines 1 - alpha_l quantile for local thresholds
  # `alpha_d` = determines 1 - alpha_d quantile for difference threshold
  # `num_run` = number of null runs
  # `use_differences` = logical whether to use the difference threshold, if FALSE
  #                     must supply `diff.thresh` value
  # `use_means` = logical whether to use the true difference means when finding the difference threshold
  # `diff.thresh` = optional user supplied difference threshold if not using data generated threshold
  
  
  
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
    vip = get.vip.pipe(filtered.trees = null_bart_run[[2]])
    
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
    diff_thresholds = rep(diff.thresh, times = ncol(conditionals))
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
    unname(quantile(vip.record, probs = (1 - alpha.g.vip), na.rm = TRUE))
  
  
  # list things to return
  if (use_differences == TRUE) {
    to_return = list(global_threshold, local_thresholds, diff_thresholds, vip.global, multipliers)
    
  } else { # if not using global SE method then don't return multipliers either
    to_return = list(global_threshold, local_thresholds, diff_thresholds, vip.global)
  }  
  
  return(to_return)
}




real_BART_run = function(formula, data,
                          num_trees = 10, num_samps = 5000,
                          num_burn = 2500, num_chains = parallel::detectCores(),
                          num_thin = 5, prior_power = 8, prior_base = .99,
                          num_threads_bart = min(num_chains, parallel::detectCores()),
                          num_threads_wrangle = parallel::detectCores(), vars) {
  
  #### for notation let x_j be col name and x_i row
  
  ##### This function will take formula, data, and BART params and
  ##### return several posterior probability df's
  ##### Return 1: posterior conditional co-inclusion df 
  ##### Return 2: posterior conditional co-inclusion minus marginal inclusion df
  ##### Return 3: posterior marginal inclusion probabilities
  # ARGS:
  # `formula` = regression formula
  # `data` = data frame with x's and y's
  # `num_trees` = param for number of BART trees
  # `num_samps` = number of MCMC samples
  # `num_burn` = number of MCMC samples to burn
  # `num_chains` = number of chains for MCMC
  # `prior_power` = power value for prior governing tree depth
  # `prior_base` = base value for prior governing tree depth
  
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



test_thresholds = function(global_thresh, local_thresholds, diff_thresholds,
                           conditional_df, diff_df, marginals, method = "global") {
  
  
  ##### This function will take conditional threshold, difference threshold, variable
  ##### specific thresholds, conditional co-inclusion df, difference df, and marginal 
  ##### inclusions and return pairs of variables that pass either the global threshold
  ##### and have conditional > marginal for both vars or the local threshold and have 
  ##### conditional > marginal for both vars
  ##### Return 1: pairs passing global threshold
  ##### Return 2: pairs passing local threshold
  # ARGS:
  # `global thresh` = 1 - alpha quantile of conditional inclusion probs from null
  # `local_thresholds` = 1 - alpha quantile of each var's conditional inclusion probs from null
  # `diff_thresholds` = conditional - marginal thresholds (either multiplier method or set)
  # `conditional_df` = conditional inclusion probs from real run
  # `diff_df` = conditional - marginal probs from real run
  # `marginals` = marginal inclusion probs from real run
  
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





## Variable selection functions
get.vip = function(bart.run) {
  # extract trees and filter out terminal nodes
  bart.trees = dbarts::extract(bart.run, type = c("trees"))
  
  # get rid of terminal nodes (THIS IS FIRST RETURN OBJECT)
  filtered.trees = bart.trees[bart.trees$var != -1,]
  num.splits = nrow(filtered.trees)
  vip = 
    filtered.trees %>% 
    group_by(var) %>% 
    summarise(prop.split = n() / num.splits)
  
  return(vip)
}



get.vip.pipe = function(filtered.trees) {

  # get vips
  num.splits = nrow(filtered.trees)
  vip = 
    as.data.frame(filtered.trees) %>% 
    group_by(var) %>% 
    summarise(prop.split = n() / num.splits)
  
  return(vip)
}






BARTselect = function(formula, data, 
                  num_trees = 10, num_samps = 5000,
                  num_burn = 5000,  num_thin = 5,
                  num_chains = dbarts::guessNumCores(),
                  alpha_g = 0.01, alpha_l = 0.01, alpha_d = 0.001,
                  alpha.g.vip = 0.01, num_null_run = 6, method = "global",
                  prior_power = 5, prior_base = .99,
                  num_threads_bart = min(parallel::detectCores(), num_chains),
                  num_threads_wrangle = parallel::detectCores(),
                  set.diff.thresh = FALSE, hierarchical = TRUE,
                  response.type = "continuous", diff.thresh) {
  
  ##### MAIN FUNCTION: This function combines all previous functions to return selected
  #####                interactions according to two methods: 1. global threshold + difference threshold
  #####                                                       2. local thresholds + difference threshold
  ##### Return 1: selected interactions and summary statistics (Global threshold)
  ##### Return 2: selected interactions and summary statistics (Local method)
  ##### Return 3: a df with variable specific conditional co-inclusion thresholds (Local thresholds)
  ##### Return 4-6: posterior probabilities from real_run (CAN GET RID OF THESE AT SOME POINT - JUST DIAGNOSTICS)
  ##### Return 7: local (1- \alpha) thresholds for each covariate generated from null runs
  
  ##### ^^ THESE NEED RENUMBERING ^^
  
  # ARGS:
  # `formula` = regression formula
  # `data` = data frame with x's and y's
  # `num_trees` = param for number of BART trees
  # `num_samps` = number of MCMC samples
  # `num_burn` = number of MCMC samples to burn
  # `num_chains` = number of chains for MCMC
  # `alpha` = determines 1 - alpha quantile for thresholds
  # `num_run` = number of null runs
  # `prior_power` = power value for prior governing tree depth
  # `prior_base` = base value for prior governing tree depth
  
  
  # grab the data frame bart used
  data_to_use = dbartsData(formula, data)
  dbarts_x = clean_names(data_to_use@x)
  dbarts_y = data_to_use@y
  vars = colnames(dbarts_x)
  response = as.character(formula[2])
  
  # print progress
  if (set.diff.thresh == FALSE) {
    # get thresholds adaptively and print result
    thresh_list = get_thresholds(formula = formula, data = data, num_trees = num_trees,
                                 num_samps = num_samps, num_burn = num_burn, num_thin = num_thin,
                                 num_chains = num_chains, num_threads_bart = num_threads_bart,
                                 num_threads_wrangle = num_threads_wrangle,
                                 alpha_g = alpha_g, alpha_l = alpha_l, alpha_d = alpha_d,
                                 alpha.g.vip = alpha.g.vip, num_run = num_null_run,
                                 use_differences = !set.diff.thresh, prior_power = prior_power,
                                 prior_base = prior_base, vars = vars)
    cat(1 - alpha_d, "difference thresholds found with global multiplier \n")
    cat("Average difference threshold =",  mean(thresh_list[[3]]), "\n")
  } else if (set.diff.thresh == TRUE) {
    # get thresholds but set differece threshold
    thresh_list = get_thresholds(formula = formula, data = data, num_trees = num_trees,
                                 num_samps = num_samps, num_burn = num_burn, num_thin = num_thin,
                                 num_chains = num_chains, num_threads_bart = num_threads_bart,
                                 num_threads_wrangle = num_threads_wrangle,
                                 alpha_g = alpha_g, alpha_l = alpha_l, alpha_d = alpha_d,
                                 alpha.g.vip = alpha.g.vip, num_run = num_null_run,
                                 use_differences = !set.diff.thresh, prior_power = prior_power,
                                 prior_base = prior_base, diff.thresh = diff.thresh,
                                 vars = vars)
    cat("Difference thresholds set:", diff.thresh, "\n")
  }
  
  cat("Variable selection VIP global threshold found:", thresh_list[[4]], "\n")
  if (method == "global") {
    cat(1 - alpha_g, "global threshold found:", thresh_list[[1]], ", starting real BART run \n")
  } else if (method == "local") {
    cat(1 - alpha_l, "variable specific local thresholds found, starting real BART run \n")
  } else if (method == "both") {
    cat(1 - alpha_g, "global threshold found:", thresh_list[[1]], "\n")
    cat(1 - alpha_l, "variable specific local thresholds found, starting real BART run \n")
  }
  
  
  
  # get conditional co-inclusion df, conditional - marginal df, and marginals vector
  real_run = real_BART_run(formula = formula, data = data, num_trees = num_trees,
                           num_samps = num_samps, num_burn = num_burn, num_chains = num_chains,
                           num_thin = num_thin, prior_power = prior_power,
                           prior_base = prior_base, num_threads_bart = num_threads_bart,
                           num_threads_wrangle = num_threads_wrangle,
                           vars = vars)
  
  vip.real = get.vip.pipe(filtered.trees = real_run[[4]])
  
  
  # compare results to thresholds
  passed_thresh = test_thresholds(thresh_list[[1]], thresh_list[[2]], thresh_list[[3]], # these 3 are threshold args from null
                                  real_run[[1]],  real_run[[2]], real_run[[3]], # these 3 are posterior probs from real run
                                  method = method)
  # first return object is global selected, second is local selected
  
  
  
  if (method == "global" & nrow(passed_thresh[[1]]) != 0) {
    
    # get rid of duplicated rows if not empty
    selected_global = distinct(passed_thresh[[1]], avg_inclusion_prob, .keep_all = TRUE)
    selected_local = NULL
    
    interact.global = unite(selected_global, "pair", var1, var2, sep = ":")$pair
    interact.local = NULL
    
  } else if (method == "global" & nrow(passed_thresh[[1]]) == 0) {
    
    selected_global = passed_thresh[[1]]
    selected_local = NULL
    
    interact.global = c()
    interact.local = NULL
    
  } else if (method == "local" & nrow(passed_thresh[[1]]) != 0) {
    
    # get rid of duplicated rows if not empty
    selected_local = distinct(passed_thresh[[1]], avg_inclusion_prob, .keep_all = TRUE)
    selected_global = NULL
    
    interact.global = NULL
    interact.local = unite(selected_local, "pair", var1, var2, sep = ":")$pair
    
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
      selected_global = distinct(passed_thresh[[1]], avg_inclusion_prob, .keep_all = TRUE)
      interact.global = unite(selected_global, "pair", var1, var2, sep = ":")$pair
      
    } else {
      
      selected_global = passed_thresh[[1]]
      interact.global = c()
      
    }
    
    if (nrow(passed_thresh[[2]]) != 0) {
      
      # get rid of duplicated rows if not empty
      selected_local = distinct(passed_thresh[[2]], avg_inclusion_prob, .keep_all = TRUE)
      interact.local = unite(selected_local, "pair", var1, var2, sep = ":")$pair
      
    } else {
      
      selected_local = passed_thresh[[2]]
      interact.local = c()
      
    }
    
  }
  
  
  # check which variables selected
  var.select = 
    vip.real %>%
    mutate(select = ifelse(prop.split >= thresh_list[[4]], 1, 0)) %>%
    filter(select == 1)
  
  var.names = vars[var.select$var]
  
  
  # now fit model using selected variables and interactions
  if (method == "global") {
    mod.global = fit.model(data = data, response.type = response.type, response = response,
                           main.effects = var.names, interactions = interact.global,
                           hierarchical = hierarchical,
                           dbarts_x = dbarts_x, dbarts_y = dbarts_y)
    mod.local = NULL
  } else if (method == "local") {
    mod.global = NULL
    mod.local = fit.model(data = data, response.type = response.type, response = response,
                          main.effects = var.names, interactions = interact.local,
                          hierarchical = hierarchical,
                          dbarts_x = dbarts_x, dbarts_y = dbarts_y)
  } else if (method == "both") {
    mod.global = fit.model(data = data, response.type = response.type, response = response,
                           main.effects = var.names, interactions = interact.global,
                           hierarchical = hierarchical,
                           dbarts_x = dbarts_x, dbarts_y = dbarts_y)
    mod.local = fit.model(data = data, response.type = response.type, response = response,
                          main.effects = var.names, interactions = interact.local,
                          hierarchical = hierarchical,
                          dbarts_x = dbarts_x, dbarts_y = dbarts_y)
  }
  
  
  
  # things to return
  to_return = list(selected_global, selected_local, 
                   real_run[[1]], real_run[[2]], real_run[[3]],
                   thresh_list[[1]], thresh_list[[2]], thresh_list[[3]],
                   method, var.select$var, var.names, interact.global, interact.local,
                   mod.global, mod.local)
  
  
  return(to_return)
  
}



