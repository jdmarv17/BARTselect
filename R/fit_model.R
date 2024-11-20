#' Function to fit user-specified GLM
#'
#' @param data A data frame of original data.
#' @param dbarts_x Matrix created using dbartsData - only covariates.
#' @param dbarts_y Matrix created using dbartsData - only outcome.
#' @param response_type Type of response - must be 'continuous' or 'binary'.
#' @param response Name of the response variable - must be character valued.
#' @param main_effects Vector of predictors to include as main effects - character valued.
#' @param interactions Vector of pairwise interactions to include - must have form 'x_i:x_j' and be character valued.
#' @param hierarchical TRUE/FALSE whether to force both main effects in model if corresponding pairwise interaction is included.
#'
#' @return A list of 2 objects. Object (1) is coefficient tibble, object (2) is the fit model object.
#' 
#' @importFrom stats binomial coef glm lm 
#' @export
#'
#' @examples
#'  library(MASS)
#'  library(dbarts)
#'  library(dplyr)
#'  library(janitor)
#'  data(birthwt)
#'  data = birthwt %>% dplyr::select(., -low)
#'  formula = bwt ~ .
#'  dbarts_data = dbartsData(formula, data)
#'  dbarts_x = clean_names(dbarts_data@x)
#'  dbarts_y = dbarts_data@y
#'  main_effects = c("smoke", "lwt", "race")
#'  interactions = c("race:smoke")
#'  model = fit_model(birthwt, dbarts_x, dbarts_y, response_type = "continuous",
#'                    response = "bwt", main_effects, interactions, hierarchical = TRUE)
#' 
fit_model = function(data, dbarts_x, dbarts_y, response_type = "continuous", response,
                     main_effects = c(), interactions = c(), hierarchical = TRUE) {
  
  # This function fits lm or glm with passed selected main effects and interactions
  # where it assumes interactions are passed in colon form
  ##### Return 1: a tibble of estimated coefficients
  ##### Return 2: lm/glm object using passed main effects and interactions
  # ARGS:
  # `data` = data used to select (same used to fit)
  # `response_type` = binary or continuous response
  # `response` = variable name of response (character valued)
  # `main_effects` = vector of predictors (character valued)
  # `interactions` = vector of interactions (character valued, w/ colon e.g. c("x1:x3", "x5:x9"))
  # `hierarchical` = TRUE/FALSE should both main effects be included if interaction is included 
  
  # do we have interactions and main effects?
  nonzero.int = ifelse(length(interactions) > 0, TRUE, FALSE)
  nonzero.var = ifelse(length(main_effects) > 0, TRUE, FALSE)
  
  # grab outcome to attach to dbarts_data
  outcome = unname(unlist(data[,response]))
  model.data = as.data.frame(cbind(dbarts_x, dbarts_y))
  colnames(model.data) = c(colnames(dbarts_x), response)
  
  
  if (hierarchical) {
    
    if (nonzero.int & nonzero.var) { # detected both main effects and interactions
      # replace ":" with "*" in interactions before formula
      interactions = stringr::str_replace_all(string = interactions, 
                                     pattern = ":",
                                     replacement = "*")
      
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main_effects, collapse = " + "), " + ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (nonzero.int & !nonzero.var) { # detected only interactions
      # replace ":" with "*" in interactions before formula
      interactions = stringr::str_replace_all(string = interactions, 
                                     pattern = ":",
                                     replacement = "*")
      
      formula = formula(paste(as.name(response), " ~ ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & nonzero.var) { # detected only main effects
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main_effects, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & !nonzero.var) { # intercept only model
      formula = formula(paste(as.name(response), " ~ ",
                              1,
                              sep = ""))
    }
    
  } else {
    if (nonzero.int & nonzero.var) { # detected both main effects and interactions
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main_effects, collapse = " + "), " + ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (nonzero.int & !nonzero.var) { # detected only interactions
      formula = formula(paste(as.name(response), " ~ ",
                              paste(interactions, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & nonzero.var) { # detected only main effects
      formula = formula(paste(as.name(response), " ~ ",
                              paste(main_effects, collapse = " + "),
                              sep = ""))
      
    } else if (!nonzero.int & !nonzero.var) { # intercept only model
      formula = formula(paste(as.name(response), " ~ ",
                              1,
                              sep = ""))
    }
  }
  
  
  # depending on response type fit a linear regression or logistic regression
  if (response_type == "continuous") {
    mod = lm(formula = formula, data = model.data)
    coefs = coef(mod)
  } else if (response_type == "binary") {
    mod = glm(formula = formula, data = model.data, family = binomial(link = "logit"))
    coefs = coef(mod)
  } else {
    cat("Expected response_type = 'continuous' or response_type = 'binary'.\n")
    cat("Other response_type values are invalid.\n")
  }
  
  coef.tibble = tibble::rownames_to_column(data.frame(coefs), "var")
  to.return = list(coef.tibble, mod)
  
  return(to.return)
}
