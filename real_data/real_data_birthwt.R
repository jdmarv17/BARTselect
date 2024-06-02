library(MASS)
library(tidyverse)
library(janitor)
library(glmnet)
library(modelsummary)
library(stargazer)


source("BARTselect/BARTselect.R")
data("birthwt")


#################################
##### MASS::birthwt Example #####
#################################



### Data pre-processing:

# get rid of binarized outcome
birth = 
  birthwt %>%
  dplyr::select(., -low) %>%
  mutate(race = factor(race))


# get number of covariates and set alpha_d
p = ncol(birth) - 1
alpha_d =  (2 / (p * (p - 1)))

# list and matrices to record results
record_list = list()
mat_interact = as.data.frame(matrix(nrow = 1, ncol = 2))
mat_var = as.data.frame(matrix(nrow = 1, ncol = 2))
colnames(mat_interact) = c("interact", "sim")
colnames(mat_var) = c("var", "sim")
set.seed(41324)

### Run BART for interaction detection and variable selection 10 times + Run Lasso for same goals:
# run method 10 times 
for (i in 1:10) {
  # run BARTselect
  record_list[[i]] = BARTid(bwt ~ ., birth, num_trees = 10, num_samps = 10000,
                            num_burn = 5000, num_null_run = 20, num_thin = 5,
                            num_chains = 4, num_threads_bart = 4, num_threads_wrangle = 8,
                            prior_power = 1, prior_base = 0.9, alpha_g = 0.05, alpha_d = alpha_d,
                            alpha.g.vip = 0.1, method = "global", set.diff.thresh = FALSE, diff.thresh = 0,
                            response.type = "continuous", hierarchical = TRUE)
  
  # record selected interactions
  new_rows_interact = as.data.frame(tibble(interact = record_list[[i]][[12]])) %>% mutate(sim = i)
  mat_interact = bind_rows(mat_interact, new_rows_interact)
  
  # record selected main effects
  new_rows_var = as.data.frame(tibble(var = record_list[[i]][[11]])) %>% mutate(sim = i)
  mat_var = bind_rows(mat_var, new_rows_var)
  cat("Sim", i, "\n")
}

# print selected main effects and interactions from BART
mat_interact %>% group_by(interact) %>% summarise(count = n()) %>% filter(count > 1)
mat_var %>% group_by(var) %>% summarise(count = n()) %>% filter(count > 1)

# lasso comparator:
set.seed(41324)
data.lasso = model.matrix(bwt ~ .^2-1, data = birth)
lasso = cv.glmnet(x = data.lasso, y = unlist(birth[,"bwt"]), intercept = TRUE)
coefs = coef(lasso, s = "lambda.min")
nonzero = colnames(data.lasso)[(coefs@i)]
nonzero




### Fit interpretable models using selected effects:

# fit a linear model with selected main effects and interactions from BART
mod_bart = lm(bwt ~ lwt + race * smoke, data = birth)
coef.bart = coef(summary(mod_bart))[, 1:2]

# fit a linear model with selected main effects and interactions from LASSO
mod_lasso = lm(bwt ~ lwt + race + age * race + age * smoke + age * ht + age * ui +
                lwt * ui + race * ptl + race * ht + ptl * ht, data = birth)
coef.lasso = coef(summary(mod_lasso))[, 1:2]

# latex code for table output
stargazer(mod_bart, mod_lasso, title = "Results", align = TRUE, no.space = TRUE, flip = FALSE)



