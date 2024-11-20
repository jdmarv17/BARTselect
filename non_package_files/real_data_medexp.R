library(Ecdat)
library(MASS)
library(tidyverse)
library(janitor)
library(mpath)
library(randomForest)
library(dbarts)

#source("BARTselect/BARTselect.R")
source("~/non_package_files/BARTselect/BARTselect.R")
data("MedExp")




#################################
##### Ecdat::MedExp Example #####
#################################

### Data pre-processing:
data = MedExp %>% mutate(med = round(med)) # round expenses to the dollar to allow for easier fitting

# set seed and split data
{set.seed(41424); train_ind = sample(1:nrow(data), size = round(0.2 * nrow(data), 0), replace = FALSE)}
train = data[train_ind, ]
remainder = data[-train_ind, ]

{set.seed(41524); predict_ind = sample(1:nrow(remainder), size = round(0.2 * nrow(data), 0), replace = FALSE)}
predict_data = remainder[predict_ind, ]
fit_data = remainder[-predict_ind,]

# list and matrices to record results
record_list = list()
mat_interact = as.data.frame(matrix(nrow = 1, ncol = 2))
mat_var = as.data.frame(matrix(nrow = 1, ncol = 2))
colnames(mat_interact) = c("interact", "sim")
colnames(mat_var) = c("var", "sim")


### Run BART for interaction detection and variable selection 10 times + Run Lasso for same goals:
for (i in 1:20) {
  
  # run BARTselect
  record_list[[i]] = BARTselect(med ~ ., train, num_trees = 10, num_samps = 5000,
                            num_burn = 5000, num_null_run = 20, num_thin = 5,
                            num_chains = 4, num_threads_bart = 4, num_threads_wrangle = 8,
                            prior_power = 4, prior_base = .99, alpha_g = 0.1, alpha_d = 0.1,
                            alpha.g.vip = 0.1, method = "global", set.diff.thresh = TRUE, diff.thresh = 0.275,
                            response.type = "continuous", hierarchical = TRUE)
  
  # record selected interactions
  new_rows_interact = as.data.frame(tibble(interact = record_list[[i]][[12]])) %>% mutate(sim = i)
  mat_interact = bind_rows(mat_interact, new_rows_interact)
  
  # record selected main effects
  new_rows_var = as.data.frame(tibble(var = record_list[[i]][[11]])) %>% mutate(sim = i)
  mat_var = bind_rows(mat_var, new_rows_var)
  
  cat("Sim", i, "\n")
}
mat_interact %>% group_by(interact) %>% summarise(count = n()) %>% filter(count > 1)
mat_var %>% group_by(var) %>% summarise(count = n()) %>% filter(count > 1)

# lasso comparator:
set.seed(51724)
data.lasso = model.matrix(med ~ .^2-1, data = fit_data)
data.lasso.comp = as.data.frame(cbind(med = fit_data[,"med"], data.lasso))
mod_lasso = glmregNB(med ~ ., data = data.lasso.comp,
                     parallel = TRUE, n.cores = 10, nlambda = 1000)
converged = mod_lasso[["converged"]]
which_min = which(mod_lasso[["resdev"]][converged] == min(mod_lasso[["resdev"]][converged]))
min_coefs = mod_lasso[["beta"]][, which_min]
min_nonzero = which(min_coefs != 0)
effects = str_replace_all(str_replace_all(names(min_nonzero), "`", ""), ":", "*")
effects



# Fit interpretable models using selected main effects and interactions:
mod_interpret = glm.nb(med ~ fmde + health + lpi * health + physlim * health, data = fit_data)
mod_lasso_final = glm.nb(med ~ ndisease + linc + child + lc * fmde + lc * lfam + idp * ndisease + 
                               idp * health + lpi * ndisease + lpi * health + fmde * lfam + 
                               physlim * health + physlim * lfam * physlim * sex + ndisease * educdec + 
                               ndisease * black + health * linc + health * educdec + health * sex + 
                               health * black + linc * age + linc * sex + age * child + child * black,
                         data = fit_data, control = list(maxit = 1000, trace = FALSE, epsilon = 1e-8))



### Predictions:

# Lasso selected model
pred_lasso = tibble(preds = predict(mod_lasso_final, predict_data),
                    observed = predict_data$med) %>%
  mutate(sq_error = (preds - observed) ^2)
(rmse_lasso = sqrt(mean(pred_lasso$sq_error)))


# BARTselect selected model
pred_interpret = tibble(preds = predict(mod_interpret, predict_data),
                        observed = predict_data$med) %>%
  mutate(sq_error = (preds - observed) ^2)
(rmse_int = sqrt(mean(pred_interpret$sq_error)))



# BART for regression
mod_bart = dbarts::bart2(med ~ ., data = fit_data, n.trees = 100, n.samples = 5000,
                         n.burn = 5000, n.thin = 5, printCutoffs = 0,
                         verbose = FALSE, keepTrees = TRUE, n.threads = 4, n.chains = 4,
                         power = 6, seed = 51724)

pred_bart = tibble(preds = apply(predict(mod_bart, predict_data, type = "ev"), 2, mean),
                   observed = predict_data$med) %>%
  mutate(sq_error = (preds - observed) ^2)
(rmse_bart = sqrt(mean(pred_bart$sq_error)))

# Random forest for regression
set.seed(51724)
rf_mod = randomForest(med ~ ., data = fit_data, ntree = 100)
pred_rf = tibble(preds = predict(rf_mod, predict_data),
                 observed = predict_data$med) %>%
  mutate(sq_error = (preds - observed) ^2)
(rmse_rf = sqrt(mean(pred_rf$sq_error)))



# get latex code
#stargazer::stargazer(mod_interpret, title = "Results", align = TRUE, no.space = TRUE)




