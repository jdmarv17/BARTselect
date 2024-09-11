## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installation-------------------------------------------------------------
#devtools::install_github("jdmarv17/BARTselect")
library(BARTselect)

## ----libraries, message=FALSE, warning=FALSE----------------------------------
library(dbarts)
library(MASS)
library(dplyr)

## ----data---------------------------------------------------------------------
birth_weights = 
  MASS::birthwt %>%
  dplyr::mutate(race = factor(race, levels = c(1, 2, 3))) %>%
  dplyr::select(., -low)

head(birth_weights)

## ----BARTselect---------------------------------------------------------------
p = ncol(birth_weights) - 1
alpha_d =  (8 / (p * (p - 1)))
chains = 4
bart_threads = 4
wrangle_threads = 8

selections = BARTselect(bwt ~ ., data = birth_weights,
                        num_trees = 10, num_samps = 10000, num_burn = 5000,
                        num_chains = chains, num_thin = 5, num_null_run = 10,
                        num_threads_bart = bart_threads, num_threads_wrangle = wrangle_threads,
                        alpha_g = 0.1, alpha_g_vip = 0.1,
                        alpha_d = alpha_d, set_diff_thresh = FALSE,
                        prior_power = 1, prior_base = 0.95, method = "global",
                        response_type = "continuous", hierarchical = TRUE)

## -----------------------------------------------------------------------------
# selected variables:
selections[[11]]
# selected interactions:
selections[[12]]

## -----------------------------------------------------------------------------
mod = lm(bwt ~ lwt + race*smoke, data = birth_weights)
summary(mod)

