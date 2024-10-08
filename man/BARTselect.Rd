% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BARTselect.R
\name{BARTselect}
\alias{BARTselect}
\title{Function to run automatic BART-based interaction and variable selection}
\usage{
BARTselect(
  formula,
  data,
  num_trees = 10,
  num_samps = 5000,
  num_burn = 5000,
  num_thin = 5,
  num_chains = min(4, dbarts::guessNumCores()),
  alpha_g = 0.01,
  alpha_l = 0.01,
  alpha_d = 0.001,
  alpha_g_vip = 0.01,
  num_null_run = 6,
  method = "global",
  prior_power = 4,
  prior_base = 0.95,
  num_threads_bart = min(dbarts::guessNumCores(), num_chains),
  num_threads_wrangle = dbarts::guessNumCores(),
  set_diff_thresh = FALSE,
  diff_thresh
)
}
\arguments{
\item{formula}{Object of class formula.}

\item{data}{A data frame of original data.}

\item{num_trees}{Number of BART trees per sample.}

\item{num_samps}{Number of posterior samples to take after burn-in.}

\item{num_burn}{Number of burn-in samples.}

\item{num_thin}{Keep every 'num_thin' draw.}

\item{num_chains}{Number of independent MCMC chains to run.}

\item{alpha_g}{Alpha level for global threshold}

\item{alpha_l}{Alpha level for local thresholds}

\item{alpha_d}{Alpha level for global SE thresholds}

\item{alpha_g_vip}{Alpha level for variable selection threshold}

\item{num_null_run}{Number of null runs to do for threshold generations}

\item{method}{A value \%in\% c("global", "local", "both") to specify which selections are returned}

\item{prior_power}{Power parameter for tree prior.}

\item{prior_base}{Base parameter for tree prior.}

\item{num_threads_bart}{Number of threads used to run BART - not recommended to be greater than 'num_chains'.}

\item{num_threads_wrangle}{Number of threads for post processing - can be larger than 'num_chains'.}

\item{set_diff_thresh}{TRUE/FALSE whether user wants to manually set the CMD threshold}

\item{diff_thresh}{A user supplied value for CMD threshold. Only used if set_diff_thresh = TRUE.}
}
\value{
A list of return objects. Object (1) is selected interactions and summary from global threshold. Object (2) is selected interactions
and summary from local thresholds. Objects (3)-(5) are the posterior CIPs, CMDs, and MIPs respectively. Objects (6)-(8)
are the global CIP threshold, local CIP thresholds, and CMD thresholds. Object (9) is a character specifying the threshold type used. Objects (10) and (11)
are selected variable indices and names from variable selection procedure. Objects (12) and (13) are the selected interactions from the global and
local thresholds respectively.
}
\description{
Function to run automatic BART-based interaction and variable selection
}
\examples{
#' library(MASS)
library(dbarts)
library(dplyr)
library(janitor)
data(birthwt)
data = birthwt \%>\% dplyr::select(., -low)
formula = bwt ~ .

select = BARTselect(formula, data, alpha_g = 0.05, alpha_g_vip = 0.1)
 
}
