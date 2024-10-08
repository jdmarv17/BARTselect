% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_posteriors.R
\name{get_posteriors}
\alias{get_posteriors}
\title{Function to get posterior MIPs, CIPs, and CMDs}
\usage{
get_posteriors(indicators, vars, num_threads_wrangle)
}
\arguments{
\item{indicators}{A data frame of indicators. This data frame is output by run_bart and has one column
for each BART covariate and one row for each posterior tree.}

\item{vars}{Column names resulting from dbartsData(data).}

\item{num_threads_wrangle}{Number of threads for processing - can be larger than 'num_chains'}
}
\value{
A list with 3 objects. First is the matrix of CIPs, second is matrix of CMDs, third is vector of MIPs
}
\description{
Function to get posterior MIPs, CIPs, and CMDs
}
\examples{
library(MASS)
library(dbarts)
library(dplyr)
library(janitor)
data(birthwt)
data = birthwt \%>\% dplyr::select(., -low)
formula = bwt ~ .
dbarts_data = dbartsData(formula, data)
dbarts_x = clean_names(dbarts_data@x)
dbarts_y = dbarts_data@y
vars = colnames(dbarts_x)
bart_run = run_bart(formula, data, vars = vars)

posteriors = get_posteriors(bart_run[[1]], vars, parallel::detectCores())

}
\keyword{internal}
