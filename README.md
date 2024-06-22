# BARTselect

This is a repository to store the function files to replicate results in the paper "Detecting Interactions Using Bayesian Additive Regression Trees." Data used for analyses comes from two R packages. The birth weight data comes from MASS::birthwt and the medical expenditures data comes from Ecdat::MedExp. Both packages are readily available for installation from CRAN.

-   BARTselect.R contains all necessary BARTselect functions to run BARTselect with real data examples
-   data_generation.R contains functions to generate data used in simulations
-   read_data_medexp.R contains analysis code for Ecdat::MedExp data analysis
    -   This file compares BARTselect to Lasso for effect selection and fit models using selected effects
    -   Also use the final selected models to predict on holdout data and compare RMSE to dbarts::bart2 and randomForest::randomForest
-   real_data_birthwt.R contains code for MASS::birthwt data analysis
    -   This file comapres BARTselect to Lasso for effect selection and fit models used selected effects
