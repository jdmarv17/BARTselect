# BARTselect

This is a repository for the R package BARTselect. Package documentation is still under development but all code is complete.

Files in this package will allow a user to replicate results in the paper "Detecting Interactions Using Bayesian Additive Regression Trees." Data used for analyses comes from two R packages. The birth weight data comes from MASS::birthwt and the medical expenditures data comes from Ecdat::MedExp. Both packages are readily available for installation from CRAN.

-   R/BARTselect.R contains the code for BARTselect function to detect variables and interactions
-   non_package_files/data_generation.R contains functions to generate data used in simulations
-   non_package_files/real_data/real_data_medexp.R contains analysis code for Ecdat::MedExp data analysis
    -   This file compares BARTselect to Lasso for effect selection and fit models using selected effects
    -   Also use the final selected models to predict on holdout data and compare RMSE to dbarts::bart2 and randomForest::randomForest
-   non_package_files/real_data/real_data_birthwt.R contains code for MASS::birthwt data analysis
    -   This file compares BARTselect to Lasso for effect selection and fit models used selected effects
