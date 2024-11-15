# BARTselect

This is a repository for the R package BARTselect. BARTselect allows users to complete interaction and variable selection using the Bayesian Additive Regression Trees (BART) model. Variable selection is completed using an adaptation of that described in *Variable selection for BART: an application to gene regulation* (Bleich et al., 2014)*.* Interaction selection is accomplished using techniques described in *Detecting Interactions Using Bayesian Additive Regression Trees* (Marvald & Love, 2024).

To install the BARTselect package a user should run the following code:

```{r}
library(devtools)
install_github("jdmarv17/BARTselect")

```

While most users will only download the R package itself, analyses from the interaction detection paper can be replicated using files in this repository. Files to do this are contained in 'non_package_files'. Data used for analyses comes from two R packages. The birth weight data comes from MASS::birthwt and the medical expenditures data comes from Ecdat::MedExp. Both packages are readily available for installation from CRAN.

-   R/BARTselect.R contains the code for BARTselect function to detect variables and interactions
-   non_package_files/data_generation.R contains functions to generate data used in simulations
-   non_package_files/real_data/real_data_medexp.R contains analysis code for Ecdat::MedExp data analysis
    -   This file compares BARTselect to Lasso for effect selection and fit models using selected effects
-   non_package_files/real_data/real_data_birthwt.R contains code for MASS::birthwt data analysis
    -   This file compares BARTselect to Lasso for effect selection and fit models used selected effects
