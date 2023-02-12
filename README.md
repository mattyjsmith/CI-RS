# CI-RS
This repository shows how to estimate causal effects of an exposure on an outcome when the cause of death is unknown. 

The essential idea is to use information on expected mortality rates (obtained from life tables) to estimate the probability that an all-cause death is due to cancer or other causes.


This is an adaptation of the work by Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernán MA. A causal framework for classical statistical estimands in failure time settings with competing events. Statistics in Medicine. 2020 https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8471


# Description of necessary files

### CI&RS 
This R script shows the code required to estimate causal effects (i.e., risk difference and risk ratios) using information from the relative survival setting.

### exprates.Rds
This is data for the expected mortality rates (stratified by age, sex, year, and comorbidity).

### Functions.R
This is the script containg the necessary functions used within the "CI&RS" script.



