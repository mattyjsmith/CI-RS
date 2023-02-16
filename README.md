# Causal Inference for the Relative Survival Setting
This repository shows how to estimate causal effects of an exposure on an outcome (in failure time settings) when the cause of death is unknown. 

Estimating causal effects in the relative survival setting involves two main steps:
  1) Use information on expected mortality rates (obtained from life tables) to estimate the probability (create weights) that an all-cause death is due to cancer or other causes
  2) Choose an estimator (e.g., g-formula, IPTW, etc.) from within the competing risk (or competing event) framework that can incorporate these weights.

This is an adaptation of the work by: 
> Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernán MA. A causal framework for classical statistical estimands in failure time settings with competing events. Statistics in Medicine. 2020 https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8471

This method can be applied for any disease-specific mortality of interest. Originally, the excess hazard model

# Description of necessary files

### CI&RS 
This R script shows the code required to estimate causal effects (i.e., risk difference and risk ratios) using information from the relative survival setting.

### exprates.Rds
This is data for the expected mortality rates that are stratified by age, sex, year, and another variable (cmb).

### Functions.R
This is the script containg the necessary functions used within the "CI&RS" script.

### Bootstrap function.R
This is the script containing the bootstrap function used within the "CI&RS" script. The bootstap function calculates confidence intervals for the Risk Difference (RD) by default. If you want confidence intervals for the Relative Risk (RR), you will need to specify this using the ***estimand="RR"*** option.

