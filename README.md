# Causal Inference for the Relative Survival Setting



This repository shows how to estimate causal effects of an exposure on an outcome (in failure-time settings) when the cause of death is unknown.[^1] 

Estimating causal effects in the relative survival setting involves two main steps:
  1) **Predict** the probability (i.e., create weights) that an all-cause death is due to the disease of interest or other causes. 
  2) **Estimate causal effects** of an exposure on an outcome using an estimator (e.g., g-formula, IPTW, etc.) from within the competing risk framework that can incorporate the weights (from Step 1).

This method can be applied for any disease-specific mortality of interest. Originally, the excess hazard model was developed as an approach to calculate net survival estimates for a sample of patients with cancer. The "excess hazard" is the additional mortality hazard of a disease of interest beyond what is expected in the general population. Thus, these weights can be applied to any disease of interest. However, careful consideration is required in terms of the causal question and whether the life tables are sufficiently stratified regarding the disease of interest. 

This is an adaptation of the work by: 
> Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernán MA. A causal framework for classical statistical estimands in failure time settings with competing events. Statistics in Medicine. 2020 https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8471

---

# Description of necessary files

### Causal Inference and Relative Survival
The [CI&RS](https://github.com/mattyjsmith/CI-RS/blob/main/CI%26RS.R) R script contains the code to estimate causal effects (i.e., risk difference and risk ratios) using information from the relative survival setting.

### Expected Mortality Rates
The [exprates.Rds](https://github.com/mattyjsmith/CI-RS/blob/main/exprates.Rds) data set contains the expected mortality rates from life tables that are stratified by age, sex, year, and another variable (cmb).

### Necessary Functions
The [Functions](https://github.com/mattyjsmith/CI-RS/blob/main/Functions.R) R script contains the necessary functions used within the "CI&RS" script:

  1) `cDataDesignOptim`[^2] is the function to simulate patient characteristics. 
  2) `cdatasimulationT1WeibOptim`[^2] is the function to simulate the failure time and vital status (only after simulating the patient characteristcs).
  3) `calculateCumInc`[^3] is the function to estimate the cumulative incidence.

### Bootstrap function
The [Bootstrap function](https://github.com/mattyjsmith/CI-RS/blob/main/Bootstrap%20function.R) R script contains the bootstrap function used within the "CI&RS" script. The bootstap function calculates confidence intervals for the Risk Difference (RD) by default. If you want confidence intervals for the Relative Risk (RR), you will need to specify this using the ***estimand="RR"*** option.

---

[^1]: This work is a collaboration between members of the [Inequalities in Cancer Outcomes Network](https://icon.lshtm.ac.uk/).

[^2]: Created by [Aurélien Belot](https://github.com/AurelienBelot). 

[^3]: Written by Young *et al* (2020).
