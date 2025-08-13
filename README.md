# Causal Inference for the Relative Survival Setting

This repository shows how to estimate causal effects of an exposure on an outcome (in failure-time settings) when the cause of death is unknown.


Estimating causal effects in the relative survival setting involves two main steps:
  1) **Predict** the probability (i.e., create weights) that an all-cause death is due to the disease of interest or other causes. 
  2) **Estimate causal effects** of an exposure on an outcome using an estimator (e.g., g-formula, IPTW, etc.) from within the competing risk framework that can incorporate the weights (from Step 1).

This method can be applied for any disease-specific mortality of interest. Originally, the excess hazard model was developed as an approach to calculate net survival estimates for a sample of patients with cancer. The "excess hazard" is the additional mortality hazard of a disease of interest beyond what is expected in the general population. Thus, these weights can be applied to any disease of interest. However, careful consideration is required in terms of the causal question and whether the life tables are sufficiently stratified regarding the disease of interest. 

This is an adaptation of the work by: 
> Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernán MA. A causal framework for classical statistical estimands in failure time settings with competing events. Statistics in Medicine. 2020 https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8471

---

# Description of necessary files

You do not need to run the simulations to get the results. The results of the simulations are provided in [results.rds](https://github.com/mattyjsmith/CI-RS/blob/main/results.rds). You can simply run the [09_performance_measures]() script, which will load the results of the simulations, to construct figures for the performance measures. The summary of the performance measures is provided in [performance_summary_with_samplesizes.csv](https://github.com/mattyjsmith/CI-RS/blob/main/performance_summary_with_samplesizes.csv).

### Running the simulation study
The [07_run_simulation](https://github.com/mattyjsmith/CI-RS/blob/main/07_run_simulation.R) R script contains the code to run simulations for estimating the total causal effects (i.e., risk difference and risk ratios) using information from the relative survival setting. Ensure you have downloaded the full repository because this script sources the other scripts containing the necessary functions. Be aware that the simulation study can take several weeks on one computer, therefore it is advised that this be run on a high performance computer.

### Expected Mortality Rates
The [LT.xlsx](https://github.com/mattyjsmith/CI-RS/blob/main/LT.xlsx) data set contains the expected mortality rates from life tables stratified by age, sex, deprivation, and calendar year.

### Notable Functions
Notable functions are:

  1) `cDataDesignOptim`[^1] is the function to simulate patient characteristics. 
  2) `cdatasimulationT1WeibOptim`[^1] is the function to simulate the failure time and vital status (only after simulating the patient characteristcs).
  3) `calculateCumInc`[^2] is the function to estimate the cumulative incidence.

---

[^1]: Created by [Aurélien Belot](https://github.com/AurelienBelot). 

[^2]: Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernán MA. A causal framework for classical statistical estimands in failure time settings with competing events. Statistics in Medicine. 2020 https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8471
