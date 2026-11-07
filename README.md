# Code for "A Weighted Estimator for Causal Inference in the Relative Survival Setting"

Simulation study (Sections 3-4 and appendices) and empirical example (Section 5) of the paper by Smith, Kostouraki, Belot, Leyrat,
Maringe, Luque-Fernandez and Rachet. The method weights each observed death by the probability that it was due to the disease
(from an excess hazard model and a population life table) and uses the weights in a g-formula for the total causal effect when the
cause of death is not recorded.

## Folder layout

| Folder | Contents |
|---|---|
| `R/` | Code shared by all simulations: life table, data generation, helper functions, the five analysis methods (bootstrap functions), the parallel runner. |
| `data/` | `LT.xlsx`, the population life table used to generate and analyse the simulated data. |
| `simulations/main/` | Main simulation: truth, bias run, coverage run, performance measures (paper Section 4, main results). |
| `simulations/other_cause/` | Further simulation: treatment effect on other-cause death, younger and older populations. |
| `simulations/cancer_type/` | Further simulation: share of deaths due to the disease (aggressive, medium, indolent). |
| `simulations/time_model/` | Further simulation: linear, interval-specific and spline time in the pooled logistic regressions. |
| `empirical_example/` | Randomised trial of estrogen therapy for prostate cancer (Byar 1980; Section 5). |
| `tables/` | Python script that turns the performance CSVs into the LaTeX table bodies of Section 4. |

Each analysis folder also keeps the small summary files behind the paper's tables (true values, performance CSVs). The large
per-replicate results (`results/`, `merged_sim_results_*.rds`) are not stored here; the run scripts regenerate them.

## Requirements

- R 4.5.1 (the version used for the paper). Packages are installed by `R/01_setup_packages.R` and the `library()` calls at the top of each
  script: `mexhaz`, `Epi`, `plyr`, `dplyr`, `MASS`, `boot`, `readxl`, `splines`, `ggplot2`, `rsimsum`, `readr`, `tidyr`, `knitr`, `here`
  (the empirical example also uses `lubridate`, `stringr`). Python 3 is needed only for the table scripts.
- The simulations run in parallel. Set the number of cores with the environment variable `N_CORES` (default 8), e.g.
  `Sys.setenv(N_CORES = 4)` before running a `07_` or `10_` script.

## How to run

Open `CI-RS.Rproj` in RStudio (or start R in this folder). Every script sets its own working directory with `here::here()`, so no
paths need editing. Restart R before running each script. Scripts are restartable: finished replicates are skipped.

### Main simulation (`simulations/main/`)
1. `08_true_values.R`: the 10,000,000-patient counterfactual truth (`true_values.csv`, included).
2. `07_run_simulation.R` with `run_type <- "bias"`: n = 1000, 3000, 5000; 1000 data sets; five methods. About 4-5 hours on 8 cores.
   Then `07b_merge_results.R` (same `run_type`) and `09a_performance_bias.R`.
3. `07_run_simulation.R` with `run_type <- "coverage"` (n = 1000, 500 data sets, 200 bootstrap resamples; about 27 hours on 8 cores),
   then `07b_merge_results.R` (`run_type <- "coverage"`) and `09_performance_measures.R`.

### Further simulations (run after the main bias run is merged)
- `simulations/other_cause/`: `08_truth_othercause.R` (1,000,000 patients per scenario), `10_run_othercause.R` (about 3 hours on 8 cores),
  `11_performance_othercause.R`. Scenarios are defined once in `00_scenarios.R`.
- `simulations/cancer_type/`: `08_truth_cancertype.R` (10,000,000 patients per scenario, about 35 minutes each), `10_run_cancertype.R`
  (medium and indolent scenarios; about an hour), `11_performance_cancertype.R`. The aggressive scenario reuses the `older_HR1` results
  of `simulations/other_cause/` (run that first).
- `simulations/time_model/`: `10_run_timemodel.R`, `11_performance_timemodel.R`.

### Empirical example
See `empirical_example/README.md`.

### Tables
`python3 tables/make_section4_tables.py` rebuilds the LaTeX table bodies of Section 4 from the performance CSVs.

## Data

- **Life table.** `data/LT.xlsx` holds population mortality rates (`exprate`) by age (`AGERT`), calendar year (`YEARRT`), sex (`SEXRT`) and
  deprivation group (`DEPRT`). It generates the other-cause deaths and provides the expected hazard in the analyses. Source: derived from the Office for National Statistics (ONS) national life tables for the United Kingdom (https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables), stratified by deprivation group.
- **Empirical example.** The Byar trial data are public and downloaded by `00_obtain_data.R`. The US Human Mortality Database rates
  need a free registration and must be downloaded by the user (see `empirical_example/README.md`).

## Reproducibility notes

- Each simulated data set uses `set.seed(j)` for replicate `j`; the truth uses `set.seed(999 + chunk)`. Results should reproduce on the same
  R version; parallel runs on other platforms can differ in the last digits.
- `R/03_data_design.R` and `R/04_simulate_survival.R` are used by the main and time-model simulations. The further simulations use the
  `_extended` versions, which add two arguments (`age_shift`, `BetaOther`) and generate other-cause death times with a shared exponential
  draw in the treated and untreated worlds. With the default arguments the data generated are statistically the same as the main
  simulation, but the random numbers differ, so the main results cannot be reproduced exactly with the extended files.
- `R/05_helpers.R` is adapted from the code of Young et al. (2020), https://github.com/CausalInference/CompetingEvents_Young_SIM_2020.

## Citation

If you use this code, please cite:

Smith MJ, Kostouraki A, Belot A, Leyrat C, Maringe C, Luque-Fernandez MA, Rachet B. A Weighted Estimator for Causal Inference in the
Relative Survival Setting. Manuscript in preparation.

## Licence and acknowledgements

The code is released under the MIT licence (see `LICENSE`). `R/05_helpers.R` is adapted from the code that accompanies Young JG,
Stensrud MJ, Tchetgen Tchetgen EJ, Hernan MA. A causal framework for classical statistical estimands in failure-time settings with
competing events. Statistics in Medicine 2020 (https://github.com/CausalInference/CompetingEvents_Young_SIM_2020, MIT licence).
