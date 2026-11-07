# Empirical example: estrogen therapy and prostate cancer mortality

The Byar (1980) randomised trial (stage III/IV prostate cancer), as analysed by Young et al. (2020). Each estrogen dose (0.2 mg, 1.0 mg,
5.0 mg) is compared with the same placebo group; follow-up is truncated at 50 months and cancer death is the outcome, with other-cause
death as the competing event. The cause-specific g-formula (cause of death known) is compared with the relative survival g-formula
(excess hazard weights). Expected mortality comes from the US Human Mortality Database (mx).

## Data you need to provide

1. Register (free) at https://www.mortality.org, download "United States, Life tables (period 1x1), Total", and save the file as
   `empirical_example/data/USHMD.txt`. The Human Mortality Database does not allow redistribution of its data, so it is not in this repo.
2. The Byar trial data are downloaded from https://github.com/IowaBiostat/data-sets by `00_obtain_data.R`.

## Order of work (run from the repo root; restart R before each script)

1. `00_obtain_data.R`: downloads the trial data and joins the life-table rates; writes `data/data_50months.csv`.
2. `01_prepare_data.R`: analysis data sets `analysis_data_0.2mg.csv`, `analysis_data_1.0mg.csv`, `analysis_data_5.0mg.csv`.
3. `04_diagnostics.R` (seconds): observed versus life-table-expected other-cause deaths by arm and subgroup.
4. `test_empirical.R` (a few minutes): the relative survival method given the true cause of death must equal the cause-specific method.
5. `03_run_empirical.R`: bootstrap for the three doses. Set `B` at the top of the script (20 for a quick test). Writes
   `results_empirical.csv`, `risks_50m.csv`, `weights_diagnostic.csv`, `cif_empirical.pdf`.
6. `05_sensitivity_lifetable.R`: relative survival with the life-table rates multiplied by a range of factors, and a calibrated factor per
   dose. Set `B_sens` at the top. Writes `sensitivity_lifetable_estimates.csv`, `sensitivity_lifetable_weights.csv`,
   `sensitivity_lifetable.pdf`. Run after step 5 (it draws the cause-specific results on the plot).
7. `python3 make_latex_tables.py`: LaTeX table bodies from the CSV outputs.

`02_functions.R` holds the functions used by steps 5 and 6. The CSV and PDF outputs from the run reported in the paper are included.
