################################################################################
# Title: Simulation of weighting approach for CI & RS
# Subtitle: Estimate true values from large dataset
# Author: Matthew J. Smith
# Date: 29th July 2025
################################################################################

# Set working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents/LSHTM/ICON/Research Fellow/CI and RS/Playground")
rm(list = ls())

# ---- Load dependencies ----
source("01_setup_packages.R")
source("02_load_lifetable.R")
source("03_data_design.R")
source("04_simulate_survival.R")
source("05_helpers.R")
source("06_bootstrap_function_M1.R")  # Contains the g-formula logic

# ---- Simulation Parameters ----
n_subjects <- 1e7  # Large N to approximate the true value (n=10,000,000)
cutTimes <- 0:23   # Monthly intervals up to 2 years

# ---- Simulate large dataset ----
set.seed(999)
base_data <- cDataDesignOptim(
  n = n_subjects,
  cens.admin = 2015,
  ydiagmin = 2010,
  ydiagmax = 2010,
  betaagecr = -0.5,
  betasex = 0.5,
  betadep = -0.5,
  betastage = 0
)

sim_data_true <- cdatasimulationT1WeibOptim.v2(
  mydat = base_data,
  dist.T1 = "Weibull",
  lambda.weib = 0.25,
  rho.weib = 0.7,
  BetaAge = 0.1,
  Betasex = -0.1,
  Betatrt = 0.5,
  Betadep = 0.4,
  Betastage = 0,
  myexpmort = expratesASD
)

# ---- Get true values by running the bootstrap function once ----
true_values <- boot_function_M1(sim_data_true, indices = 1:nrow(sim_data_true))

# ---- Assign labels to the values for clarity ----
names(true_values) <- c("RD", "logRR", "RD_comp", "logRR_comp")

# ---- Display or save ----
print(true_values)

# Optional: save to file
write.csv(as.data.frame(t(true_values)), "true_values_from_large_data.csv", row.names = FALSE)