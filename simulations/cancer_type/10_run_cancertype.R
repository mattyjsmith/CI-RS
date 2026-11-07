################################################################################
# Cancer aggressiveness: Method 1 (cause-specific) and Method 2 (relative survival, correct model), spline time model, point
# estimates only, n = 5000, 1000 simulated data sets per scenario, no treatment effect on other-cause death, older population.
# Runs the MEDIUM and INDOLENT scenarios; the aggressive one is the older_HR1 scenario of the OtherCause folder (same data and seeds).
# One file per method and replicate in results/<scenario>/. Restartable; a watchdog stops any worker that hangs.
# A method that fails on a data set is logged in errors.log and saved as NA. For a quick test set sim_ids <- 1:2 after the settings.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "cancer_type"))
rm(list = ls())
RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # R's default generator, as in the other runs: a session that has used
                                                          # RNGkind("L'Ecuyer-CMRG") would otherwise give different data from the same seed
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")   # restart R before sourcing so this takes effect

source(here::here("R", "01_setup_packages.R"))
library(splines)
source("00_scenarios.R")
source(here::here("R", "02_load_lifetable.R"))
source(here::here("R", "03_data_design_extended.R"))
source(here::here("R", "04_simulate_survival_extended.R"))
source(here::here("R", "05_helpers.R"))
source(here::here("R", "06_bootstrap_function_M1.R"))
source(here::here("R", "06_bootstrap_function_M2.R"))
source(here::here("R", "07a_watchdog.R"))

requested_cores <- as.integer(Sys.getenv("N_CORES", "8"))
time_model <- "spline"
cutTimes   <- 0:23
n_subjects <- 5000
sim_ids    <- 1:1000
to_run     <- scenarios[scenarios$scenario != "aggressive", ]

methods <- list(
  "Cause-specific"                    = boot_function_M1,
  "Relative survival – Perfect world" = boot_function_M2
)

run_replicate <- function(j, sc) {
  message(sprintf("%s: simulation %d / %d", sc$scenario, j, max(sim_ids)))
  set.seed(j)
  invisible(capture.output(base_data <- cDataDesignOptim(
    n = n_subjects, cens.admin = 2015, ydiagmin = 2010, ydiagmax = 2010,
    betaagecr = -0.5, betasex = 0.5, betadep = -0.5, betastage = 0, age_shift = sc$age_shift)))
  invisible(capture.output(sim_data <- cdatasimulationT1WeibOptim.v2(
    mydat = base_data, dist.T1 = "Weibull", lambda.weib = sc$lambda, rho.weib = 0.7,
    BetaAge = 0.1, Betasex = -0.1, Betatrt = 0.5, Betadep = 0.4, Betastage = 0,
    BetaOther = 0, myexpmort = expratesASD)))
  sim_data$patno <- 1:n_subjects
  
  for (m in names(methods)) {
    out_file <- file.path("results", sc$scenario, sprintf("sim_%s_n%d_rep%d.rds", gsub(" ", "_", m), n_subjects, j))
    if (file.exists(out_file)) next
    est <- tryCatch(methods[[m]](sim_data, 1:n_subjects), error = function(e) {
      cat(sprintf("%s | rep %d | %s | %s\n", sc$scenario, j, m, conditionMessage(e)), file = "errors.log", append = TRUE)
      rep(NA_real_, 4)
    })
    saveRDS(data.frame(Scenario = sc$scenario, SampleSize = n_subjects, Setting = m, Rep = j,
                       RD = est[1], logRR = est[2], RD_comp = est[3], logRR_comp = est[4]), out_file)
  }
  invisible(NULL)
}

start_time <- Sys.time()
for (s in seq_len(nrow(to_run))) {
  sc <- to_run[s, ]
  dir.create(file.path("results", sc$scenario), recursive = TRUE, showWarnings = FALSE)
  replicate_done <- function(j) all(file.exists(file.path("results", sc$scenario,
    sprintf("sim_%s_n%d_rep%d.rds", gsub(" ", "_", names(methods)), n_subjects, j))))
  message(sprintf("\n===== Scenario %s (lambda %.3f, age shift %d) =====", sc$scenario, sc$lambda, sc$age_shift))
  for (attempt in 1:3) {
    todo <- sim_ids[!vapply(sim_ids, replicate_done, logical(1))]
    if (length(todo) == 0) break
    message(sprintf("%s, pass %d: %d replicates to run", sc$scenario, attempt, length(todo)))
    run_with_watchdog(todo, function(j) run_replicate(j, sc), n_workers = requested_cores, timeout_secs = 1800)
  }
}
message(sprintf("Total time: %.1f hours", as.numeric(difftime(Sys.time(), start_time, units = "hours"))))
