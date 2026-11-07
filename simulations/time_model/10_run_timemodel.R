################################################################################
# TIME-MODEL CHECK. Methods 1 (cause-specific) and 2 (relative survival, correct model) at n = 5000, with the time model of
# the pooled logistic regressions varied:  interval-specific intercepts ("factor")  and a natural cubic spline ("spline").
# The linear time term is the main simulation (Paper folder): same data, same seeds, so the comparison is paired.
# Point estimates only (no bootstrap). The slow excess hazard step is done once per data set and shared by the variants.
# Restartable; a watchdog stops any worker that hangs. For a quick test use sim_ids <- 1:2 on a line after the settings.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "time_model"))
rm(list = ls())
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")   # restart R before sourcing so this takes effect

source(here::here("R", "01_setup_packages.R"))
library(splines)
source(here::here("R", "02_load_lifetable.R"))
source(here::here("R", "03_data_design.R"))
source(here::here("R", "04_simulate_survival.R"))
source(here::here("R", "05_helpers.R"))
source(here::here("R", "06_M1_timemodel.R"))
source(here::here("R", "06_M2_timemodel.R"))
source(here::here("R", "07a_watchdog.R"))

requested_cores <- as.integer(Sys.getenv("N_CORES", "8"))
cutTimes    <- 0:23
n_subjects  <- 5000
sim_ids     <- 1:1000
time_models <- c("linear", "factor", "spline")

# Excess hazard of the correct model (same code as inside Method 2), done once per data set
excess_hazard_perfect <- function(sim_data) {
  invisible(capture.output(ehmodel <- suppressMessages(
    mexhaz(Surv(Timey, vstatus) ~ agecr + IsexH + dep + trt, data = sim_data, base = "weibull", expected = "exprate"))))
  vapply(seq_len(nrow(sim_data)), function(k) {
    predict(ehmodel, time.pts = sim_data$Timey[k], data.val = sim_data[k, , drop = FALSE])$results$hazard
  }, numeric(1))
}

out_file <- function(tm, method_slug, j) file.path("results", tm, sprintf("sim_%s_n%d_rep%d.rds", method_slug, n_subjects, j))

run_replicate <- function(j) {
  message(sprintf("simulation %d / %d", j, max(sim_ids)))
  set.seed(j)
  invisible(capture.output(base_data <- cDataDesignOptim(
    n = n_subjects, cens.admin = 2015, ydiagmin = 2010, ydiagmax = 2010,
    betaagecr = -0.5, betasex = 0.5, betadep = -0.5, betastage = 0)))
  invisible(capture.output(sim_data <- cdatasimulationT1WeibOptim.v2(
    mydat = base_data, dist.T1 = "Weibull", lambda.weib = 0.25, rho.weib = 0.7,
    BetaAge = 0.1, Betasex = -0.1, Betatrt = 0.5, Betadep = 0.4, Betastage = 0,
    myexpmort = expratesASD)))
  sim_data$patno <- 1:n_subjects
  sim_data$excesshazard <- excess_hazard_perfect(sim_data)
  rows <- 1:n_subjects
  
  for (tm in time_models) {
    est1 <- boot_function_M1_tm(sim_data, rows, time_model = tm)
    est2 <- boot_function_M2_tm(sim_data, rows, time_model = tm)
    for (m in list(list("Cause-specific", "Cause-specific", est1),
                   list("Relative survival – Perfect world", "Relative_survival_–_Perfect_world", est2))) {
      saveRDS(data.frame(TimeModel = tm, SampleSize = n_subjects, Setting = m[[1]], Rep = j,
                         RD = m[[3]][1], logRR = m[[3]][2], RD_comp = m[[3]][3], logRR_comp = m[[3]][4]),
              out_file(tm, m[[2]], j))
    }
  }
  invisible(NULL)
}

replicate_done <- function(j) all(file.exists(unlist(lapply(time_models, function(tm)
  c(out_file(tm, "Cause-specific", j), out_file(tm, "Relative_survival_–_Perfect_world", j))))))

for (tm in time_models) dir.create(file.path("results", tm), recursive = TRUE, showWarnings = FALSE)
start_time <- Sys.time()
for (attempt in 1:3) {
  todo <- sim_ids[!vapply(sim_ids, replicate_done, logical(1))]
  if (length(todo) == 0) break
  message(sprintf("pass %d: %d replicates to run", attempt, length(todo)))
  run_with_watchdog(todo, run_replicate, n_workers = requested_cores, timeout_secs = 1800)
}
message(sprintf("Total time: %.1f hours", as.numeric(difftime(Sys.time(), start_time, units = "hours"))))
