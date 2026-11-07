################################################################################
# Title: Simulation of weighting approach for CI & RS
# Author: Matthew J. Smith
# Date: 29th July 2025
# CORRECTED VERSION OF THE PAPER'S SIMULATION (see README_paper.md): same scenarios and methods as the original,
# with the corrections found in October 2026 applied.
################################################################################

# Set working directory
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "main"))
rm(list = ls())
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1")   # restart R before sourcing so this takes effect

# ---- Load dependencies ----
{
source(here::here("R", "01_setup_packages.R"))        # Install and load packages
library(splines)                     # natural splines for the time model
source(here::here("R", "02_load_lifetable.R"))        # Load the life table
source(here::here("R", "03_data_design.R"))           # Simulate covariates and treatment assignment
source(here::here("R", "04_simulate_survival.R"))     # Simulate survival time
source(here::here("R", "05_helpers.R"))               # Functions to calculate cumulative incidence
source(here::here("R", "06_bootstrap_function_M1.R")) # Bootstrap function for cause-specific
source(here::here("R", "06_bootstrap_function_M2.R")) # Bootstrap function for RS – perfect world
source(here::here("R", "06_bootstrap_function_M3.R")) # Bootstrap function for RS – No Dep in EHM
source(here::here("R", "06_bootstrap_function_M4.R")) # Bootstrap function for RS – No Dep in LT
source(here::here("R", "06_bootstrap_function_M5.R")) # Bootstrap function for RS – No Dep in EHM or LT
source(here::here("R", "07a_watchdog.R"))             # Runs replicates in parallel; stops any worker that hangs
}

# ---- Core allocation ----
requested_cores <- as.integer(Sys.getenv("N_CORES", "8"))  

# ---- Settings ----
# ---- Run type ----
# "bias"     : MAIN RESULTS. Point estimates only (no bootstrap); n = 1000, 3000, 5000; 1000 simulated data sets.
#              Gives bias, relative bias and empirical SE. Much faster, and the replicates run in parallel.
# "coverage" : APPENDIX. Bootstrap intervals (boot_reps resamples); n = 1000; 500 simulated data sets.
#              Gives coverage, bias-eliminated coverage and relative error.
# The two run types write to separate folders (results/bias, results/coverage).
run_type <- "coverage"
if (run_type == "bias")     { sample_sizes <- c(1000, 3000, 5000); sim_ids <- 1:1000; do_bootstrap <- FALSE }
if (run_type == "coverage") { sample_sizes <- c(1000);             sim_ids <- 1:500;  do_bootstrap <- TRUE  }
# For a quick test, set sim_ids <- 1:2 (and e.g. sample_sizes <- 1000) on a line after these.
results_dir <- file.path("results", run_type)

# ---- Time model of the pooled logistic regressions ----
# "spline": natural cubic spline in time (4 df);  "linear": linear time (the earlier Paper folder);  "factor": interval-specific intercepts
time_model <- "spline"
#n_sims     <- 2 # Required at least 650
cutTimes   <- 0:23 # We truncate at 2-years-since-diagnosis, then split into months (i.e., 24 intervals)
boot_reps  <- 200 # Number of bootstrap replicates. This simulation is extremely computationally heavy!
bootstrap_to_run <- c(
  "Cause-specific",
  "Relative survival – Perfect world",
  "Relative survival – Imperfect world (without Dep in EHM)",
  "Relative survival – Imperfect world (without Dep in LT)",
  "Relative survival – Imperfect world (without Dep in EHM or LT)"  
)

# ---- Bootstrap Configurations ----
{
bootstrap_configs <- list(
  "Cause-specific" = list(
    fn = boot_function_M1,
    pre_boot_summary = function(sim_data) {

      sim_data$cause   <- ifelse(sim_data$Timey >= 2, 0, sim_data$cause)
      sim_data$vstatus <- ifelse(sim_data$Timey >= 2, 0, sim_data$vstatus)
      sim_data$Timey   <- ifelse(sim_data$Timey >= 2, 2, sim_data$Timey)
      
      list(
        SumW = sum(sim_data$cause == 1),  
        NY   = sum(sim_data$cause == 1),
        SumP = sum(sim_data$cause == 2),
        ND   = sum(sim_data$cause == 2)
      )
    }
  ),
  "Relative survival – Perfect world" = list(
    fn = boot_function_M2,
    pre_boot_summary = function(sim_data) {
      
      invisible(capture.output(ehmodel <- suppressMessages(
        mexhaz(Surv(Timey, vstatus) ~ agecr + IsexH + dep + trt, 
               data = sim_data, base = "weibull", expected = "exprate")
        )
      ))
      
      sim_data$excesshazard <- vapply(
        seq_len(nrow(sim_data)),
        function(k) {
          predict(ehmodel,
                  time.pts = sim_data$Timey[k],
                  data.val = sim_data[k, , drop = FALSE]
          )$results$hazard
        },
        numeric(1)
      )
      
      sim_data$weightsE <- sim_data$excesshazard / (sim_data$excesshazard + sim_data$exprate)
      sim_data$weightsP <- sim_data$exprate / (sim_data$excesshazard + sim_data$exprate)
      
      sim_data$cause <- ifelse(sim_data$Timey >= 2, 0, sim_data$cause)
      sim_data$vstatus <- ifelse(sim_data$Timey >= 2, 0, sim_data$vstatus)
      sim_data$Timey <- ifelse(sim_data$Timey >= 2, 2, sim_data$Timey)
      
      list(
        SumW = sum(sim_data$weightsE[sim_data$vstatus == 1]),
        NY   = sum(sim_data$cause == 1),
        SumP = sum(sim_data$weightsP[sim_data$vstatus == 1]),
        ND   = sum(sim_data$cause == 2)
      )
    }
  ),
  "Relative survival – Imperfect world (without Dep in EHM)" = list(
    fn = boot_function_M3,
    pre_boot_summary = function(sim_data) {
      
      invisible(capture.output(ehmodel <- suppressMessages(
        mexhaz(Surv(Timey, vstatus) ~ agecr + IsexH + trt, 
               data = sim_data, base = "weibull", expected = "exprate")
        )
      ))
      
      sim_data$excesshazard <- vapply(
        seq_len(nrow(sim_data)),
        function(k) {
          predict(ehmodel,
                  time.pts = sim_data$Timey[k],
                  data.val = sim_data[k, , drop = FALSE]
          )$results$hazard
        },
        numeric(1)
      )
      
      sim_data$weightsE <- sim_data$excesshazard / (sim_data$excesshazard + sim_data$exprate)
      sim_data$weightsP <- sim_data$exprate / (sim_data$excesshazard + sim_data$exprate)
      
      sim_data$cause <- ifelse(sim_data$Timey >= 2, 0, sim_data$cause)
      sim_data$vstatus <- ifelse(sim_data$Timey >= 2, 0, sim_data$vstatus)
      sim_data$Timey <- ifelse(sim_data$Timey >= 2, 2, sim_data$Timey)
      
      list(
        SumW = sum(sim_data$weightsE[sim_data$vstatus == 1]),
        NY   = sum(sim_data$cause == 1),
        SumP = sum(sim_data$weightsP[sim_data$vstatus == 1]),
        ND   = sum(sim_data$cause == 2)
      )
    }
  ),
  "Relative survival – Imperfect world (without Dep in LT)" = list(
    fn = boot_function_M4,
    pre_boot_summary = function(sim_data) {
      
      invisible(capture.output(ehmodel <- suppressMessages(
        mexhaz(Surv(Timey, vstatus) ~ agecr + IsexH + dep + trt, 
               data = sim_data, base = "weibull", expected = "expratesMean")
      )))
      
      sim_data$excesshazard <- vapply(
        seq_len(nrow(sim_data)),
        function(k) {
          predict(ehmodel,
                  time.pts = sim_data$Timey[k],
                  data.val = sim_data[k, , drop = FALSE]
          )$results$hazard
        },
        numeric(1)
      )
      
      sim_data$weightsE <- sim_data$excesshazard / (sim_data$excesshazard + sim_data$expratesMean)
      sim_data$weightsP <- sim_data$expratesMean / (sim_data$excesshazard + sim_data$expratesMean)
      
      sim_data$cause   <- ifelse(sim_data$Timey >= 2, 0, sim_data$cause)
      sim_data$vstatus <- ifelse(sim_data$Timey >= 2, 0, sim_data$vstatus)
      sim_data$Timey   <- ifelse(sim_data$Timey >= 2, 2, sim_data$Timey)
      
      list(
        SumW = sum(sim_data$weightsE[sim_data$vstatus == 1]),
        NY   = sum(sim_data$cause == 1),
        SumP = sum(sim_data$weightsP[sim_data$vstatus == 1]),
        ND   = sum(sim_data$cause == 2)
      )
    }
  ),
  "Relative survival – Imperfect world (without Dep in EHM or LT)" = list(
    fn = boot_function_M5,
    pre_boot_summary = function(sim_data) {
      
      invisible(capture.output(ehmodel <- suppressMessages(
        mexhaz(Surv(Timey, vstatus) ~ agecr + IsexH + trt,
               data = sim_data, base = "weibull", expected = "expratesMean")
      )))
      
      sim_data$excesshazard <- vapply(
        seq_len(nrow(sim_data)),
        function(k) {
          predict(ehmodel,
                  time.pts = sim_data$Timey[k],
                  data.val = sim_data[k, , drop = FALSE]
          )$results$hazard
        },
        numeric(1)
      )
      
      sim_data$weightsE <- sim_data$excesshazard / (sim_data$excesshazard + sim_data$expratesMean)
      sim_data$weightsP <- sim_data$expratesMean / (sim_data$excesshazard + sim_data$expratesMean)
      
      sim_data$cause   <- ifelse(sim_data$Timey >= 2, 0, sim_data$cause)
      sim_data$vstatus <- ifelse(sim_data$Timey >= 2, 0, sim_data$vstatus)
      sim_data$Timey   <- ifelse(sim_data$Timey >= 2, 2, sim_data$Timey)
      
      list(
        SumW = sum(sim_data$weightsE[sim_data$vstatus == 1]),
        NY   = sum(sim_data$cause == 1),
        SumP = sum(sim_data$weightsP[sim_data$vstatus == 1]),
        ND   = sum(sim_data$cause == 2)
      )
    }
  )
)
}

# ---- Simulation Function ----
{
run_simulation <- function(boot_function, summary_fn, setting_label, sim_id, sim_data,
                           parallel_type = "multicore", n_cores = 6) {
  if (!do_bootstrap) {
    # Point estimates only (no bootstrap): used for the bias results
    t0 <- boot_function(sim_data, indices = 1:nrow(sim_data))
    return(data.frame(
      SampleSize = n_subjects, Setting = setting_label, Rep = sim_id,
      RD = t0[1], SE_RD = NA, CI_RD_L = NA, CI_RD_U = NA, BS_RD = NA,
      logRR = t0[2], SE_logRR = NA, CI_logRR_L = NA, CI_logRR_U = NA, BS_logRR = NA,
      RD_comp = t0[3], SE_RDcomp = NA, CI_RDcomp_L = NA, CI_RDcomp_U = NA, BS_RDcomp = NA,
      logRR_comp = t0[4], SE_logRRcomp = NA, CI_logRRcomp_L = NA, CI_logRRcomp_U = NA, BS_logRRcomp = NA,
      SumW = NA, NY = NA, SumP = NA, ND = NA
    ))
  }
  
  summary_values <- summary_fn(sim_data)
  
  boot_progress_env$counter <- 0
  
  # Handle multicore vs snow
  if (parallel_type == "snow") {
    cl <- parallel::makeCluster(n_cores)
    bs_result <- boot(sim_data, statistic = boot_function, R = boot_reps,
                      parallel = "snow", ncpus = n_cores, cl = cl)
    stopCluster(cl)
  } else {
    bs_result <- boot(sim_data, statistic = boot_function, R = boot_reps,
                      parallel = "multicore", ncpus = n_cores)
  }
  
  # Extract confidence intervals
  ci1 <- boot.ci(bs_result, type = "norm", index = 1) # CI for RD (risk of disease-related death)
  ci2 <- boot.ci(bs_result, type = "norm", index = 2) # CI for RR (risk of disease-related death)
  ci3 <- boot.ci(bs_result, type = "norm", index = 3) # CI for RD (risk of other-cause death)
  ci4 <- boot.ci(bs_result, type = "norm", index = 4) # CI for RR (risk of other-cause death)
  
  data.frame(
    SampleSize = n_subjects,
    Setting = setting_label,
    Rep = sim_id,
    RD = bs_result$t0[1],
    SE_RD = sd(bs_result$t[, 1]),
    CI_RD_L = ci1$normal[2],  # CI for RD (risk of disease-related death) - Lower
    CI_RD_U = ci1$normal[3],  # CI for RD (risk of disease-related death) - Upper
    BS_RD = mean(bs_result$t[, 1]), # RD (risk of disease-related death)
    logRR = bs_result$t0[2],
    SE_logRR = sd(bs_result$t[, 2]),
    CI_logRR_L = ci2$normal[2],
    CI_logRR_U = ci2$normal[3],
    BS_logRR = mean(bs_result$t[, 2]), # RR (risk of disease-related death)
    RD_comp = bs_result$t0[3],
    SE_RDcomp = sd(bs_result$t[, 3]),
    CI_RDcomp_L = ci3$normal[2],
    CI_RDcomp_U = ci3$normal[3],
    BS_RDcomp = mean(bs_result$t[, 3]), # RD (risk of other-cause death)
    logRR_comp = bs_result$t0[4],
    SE_logRRcomp = sd(bs_result$t[, 4]),
    CI_logRRcomp_L = ci4$normal[2],
    CI_logRRcomp_U = ci4$normal[3],
    BS_logRRcomp = mean(bs_result$t[, 4]), # RR (risk of other-cause death)
    SumW = summary_values$SumW,
    NY   = summary_values$NY,
    SumP = summary_values$SumP,
    ND   = summary_values$ND
  )
}
}

# ---- Main Simulation Loop ----
{
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  boot_progress_env <- new.env()
  boot_progress_env$counter <- 0

  start_time <- Sys.time()
  
  # One simulated data set (seed = j): all five methods, one results file per method
  run_replicate <- function(j, n_cores_boot) {
    message(sprintf("Simulation %d / %d (n = %d)", j, max(sim_ids), n_subjects))
    set.seed(j)  # Set the seed once for this replicate
    
    # Simulate the data
    invisible(capture.output(base_data <- cDataDesignOptim(
      n = n_subjects,
      cens.admin = 2015,
      ydiagmin = 2010,
      ydiagmax = 2010,
      betaagecr = -0.5,
      betasex = 0.5,
      betadep = -0.5,
      betastage = 0
    )))
    
    invisible(capture.output(sim_data <- cdatasimulationT1WeibOptim.v2(
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
    )))
    
    sim_data$patno <- 1:n_subjects
    
    for (setting in names(bootstrap_configs)) {
      if (!(setting %in% bootstrap_to_run)) next
      cfg <- bootstrap_configs[[setting]]
      
      # Skip replicates already saved (lets the run be stopped and restarted)
      out_file <- file.path(results_dir, sprintf("sim_%s_n%d_rep%d.rds", gsub(" ", "_", setting), n_subjects, j))
      if (file.exists(out_file)) next
      
      res <- run_simulation(
        boot_function = cfg$fn,
        summary_fn = cfg$pre_boot_summary,
        setting_label = setting,
        sim_id = j,
        sim_data = sim_data,
        parallel_type = "multicore",
        n_cores = n_cores_boot
      )
      saveRDS(res, file = out_file)
      rm(res); gc()
    }
    invisible(NULL)
  }
  
  # TRUE when the files of all five methods exist for replicate j (at the current n_subjects)
  replicate_done <- function(j) {
    all(file.exists(file.path(results_dir,
        sprintf("sim_%s_n%d_rep%d.rds", gsub(" ", "_", bootstrap_to_run), n_subjects, j))))
  }
  
  for (n_subjects in sample_sizes) {
    message(sprintf("Running simulations for n = %d (run type: %s)", n_subjects, run_type))
    if (do_bootstrap) {
      # Replicates one after another; the bootstrap resamples use all the cores
      for (j in sim_ids) run_replicate(j, requested_cores)
    } else {
      # No bootstrap: the replicates themselves run in parallel. A watchdog stops any worker that hangs, and
      # replicates that did not finish are run again in the next pass (finished ones are skipped without being re-simulated).
      for (attempt in 1:3) {
        todo <- sim_ids[!vapply(sim_ids, replicate_done, logical(1))]
        if (length(todo) == 0) break
        message(sprintf("n = %d, pass %d: %d replicates to run", n_subjects, attempt, length(todo)))
        run_with_watchdog(todo, function(j) run_replicate(j, 1), n_workers = requested_cores, timeout_secs = 1800)
      }
    }
  }

#results <- do.call(rbind, all_results)
#results <- results[order(results$Setting, results$Rep), ]

end_time <- Sys.time()


# ---- Print Total Time in Clear Format ----
total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))

days    <- total_secs %/% (24 * 3600)
hours   <- (total_secs %% (24 * 3600)) %/% 3600
minutes <- (total_secs %% 3600) %/% 60
seconds <- round(total_secs %% 60)

time_str <- sprintf(
  "Total time: %dd %02dh %02dm %02ds",
  days, hours, minutes, seconds
)
message(time_str)
}



# ---- Combine results into one data set ----
# 1. List all files in the folder
rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)

# 2. Optionally filter files by pattern (e.g. only n=1000)
#rds_files <- rds_files[grepl("n1000", rds_files)]

# 3. Read all .rds files into a list
all_results_list <- lapply(rds_files, readRDS)

# 4. Combine into a single data frame
combined_results <- do.call(rbind, all_results_list)

# Save results
# save(results, file = "results/sim_results.RData")
# saveRDS(results, file = "results/sim_results.rds")


