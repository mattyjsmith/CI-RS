################################################################################
# Sensitivity analysis: how much do the relative survival results depend on the level of the expected (life-table) mortality?
# The life-table rates are multiplied by a factor and the relative survival g-formula is repeated, for each dose against placebo.
# 1 = the life table as it is. A "calibrated" factor is also found for each dose: the factor at which the other-cause weights add
# up to the observed number of other-cause deaths (interpolated between the runs). This is an ILLUSTRATION: a factor chosen that way
# uses the cause-of-death information the relative survival method is meant to do without. The excess hazard model cannot be fitted
# when the inflated life-table rates approach the all-cause hazard, so large factors fail: failures are logged in sensitivity_errors.log.
# Run after 01_prepare_data.R (and 03_run_empirical.R, if you want the cause-specific results drawn on the plot).
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("empirical_example"))
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(survival); library(mexhaz); library(boot); library(cmprsk); library(ggplot2)
})
source("02_functions.R")

B_sens         <- 100                         # bootstrap resamples per dose and factor (use 20 for a quick test)
ncpus          <- as.integer(Sys.getenv("N_CORES", "8"))
mx_multipliers <- c(1, 1.5, 2, 2.5, 3, 3.5, 4)
doses  <- c("0.2 mg estrogen", "1.0 mg estrogen", "5.0 mg estrogen")
tag_of <- function(dose) gsub(" ", "", sub(" estrogen$", "", dose))
est_names <- c("RD", "logRR", "RD_comp", "logRR_comp")
RNGkind("L'Ecuyer-CMRG")

one_run <- function(dose, mult, seed, r, r_total) {
  message(sprintf("  [run %d of %d] %s, life-table rates x %s  (started %s)", r, r_total, dose, mult, format(Sys.time(), "%H:%M")))
  dat <- as.data.frame(read_csv(sprintf("analysis_data_%s.csv", tag_of(dose)), show_col_types = FALSE))
  set.seed(seed)
  log_error <- function(what, e) cat(sprintf("%s | x%s | %s | %s\n", dose, mult, what, conditionMessage(e)), file = "sensitivity_errors.log", append = TRUE)
  stat <- function(d) tryCatch(gformula_est(d, "rs", mx_mult = mult), error = function(e) { log_error("resample", e); rep(NA_real_, 4) })
  b <- boot(data = dat, statistic = function(data, indices) stat(data[indices, ]), R = B_sens, strata = dat$rx,
            parallel = if (ncpus > 1) "multicore" else "no", ncpus = ncpus)
  w <- tryCatch(gformula_fit(dat, "rs", mx_mult = mult)$weights, error = function(e) { log_error("full data", e); NULL })
  list(
    est = data.frame(dose = dose, mx_multiplier = mult, estimand = est_names, estimate = b$t0,
                     lower = apply(b$t, 2, quantile, probs = 0.025, na.rm = TRUE),
                     upper = apply(b$t, 2, quantile, probs = 0.975, na.rm = TRUE),
                     failed_resamples = colSums(!is.finite(b$t)), row.names = NULL),
    wts = if (is.null(w)) NULL else {
      d <- w[w$died == 1, ]
      data.frame(dose = dose, mx_multiplier = mult, deaths = nrow(d),
                 observed_prostate = sum(d$eventType_int == 1), sum_wE = sum(d$wE),
                 observed_other = sum(d$eventType_int == 2),    sum_wP = sum(d$wP))
    })
}

grid <- expand.grid(i = seq_along(doses), j = seq_along(mx_multipliers))
grid$id <- seq_len(nrow(grid))

# Reuse the results saved by an earlier run (same data and code): only the dose/factor combinations that are missing, or whose
# point estimate failed, are run. Delete the two sensitivity_lifetable_*.csv files to start from scratch.
est_file <- "sensitivity_lifetable_estimates.csv"; wts_file <- "sensitivity_lifetable_weights.csv"
old_est <- old_wts <- NULL
if (file.exists(est_file) && file.exists(wts_file)) {
  old_est <- read_csv(est_file, show_col_types = FALSE) %>% dplyr::filter(mx_multiplier %in% mx_multipliers)
  old_wts <- read_csv(wts_file, show_col_types = FALSE) %>% dplyr::filter(mx_multiplier %in% mx_multipliers)
  old_est$calibrated <- NULL
  ok_key  <- old_est %>% dplyr::filter(estimand == "RD", !is.na(estimate)) %>%
    dplyr::transmute(key = paste(dose, mx_multiplier)) %>% dplyr::pull(key)
  old_est <- old_est %>% dplyr::filter(paste(dose, mx_multiplier) %in% ok_key)
  old_wts <- old_wts %>% dplyr::filter(paste(dose, mx_multiplier) %in% ok_key)
}
done_key <- if (is.null(old_est)) character(0) else unique(paste(old_est$dose, old_est$mx_multiplier))
todo <- grid[!(paste(doses[grid$i], mx_multipliers[grid$j]) %in% done_key), ]
message(nrow(grid) - nrow(todo), " of ", nrow(grid), " dose/factor combinations already done; running ", nrow(todo))
out  <- lapply(seq_len(nrow(todo)), function(r)
  one_run(doses[todo$i[r]], mx_multipliers[todo$j[r]], seed = 20261102 + todo$id[r], r = r, r_total = nrow(todo)))
est  <- bind_rows(old_est, bind_rows(lapply(out, `[[`, "est")))
wts  <- bind_rows(old_wts, bind_rows(lapply(out, `[[`, "wts")))

# ---- Calibrated factor: where the other-cause weights add up to the observed number of other-cause deaths ----
cal <- wts %>% dplyr::group_by(dose) %>% dplyr::arrange(mx_multiplier, .by_group = TRUE) %>%
  dplyr::summarise(observed_other = dplyr::first(observed_other),
                   theta_cal = if (sum(!is.na(sum_wP)) >= 2) suppressWarnings(approx(sum_wP, mx_multiplier, xout = dplyr::first(observed_other))$y) else NA_real_,
                   .groups = "drop")
cat("\nFactor at which the other-cause weights equal the observed number of other-cause deaths (NA = outside the runs that could be fitted):\n")
print(as.data.frame(cal %>% dplyr::mutate(theta_cal = round(theta_cal, 2))), row.names = FALSE)
cal_runs <- lapply(seq_len(nrow(cal)), function(r) {
  if (is.na(cal$theta_cal[r])) return(NULL)
  one_run(cal$dose[r], round(cal$theta_cal[r], 2), seed = 20261200 + r, r = r, r_total = nrow(cal))
})
cal_runs <- Filter(Negate(is.null), cal_runs)
if (length(cal_runs) > 0) {
  est <- bind_rows(est, bind_rows(lapply(cal_runs, `[[`, "est")))
  wts <- bind_rows(wts, bind_rows(lapply(cal_runs, `[[`, "wts")))
}
est$calibrated <- !(est$mx_multiplier %in% mx_multipliers)

cat("\nRelative survival risk differences at 50 months by the factor applied to the life-table rates:\n")
print(as.data.frame(est %>% dplyr::filter(estimand %in% c("RD", "RD_comp")) %>%
                      dplyr::mutate(across(c(estimate, lower, upper), ~ round(.x, 3)))), row.names = FALSE)
cat("\nWeights (sum over deaths) against the observed causes of death:\n")
print(as.data.frame(wts %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

# ---- Plot: relative survival risk differences against the factor, with the cause-specific results for reference ----
lab <- c(RD = "Prostate cancer death: risk difference", RD_comp = "Other-cause death: risk difference")
res_plot <- est %>% dplyr::filter(estimand %in% names(lab)) %>%
  dplyr::mutate(estimand = unname(lab[estimand]), dose = factor(dose, levels = doses))
p <- ggplot(res_plot, aes(mx_multiplier, estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + geom_line() + geom_point()
if (file.exists("results_empirical.csv")) {
  ref <- read_csv("results_empirical.csv", show_col_types = FALSE) %>%
    dplyr::filter(method %in% c("Cause-specific g-formula", "Aalen-Johansen (non-parametric)"), estimand %in% unname(lab)) %>%
    dplyr::mutate(dose = factor(dose, levels = doses))
  p <- p + geom_hline(data = ref, aes(yintercept = estimate, colour = method), linetype = 2, inherit.aes = FALSE)
}
p <- p + facet_grid(estimand ~ dose, scales = "free_y") +
  labs(x = "Factor applied to the life-table rates", y = "Risk difference at 50 months", colour = "Reference (dashed)",
       title = "Relative survival g-formula against the assumed level of expected mortality") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
print(p)

write_csv(est, "sensitivity_lifetable_estimates.csv")
write_csv(wts, "sensitivity_lifetable_weights.csv")
ggsave("sensitivity_lifetable.pdf", p, width = 12, height = 6)
message("Saved sensitivity_lifetable_estimates.csv, sensitivity_lifetable_weights.csv, sensitivity_lifetable.pdf")
