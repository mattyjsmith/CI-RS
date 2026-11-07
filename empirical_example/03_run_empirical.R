################################################################################
# Empirical example (Byar 1980 trial): for each estrogen dose against placebo, the non-parametric (Aalen-Johansen),
# cause-specific g-formula (Young et al.) and relative survival g-formula, with 95% percentile bootstrap intervals
# (resampling within each arm), 50 months of follow-up.
# Run 01_prepare_data.R first. For a quick test set B <- 20.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("empirical_example"))
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(survival); library(mexhaz); library(boot); library(cmprsk); library(ggplot2)
})
source("02_functions.R")

B      <- 200                                 # bootstrap resamples (use 20 for a quick test)
batch  <- 10                                  # resamples per batch: a progress line is printed after each batch
ncpus  <- as.integer(Sys.getenv("N_CORES", "8"))
doses  <- c("0.2 mg estrogen", "1.0 mg estrogen", "5.0 mg estrogen")
tag_of <- function(dose) gsub(" ", "", sub(" estrogen$", "", dose))
RNGkind("L'Ecuyer-CMRG")

est_names <- c("RD", "logRR", "RD_comp", "logRR_comp")
methods <- list(
  "Aalen-Johansen (non-parametric)" = aj_est,
  "Cause-specific g-formula"        = function(d) gformula_est(d, "young"),
  "Relative survival g-formula"     = function(d) gformula_est(d, "rs")
)

analyse_dose <- function(dose, seed, i) {
  message("\n===== Dose ", i, " of ", length(doses), ": ", dose, " vs placebo  (started ", format(Sys.time(), "%H:%M"), ") =====")
  dat <- as.data.frame(read_csv(sprintf("analysis_data_%s.csv", tag_of(dose)), show_col_types = FALSE))
  set.seed(seed)

  # Bootstrap in batches, printing progress (resamples done, minutes elapsed, estimated minutes left) after each batch
  run_boot <- function(stat, label) {
    # A failed resample is recorded in bootstrap_errors.log (dose | method | message) and counted in failed_resamples
    safe <- function(d) tryCatch(stat(d), error = function(e) {
      cat(sprintf("%s | %s | %s\n", dose, label, conditionMessage(e)), file = "bootstrap_errors.log", append = TRUE)
      rep(NA_real_, 4)
    })
    sizes <- rep(batch, B %/% batch); if (B %% batch > 0) sizes <- c(sizes, B %% batch)
    t_all <- NULL; t0 <- NULL; start <- Sys.time()
    for (sz in sizes) {
      b <- boot(data = dat, statistic = function(data, indices) safe(data[indices, ]), R = sz, strata = dat$rx,
                parallel = if (ncpus > 1) "multicore" else "no", ncpus = ncpus)
      t_all <- rbind(t_all, b$t); if (is.null(t0)) t0 <- b$t0
      done <- nrow(t_all); el <- as.numeric(difftime(Sys.time(), start, units = "mins"))
      message(sprintf("      %s: %d / %d resamples, %.1f min elapsed, about %.1f min left",
                      label, done, B, el, el * (B - done) / done))
    }
    list(t0 = t0, t = t_all)
  }

  # ---- Estimates and bootstrap intervals ----
  res <- list()
  for (m in names(methods)) {
    message("    ", m)
    b <- run_boot(methods[[m]], m)
    res[[m]] <- data.frame(dose = dose, method = m, estimand = est_names, estimate = b$t0,
                           lower = apply(b$t, 2, quantile, probs = 0.025, na.rm = TRUE),
                           upper = apply(b$t, 2, quantile, probs = 0.975, na.rm = TRUE),
                           failed_resamples = colSums(!is.finite(b$t)), row.names = NULL)
  }
  results <- bind_rows(res)
  rr <- results %>% dplyr::filter(estimand %in% c("logRR", "logRR_comp")) %>%
    dplyr::mutate(estimand = ifelse(estimand == "logRR", "RR", "RR_comp"),
                  estimate = exp(estimate), lower = exp(lower), upper = exp(upper))
  results_all <- bind_rows(results %>% dplyr::filter(!estimand %in% c("logRR", "logRR_comp")), rr)

  # ---- 50-month risks by arm, weights diagnostic, curves ----
  k  <- length(ends)
  fy <- gformula_fit(dat, "young")
  fr <- gformula_fit(dat, "rs")
  ciaj <- aj_fit(dat)
  tp <- cmprsk::timepoints(ciaj, times = 50)$est
  risks <- tibble(dose = dose, method = rep(names(methods), each = 2), arm = rep(c("Placebo", "Estrogen"), 3),
    prostate_cancer_death = c(tp["0 1", 1], tp["1 1", 1], fy$cifP_placebo[k], fy$cifP_treated[k], fr$cifP_placebo[k], fr$cifP_treated[k]),
    other_cause_death     = c(tp["0 2", 1], tp["1 2", 1], fy$cifO_placebo[k], fy$cifO_treated[k], fr$cifO_placebo[k], fr$cifO_treated[k]))

  wd <- fr$weights %>% dplyr::filter(died == 1) %>% dplyr::group_by(rx) %>%
    dplyr::summarise(deaths = dplyr::n(),
                     observed_prostate = sum(eventType_int == 1), sum_wE = sum(wE),
                     observed_other = sum(eventType_int == 2),    sum_wP = sum(wP), .groups = "drop")
  exp_oc <- fr$weights %>% dplyr::group_by(rx) %>% dplyr::summarise(expected_other_life_table = sum(mx * Timey), .groups = "drop")
  wd <- dplyr::left_join(wd, exp_oc, by = "rx") %>%
    dplyr::mutate(dose = dose, arm = ifelse(rx == 1, "Estrogen", "Placebo"),
                  observed_over_expected_other = observed_other / expected_other_life_table) %>% dplyr::select(dose, arm, everything(), -rx)

  aj_curve <- function(nm, arm, outcome) tibble(dose = dose, month = ciaj[[nm]]$time, cif = ciaj[[nm]]$est, arm = arm, outcome = outcome, method = "Aalen-Johansen")
  g_curve <- function(f, method) bind_rows(
    tibble(dose = dose, month = c(0, ends), cif = c(0, f$cifP_placebo), arm = "Placebo",  outcome = "Prostate cancer death", method = method),
    tibble(dose = dose, month = c(0, ends), cif = c(0, f$cifP_treated), arm = "Estrogen", outcome = "Prostate cancer death", method = method),
    tibble(dose = dose, month = c(0, ends), cif = c(0, f$cifO_placebo), arm = "Placebo",  outcome = "Other-cause death",     method = method),
    tibble(dose = dose, month = c(0, ends), cif = c(0, f$cifO_treated), arm = "Estrogen", outcome = "Other-cause death",     method = method))
  curves <- bind_rows(aj_curve("0 1", "Placebo", "Prostate cancer death"), aj_curve("1 1", "Estrogen", "Prostate cancer death"),
                      aj_curve("0 2", "Placebo", "Other-cause death"),     aj_curve("1 2", "Estrogen", "Other-cause death"),
                      g_curve(fy, "Cause-specific g-formula"), g_curve(fr, "Relative survival g-formula"))
  list(results = results_all, risks = risks, weights = wd, curves = curves)
}

out <- lapply(seq_along(doses), function(i) analyse_dose(doses[i], seed = 20261101 + i, i = i))

results <- bind_rows(lapply(out, `[[`, "results")) %>%
  dplyr::mutate(estimand = factor(estimand, levels = c("RD", "RR", "RD_comp", "RR_comp"),
                                  labels = c("Prostate cancer death: risk difference", "Prostate cancer death: risk ratio",
                                             "Other-cause death: risk difference", "Other-cause death: risk ratio")),
                dose = factor(dose, levels = doses)) %>%
  dplyr::arrange(dose, estimand, method)
risks   <- bind_rows(lapply(out, `[[`, "risks"))
weights <- bind_rows(lapply(out, `[[`, "weights"))
curves  <- bind_rows(lapply(out, `[[`, "curves")) %>% dplyr::mutate(dose = factor(dose, levels = doses))

cat("\n\nEffect of each estrogen dose vs placebo at 50 months (95% percentile bootstrap intervals, B =", B, "):\n")
print(as.data.frame(results %>% dplyr::mutate(across(c(estimate, lower, upper), ~ round(.x, 3)))), row.names = FALSE)
cat("\nRisk at 50 months by arm:\n")
print(as.data.frame(risks %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat("\nRelative survival weights (sum over deaths) against the observed causes of death; expected_other_life_table = approximate\n",
    "number of other-cause deaths if the patients died at the life-table rate:\n", sep = "")
print(as.data.frame(weights %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

p <- ggplot(curves, aes(month, cif, colour = arm, linetype = method)) +
  geom_line(data = dplyr::filter(curves, method != "Aalen-Johansen")) +
  geom_step(data = dplyr::filter(curves, method == "Aalen-Johansen")) +
  facet_grid(outcome ~ dose) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Months since randomisation", y = "Cumulative incidence", colour = "Arm", linetype = "Method") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
print(p)

write_csv(results, "results_empirical.csv")
write_csv(risks, "risks_50m.csv")
write_csv(weights, "weights_diagnostic.csv")
ggsave("cif_empirical.pdf", p, width = 12, height = 6)
saveRDS(list(results = results, risks = risks, weights = weights, curves = curves), "empirical_results.rds")
message("Saved results_empirical.csv, risks_50m.csv, weights_diagnostic.csv, cif_empirical.pdf, empirical_results.rds")
