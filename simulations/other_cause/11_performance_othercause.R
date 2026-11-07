################################################################################
# Performance (bias, relative bias, empirical SE) of Methods 1 and 2 by scenario. n = 5000.
# The no-effect baseline is taken from the main simulation (Paper folder) if it is available.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "other_cause"))
rm(list = ls())
library(dplyr); library(rsimsum); library(readr); library(tidyr); library(ggplot2); library(knitr)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; summarise <- dplyr::summarise; recode <- dplyr::recode

source("00_scenarios.R")

time_model   <- "spline"     # "spline", "factor" or "linear": must match the run being analysed
results_root <- if (time_model == "linear") "results" else paste0("results_", time_model)

# ---- Results of the five new scenarios ----
read_scenario <- function(sn) {
  files <- list.files(file.path(results_root, sn), pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  do.call(rbind, lapply(files, readRDS))
}
results <- bind_rows(lapply(scenarios$scenario, read_scenario))
truths  <- read_csv("true_values_othercause.csv", show_col_types = FALSE) %>%
  select(Scenario, RD, logRR, RD_comp, logRR_comp)

# ---- Baseline (no effect on other-cause death, ordinary age) from the main simulation ----
main_file  <- if (time_model == "linear") "../main/merged_sim_results_bias_n5000.rds" else "../main/merged_sim_results_bias_n5000.rds"
main_truth <- "../main/true_values.csv"
if (file.exists(main_file) && file.exists(main_truth)) {
  b <- readRDS(main_file)
  b <- b %>% filter(Setting %in% c("Cause-specific", "Relative survival – Perfect world")) %>%
    transmute(Scenario = "base_HR1", SampleSize = 5000, Setting, Rep, RD, logRR, RD_comp, logRR_comp)
  results <- bind_rows(b, results)
  truths  <- bind_rows(read_csv(main_truth, show_col_types = FALSE) %>% mutate(Scenario = "base_HR1"), truths)
} else cat("Main simulation (n = 5000) not found: baseline scenario left out.\n")

truth_long <- truths %>% pivot_longer(c(RD, logRR, RD_comp, logRR_comp), names_to = "estimand", values_to = "theta") %>%
  mutate(estimand = recode(estimand, RD = "Total RD", logRR = "Total logRR", RD_comp = "Competing RD", logRR_comp = "Competing logRR"))

rsim_data <- bind_rows(
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Total RD",        estimate = RD),
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Total logRR",     estimate = logRR),
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Competing RD",    estimate = RD_comp),
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Competing logRR", estimate = logRR_comp)
) %>% left_join(truth_long, by = c("Scenario", "estimand")) %>%
  mutate(method = recode(method, "Cause-specific" = "Cause-specific", "Relative survival – Perfect world" = "Relative survival"))

rsum <- simsum(data = rsim_data, estvarname = "estimate", true = "theta", methodvar = "method", by = c("estimand", "Scenario"))
ss   <- summary(rsum)$summ
write_csv(ss, sprintf("performance_othercause_%s.csv", time_model))

scen_levels <- c("base_HR1", "base_HR1.3", "base_HR1.65", "older_HR1", "older_HR1.3", "older_HR1.65")
scen_labels <- c("Usual age, HR 1.0 (main run)", "Usual age, HR 1.3", "Usual age, HR 1.65",
                 "Older, HR 1.0", "Older, HR 1.3", "Older, HR 1.65")

tab <- ss %>% filter(stat == "rbias") %>%
  mutate(rbias_pct = round(100 * est, 1), mcse_pct = round(100 * mcse, 1)) %>%
  select(estimand, Scenario, method, rbias_pct, mcse_pct) %>%
  arrange(estimand, factor(Scenario, levels = scen_levels), method)
print(kable(tab, format = "markdown"))

plot_df <- ss %>% filter(stat == "rbias") %>%
  mutate(Scenario = factor(Scenario, levels = scen_levels, labels = scen_labels),
         estimand = factor(estimand, levels = c("Total RD", "Total logRR", "Competing RD", "Competing logRR")))
p <- ggplot(plot_df, aes(x = Scenario, y = 100 * est, colour = method, group = method)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = 100 * lower, ymax = 100 * upper), width = 0.25, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ estimand, scales = "free_y") +
  labs(x = NULL, y = "Relative bias (%)", colour = "Method") +
  theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(sprintf("plot_othercause_%s.pdf", time_model), p, width = 11, height = 7)
print(p)

# ---- ABSOLUTE bias (use this for the competing estimands: the true value is close to zero in some scenarios, so a
# ---- relative bias is huge and unstable there) ----
tab_abs <- ss %>% filter(stat == "bias") %>%
  mutate(bias = signif(est, 3), mcse = signif(mcse, 2)) %>%
  select(estimand, Scenario, method, bias, mcse) %>%
  arrange(estimand, factor(Scenario, levels = scen_levels), method)
print(kable(tab_abs, format = "markdown"))

plot_abs <- ss %>% filter(stat == "bias") %>%
  mutate(Scenario = factor(Scenario, levels = scen_levels, labels = scen_labels),
         estimand = factor(estimand, levels = c("Total RD", "Total logRR", "Competing RD", "Competing logRR")))
p2 <- ggplot(plot_abs, aes(x = Scenario, y = est, colour = method, group = method)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ estimand, scales = "free_y") +
  labs(x = NULL, y = "Bias (estimate minus truth)", colour = "Method") +
  theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(sprintf("plot_othercause_absolute_%s.pdf", time_model), p2, width = 11, height = 7)
print(p2)
