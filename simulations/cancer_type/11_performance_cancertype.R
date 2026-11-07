################################################################################
# Performance (bias, relative bias) of Methods 1 and 2 by cancer aggressiveness. n = 5000.
# Aggressive: the older_HR1 results of the OtherCause folder (spline time model); medium and indolent: results/ here.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "cancer_type"))
rm(list = ls())
library(dplyr); library(rsimsum); library(readr); library(tidyr); library(ggplot2); library(knitr)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; summarise <- dplyr::summarise; recode <- dplyr::recode

read_dir <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(files, readRDS))
}
agg <- read_dir("../other_cause/results_spline/older_HR1") %>% mutate(Scenario = "aggressive")
new <- bind_rows(lapply(c("medium", "indolent"), function(sn) read_dir(file.path("results", sn))))
results <- bind_rows(agg, new)
cat("Replicates (rows) per scenario and method; failed = NA estimates:\n")
print(results %>% group_by(Scenario, Setting) %>% summarise(rows = n(), failed = sum(is.na(RD)), .groups = "drop"))
results <- results %>% filter(!is.na(RD))

truths <- read_csv("true_values_cancertype.csv", show_col_types = FALSE)
cat("\nTrue values and the cancer share of deaths (untreated):\n")
print(as.data.frame(truths %>% select(Scenario, lambda, patients, cancer_CIF_untreated, other_CIF_untreated,
                                      cancer_share_of_deaths_untreated, RD, RD_comp) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
truth_long <- truths %>% select(Scenario, RD, logRR, RD_comp, logRR_comp) %>%
  pivot_longer(c(RD, logRR, RD_comp, logRR_comp), names_to = "estimand", values_to = "theta") %>%
  mutate(estimand = recode(estimand, RD = "Total RD", logRR = "Total logRR", RD_comp = "Competing RD", logRR_comp = "Competing logRR"))

rsim_data <- bind_rows(
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Total RD",        estimate = RD),
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Total logRR",     estimate = logRR),
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Competing RD",    estimate = RD_comp),
  results %>% transmute(sim = Rep, method = Setting, Scenario, estimand = "Competing logRR", estimate = logRR_comp)
) %>% left_join(truth_long, by = c("Scenario", "estimand")) %>%
  mutate(method = recode(method, "Relative survival – Perfect world" = "Relative survival"))

rsum <- simsum(data = rsim_data, estvarname = "estimate", true = "theta", methodvar = "method", by = c("estimand", "Scenario"))
ss   <- summary(rsum)$summ
write_csv(ss, "performance_cancertype.csv")

lv <- c("aggressive", "medium", "indolent")
tab <- ss %>% filter(stat == "rbias") %>%
  mutate(rbias_pct = round(100 * est, 1), mcse_pct = round(100 * mcse, 1)) %>%
  select(estimand, Scenario, method, rbias_pct, mcse_pct) %>% arrange(estimand, factor(Scenario, levels = lv), method)
cat("\nRelative bias (%):\n"); print(kable(tab, format = "markdown"))
tab_abs <- ss %>% filter(stat == "bias") %>%
  mutate(bias = signif(est, 3), mcse = signif(mcse, 2)) %>%
  select(estimand, Scenario, method, bias, mcse) %>% arrange(estimand, factor(Scenario, levels = lv), method)
cat("\nAbsolute bias:\n"); print(kable(tab_abs, format = "markdown"))
emp <- ss %>% filter(stat == "empse") %>% mutate(empse = signif(est, 3)) %>%
  select(estimand, Scenario, method, empse) %>% arrange(estimand, factor(Scenario, levels = lv), method)
cat("\nEmpirical SE of the estimates:\n"); print(kable(emp, format = "markdown"))

mk_plot <- function(stat_name, ylab, scale) {
  d <- ss %>% filter(stat == stat_name) %>%
    mutate(Scenario = factor(Scenario, levels = lv), estimand = factor(estimand, levels = c("Total RD", "Total logRR", "Competing RD", "Competing logRR")))
  ggplot(d, aes(x = Scenario, y = scale * est, colour = method, group = method)) +
    geom_point(position = position_dodge(0.4)) +
    geom_errorbar(aes(ymin = scale * lower, ymax = scale * upper), width = 0.25, position = position_dodge(0.4)) +
    geom_hline(yintercept = 0, linetype = 2) + facet_wrap(~ estimand, scales = "free_y") +
    labs(x = "Cancer aggressiveness", y = ylab, colour = "Method") + theme_bw() + theme(legend.position = "bottom")
}
p1 <- mk_plot("rbias", "Relative bias (%)", 100); p2 <- mk_plot("bias", "Bias (estimate minus truth)", 1)
ggsave("plot_cancertype_relative.pdf", p1, width = 10, height = 7); ggsave("plot_cancertype_absolute.pdf", p2, width = 10, height = 7)
print(p1); print(p2)
