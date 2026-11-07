################################################################################
# Does the time model matter? Relative bias of Methods 1 and 2 at n = 5000 for linear time,
# interval-specific intercepts and a spline, and the paired differences between them. Truth: ../main/true_values.csv
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "time_model"))
rm(list = ls())
library(dplyr); library(rsimsum); library(readr); library(tidyr); library(ggplot2); library(knitr)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; summarise <- dplyr::summarise; recode <- dplyr::recode

read_tm <- function(tm) {
  files <- list.files(file.path("results", tm), pattern = "\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(files, readRDS))
}
results <- bind_rows(lapply(c("linear", "factor", "spline"), read_tm))
cat("Replicates per time model and method:\n"); print(table(results$TimeModel, results$Setting))

truth <- read_csv("../main/true_values.csv", show_col_types = FALSE) %>%
  pivot_longer(everything(), names_to = "estimand", values_to = "theta") %>%
  mutate(estimand = recode(estimand, RD = "Total RD", logRR = "Total logRR", RD_comp = "Competing RD", logRR_comp = "Competing logRR"))

long <- bind_rows(
  results %>% transmute(sim = Rep, TimeModel, method = Setting, estimand = "Total RD",        estimate = RD),
  results %>% transmute(sim = Rep, TimeModel, method = Setting, estimand = "Total logRR",     estimate = logRR),
  results %>% transmute(sim = Rep, TimeModel, method = Setting, estimand = "Competing RD",    estimate = RD_comp),
  results %>% transmute(sim = Rep, TimeModel, method = Setting, estimand = "Competing logRR", estimate = logRR_comp)
) %>% left_join(truth, by = "estimand") %>%
  mutate(method = recode(method, "Relative survival – Perfect world" = "Relative survival"))

rsum <- simsum(data = long, estvarname = "estimate", true = "theta", methodvar = "method", by = c("estimand", "TimeModel"))
ss <- summary(rsum)$summ
write_csv(ss, "performance_timemodel.csv")

tab <- ss %>% filter(stat == "rbias") %>%
  mutate(rbias_pct = round(100 * est, 2), mcse_pct = round(100 * mcse, 2)) %>%
  select(estimand, method, TimeModel, rbias_pct, mcse_pct) %>%
  arrange(estimand, method, factor(TimeModel, levels = c("linear", "factor", "spline")))
cat("\nRelative bias (%):\n"); print(kable(tab, format = "markdown"))

# Paired differences (same data sets): new time model minus linear, in absolute terms
wide <- long %>% select(sim, method, estimand, TimeModel, estimate) %>%
  pivot_wider(names_from = TimeModel, values_from = estimate)
paired <- bind_rows(
  wide %>% mutate(diff = factor - linear, comparison = "interval-specific intercepts minus linear"),
  wide %>% mutate(diff = spline - linear, comparison = "spline minus linear")
) %>% group_by(estimand, method, comparison) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE), se = sd(diff, na.rm = TRUE) / sqrt(sum(!is.na(diff))), .groups = "drop")
cat("\nPaired differences in the estimates (mean and SE over data sets):\n")
print(kable(paired %>% mutate(mean_diff = signif(mean_diff, 3), se = signif(se, 2)), format = "markdown"))

plot_df <- ss %>% filter(stat == "rbias") %>%
  mutate(TimeModel = factor(TimeModel, levels = c("linear", "factor", "spline"),
                            labels = c("Linear time", "Interval-specific intercepts", "Spline (4 df)")),
         estimand = factor(estimand, levels = c("Total RD", "Total logRR", "Competing RD", "Competing logRR")))
p <- ggplot(plot_df, aes(x = TimeModel, y = 100 * est, colour = method, group = method)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = 100 * lower, ymax = 100 * upper), width = 0.25, position = position_dodge(0.4)) +
  geom_hline(yintercept = 0, linetype = 2) + facet_wrap(~ estimand, scales = "free_y") +
  labs(x = NULL, y = "Relative bias (%)", colour = "Method") + theme_bw() + theme(legend.position = "bottom")
ggsave("plot_timemodel.pdf", p, width = 10, height = 7)
print(p)
