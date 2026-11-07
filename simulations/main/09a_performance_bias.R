################################################################################
# Performance measures for the BIAS run (main results): bias, relative bias and empirical SE,
# by method and sample size. Truth: true_values.csv (counterfactual outcomes, 08_true_values.R)
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "main"))
rm(list = ls())

library(dplyr)
library(rsimsum)
library(readr)
library(tidyr)
library(ggplot2)
library(knitr)

# Make sure the dplyr versions are used (MASS and plyr, loaded by 01_setup_packages.R, mask some of them)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; summarise <- dplyr::summarise; recode <- dplyr::recode

# ---- Load the merged results (made by 07b_merge_results.R with run_type <- "bias") ----
sample_sizes <- c(1000, 3000, 5000)
sample_sizes <- sample_sizes[file.exists(sprintf("merged_sim_results_bias_n%d.rds", sample_sizes))]   # those available so far
cat("Sample sizes found:", sample_sizes, "\n")
results <- bind_rows(lapply(sample_sizes, function(n) {
  r <- readRDS(sprintf("merged_sim_results_bias_n%d.rds", n))
  r$SampleSize <- n
  r
}))

# ---- True values ----
truths <- read_csv("true_values.csv") %>%
  pivot_longer(everything(), names_to = "estimand", values_to = "theta") %>%
  mutate(estimand = recode(estimand,
                           RD         = "Total RD",
                           logRR      = "Total logRR",
                           RD_comp    = "Competing RD",
                           logRR_comp = "Competing logRR"))

# ---- Long format ----
rsim_data <- bind_rows(
  results %>% transmute(sim = Rep, method = Setting, estimand = "Total RD",        SampleSize, estimate = RD),
  results %>% transmute(sim = Rep, method = Setting, estimand = "Total logRR",     SampleSize, estimate = logRR),
  results %>% transmute(sim = Rep, method = Setting, estimand = "Competing RD",    SampleSize, estimate = RD_comp),
  results %>% transmute(sim = Rep, method = Setting, estimand = "Competing logRR", SampleSize, estimate = logRR_comp)
) %>%
  left_join(truths, by = "estimand") %>%
  mutate(method = recode(method,
                         "Cause-specific" = "1",
                         "Relative survival – Perfect world" = "2",
                         "Relative survival – Imperfect world (without Dep in EHM)" = "3",
                         "Relative survival – Imperfect world (without Dep in LT)" = "4",
                         "Relative survival – Imperfect world (without Dep in EHM or LT)" = "5"))

# ---- rsimsum (no standard errors: bias measures only) ----
rsum <- simsum(
  data       = rsim_data,
  estvarname = "estimate",
  true       = "theta",
  methodvar  = "method",
  by         = c("estimand", "SampleSize")
)
sumstats <- summary(rsum)
write_csv(sumstats$summ, "performance_bias.csv")

# ---- Table: relative bias (%) with Monte Carlo SE ----
tab <- sumstats$summ %>%
  filter(stat == "rbias") %>%
  mutate(rbias_pct = round(100 * est, 1), mcse_pct = round(100 * mcse, 1)) %>%
  select(estimand, SampleSize, method, rbias_pct, mcse_pct) %>%
  arrange(estimand, SampleSize, method)
print(kable(tab, format = "markdown"))

# ---- Plot: relative bias ----
method_labels <- c("1" = "Cause-specific",
                   "2" = "Relative survival (correct EHM & LT)",
                   "3" = "Relative survival (no Dep in EHM)",
                   "4" = "Relative survival (no Dep in LT)",
                   "5" = "Relative survival (no Dep in EHM or LT)")
plot_df <- sumstats$summ %>%
  filter(stat == "rbias") %>%
  mutate(estimand = factor(estimand, levels = c("Total RD", "Total logRR", "Competing RD", "Competing logRR")),
         method   = factor(method, levels = as.character(1:5), labels = method_labels))

p <- ggplot(plot_df, aes(x = factor(SampleSize), y = 100 * est, colour = method, group = method)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = 100 * lower, ymax = 100 * upper), width = 0.25, position = position_dodge(0.6)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ estimand, scales = "free_y") +
  labs(x = "Sample size", y = "Relative bias (%)", colour = "Method") +
  theme_bw() + theme(legend.position = "bottom")
ggsave("plot_bias_main.pdf", p, width = 10, height = 7)
print(p)
