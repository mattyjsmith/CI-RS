################################################################################
# Title: Performance Measures Using rsimsum
# Author: Matthew J. Smith
# Note: Truth values derived from large simulated dataset (boot_function_M1)
################################################################################

# ---- Setup ----
setwd("~/Playground")
rm(list = ls())

# ---- Load libraries ----
library(dplyr)
library(rsimsum)
library(readr)
library(knitr)
library(tidyr)
library(ggplot2)

# ---- Load results ----
results <- readRDS("results.rds")

# ---- Step 1: Load true values from large dataset ----
truths <- read_csv("true_values_from_large_data.csv") %>%
  pivot_longer(everything(), names_to = "estimand", values_to = "theta") %>%
  mutate(estimand = recode(estimand,
                           RD         = "Total RD",
                           logRR      = "Total logRR",
                           RD_comp    = "Competing RD",
                           logRR_comp = "Competing logRR"
  ))

# ---- Step 2: Reshape results to long format ----
rsim_data <- bind_rows(
  results %>%
    transmute(
      sim = Rep, method = Setting, estimand = "Total RD", SampleSize,
      estimate = RD, se = SE_RD, lower = CI_RD_L, upper = CI_RD_U
    ),
  results %>%
    transmute(
      sim = Rep, method = Setting, estimand = "Total logRR", SampleSize,
      estimate = logRR, se = SE_logRR, lower = CI_logRR_L, upper = CI_logRR_U
    ),
  results %>%
    transmute(
      sim = Rep, method = Setting, estimand = "Competing RD", SampleSize,
      estimate = RD_comp, se = SE_RDcomp, lower = CI_RDcomp_L, upper = CI_RDcomp_U
    ),
  results %>%
    transmute(
      sim = Rep, method = Setting, estimand = "Competing logRR", SampleSize,
      estimate = logRR_comp, se = SE_logRRcomp, lower = CI_logRRcomp_L, upper = CI_logRRcomp_U
    )
) %>%
  left_join(truths, by = "estimand")

# Save the merged data
saveRDS(rsim_data, file = "rsim_data.rds")
  # rsim_data <- readRDS("rsim_data.rds")

# Relabel long method names with numeric codes
rsim_data <- rsim_data %>%
  mutate(method = recode(method,
                         "Cause-specific" = "1",
                         "Relative survival – Perfect world" = "2",
                         "Relative survival – Imperfect world (without Dep in EHM)" = "3",
                         "Relative survival – Imperfect world (without Dep in LT)" = "4",
                         "Relative survival – Imperfect world (without Dep in EHM or LT)" = "5"
  ))

# ---- Step 3: Run rsimsum ----
rsum <- simsum(
  data = rsim_data,
  estvarname = "estimate",
  se = "se",
  true = "theta",
  methodvar = "method",
  by = c("estimand", "SampleSize")  
)

# ---- Step 4: Output summary ----
sumstats <- summary(rsum)
print(sumstats)

# Save results
write_csv(sumstats$summ, "performance_summary_rsimsum.csv")
saveRDS(sumstats, "performance_summary_rsimsum.Rds")
write_csv(summary(rsum)$summ, "performance_summary_with_samplesizes.csv")

# Markdown-friendly table
cat("\n\n--- Performance Summary Table ---\n")
print(kable(sumstats$summ, format = "markdown"))


# ---- Step 5: Plotting ----
df_summary <- sumstats$summ

# Reorder estimands
df_summary <- df_summary %>%
  mutate(estimand = factor(estimand,
                           levels = c("Total RD", "Total logRR", "Competing RD", "Competing logRR")
  ))

# Map method codes to descriptive labels and order
method_labels <- c(
  "1" = "Cause-specific",
  "2" = "Relative survival (Correct EHM & LT)",
  "3" = "Relative survival (no Dep in EHM)",
  "4" = "Relative survival (no Dep in LT)",
  "5" = "Relative survival (no Dep in EHM or LT)"
)

df_summary <- df_summary %>%
  mutate(
    method = factor(method,
                    levels = as.character(1:5),
                    labels = method_labels
    )
  )

# ---- Relative Bias Plot (Shape-only aesthetic) ----
plot_bias <- df_summary %>%
  filter(stat == "rbias", SampleSize %in% c(1000, 3000, 5000)) %>%
  ggplot(aes(x = estimand, y = est, shape = method, group = method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~SampleSize, ncol = 1) +  
  labs(
    x = "Estimand", y = "Relative Bias", shape = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    axis.title       = element_text(size = 13),  # bigger axis titles
    strip.text       = element_text(size = 13, face = "bold"), # bigger facet labels
    legend.position  = "bottom",
    legend.justification = "center",
    legend.box       = "vertical",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8.5),
    legend.title     = element_blank(),
    plot.margin      = margin(t = 5, r = 5, b = 30, l = 5)
  ) +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))

# ---- Coverage Plot (Shape-only aesthetic) ----
plot_coverage <- df_summary %>%
  filter(stat == "cover", SampleSize %in% c(1000, 3000, 5000)) %>%
  ggplot(aes(x = estimand, y = est, shape = method, group = method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  facet_wrap(~SampleSize, ncol = 1) +  
  scale_y_continuous(limits = c(0.35, 1.00), breaks = seq(0.35, 1.00, by = 0.1)) +
  labs(
    x = "Estimand", y = "Coverage", shape = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    axis.title       = element_text(size = 13),
    strip.text       = element_text(size = 13, face = "bold"),
    legend.position  = "bottom",
    legend.justification = "center",
    legend.box       = "vertical",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8.5),
    legend.title     = element_blank(),
    plot.margin      = margin(t = 5, r = 5, b = 30, l = 5)
  ) +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))

# ---- Bias-Eliminated Coverage Plot (Shape-only aesthetic) ----
plot_be_coverage <- df_summary %>%
  filter(stat == "becover", SampleSize %in% c(1000, 3000, 5000)) %>%
  ggplot(aes(x = estimand, y = est, shape = method, group = method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  facet_wrap(~SampleSize, ncol = 1) +  
  scale_y_continuous(limits = c(0.85, 1.00), breaks = seq(0.85, 1.00, by = 0.05)) +
  labs(
    x = "Estimand", y = "Bias-Eliminated Coverage", shape = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    axis.title       = element_text(size = 13),
    strip.text       = element_text(size = 13, face = "bold"),
    legend.position  = "bottom",
    legend.justification = "center",
    legend.box       = "vertical",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8.5),
    legend.title     = element_blank(),
    plot.margin      = margin(t = 5, r = 5, b = 30, l = 5)
  ) +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))

# ---- Relative Bias Plot (Shape-only aesthetic) ----
plot_relerror <- df_summary %>%
  filter(stat == "relerror", SampleSize %in% c(1000, 3000, 5000)) %>%
  ggplot(aes(x = estimand, y = est, shape = method, group = method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~SampleSize, ncol = 1) +  
  labs(
    x = "Estimand", y = "Relative Error", shape = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5),
    axis.title       = element_text(size = 13),  # bigger axis titles
    strip.text       = element_text(size = 13, face = "bold"), # bigger facet labels
    legend.position  = "bottom",
    legend.justification = "center",
    legend.box       = "vertical",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8.5),
    legend.title     = element_blank(),
    plot.margin      = margin(t = 5, r = 5, b = 30, l = 5)
  ) +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))


# ---- Show and Save Plots ----
print(plot_bias)
print(plot_coverage)
print(plot_be_coverage)
print(plot_relerror)

library(ggplot2)

# ---- Save A4 Landscape PDFs ----
# A4 dimensions in inches (landscape)
a4_width  <- 11.69
a4_height <- 8.27

# Relative Bias
ggsave(
  filename = "plot_bias.pdf",
  plot     = plot_bias,
  width    = a4_width,
  height   = a4_height,
  units    = "in"
)

# Coverage
ggsave(
  filename = "plot_coverage.pdf",
  plot     = plot_coverage,
  width    = a4_width,
  height   = a4_height,
  units    = "in"
)

# Bias-Eliminated Coverage
ggsave(
  filename = "plot_be_coverage.pdf",
  plot     = plot_be_coverage,
  width    = a4_width,
  height   = a4_height,
  units    = "in"
)

# Relative error
ggsave(
  filename = "plot_relerror.pdf",
  plot     = plot_relerror,
  width    = a4_width,
  height   = a4_height,
  units    = "in"
)








##########################################################################################
# Appendix
##########################################################################################
# ---- Check intervals are based on MCSE ----

# Function to check CI = est ± 1.96 * mcse
verify_mcse_ci <- function(df, tolerance = 1e-6) {
  # Make sure dplyr is loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The dplyr package is required but not installed.")
  }
  
  # Filter rows with non-NA mcse and compute expected CI bounds
  df_check <- df %>%
    dplyr::filter(!is.na(mcse)) %>%
    dplyr::mutate(
      lower_calc = est - 1.96 * mcse,
      upper_calc = est + 1.96 * mcse,
      lower_match = abs(lower - lower_calc) < tolerance,
      upper_match = abs(upper - upper_calc) < tolerance,
      ci_match = lower_match & upper_match
    )
  
  # Report summary
  n_checked <- nrow(df_check)
  n_passed <- sum(df_check$ci_match, na.rm = TRUE)
  
  cat(sprintf(
    "✅ Verified %d rows with non-NA MCSE:\n  - %d matched within tolerance (%.0e).\n  - %d failed.\n\n",
    n_checked, n_passed, tolerance, n_checked - n_passed
  ))
  
  # Return mismatches (if any)
  mismatches <- df_check %>% dplyr::filter(!ci_match)
  if (nrow(mismatches) > 0) {
    cat("⚠️ Mismatches found:\n")
    print(mismatches %>%
            dplyr::select(stat, est, mcse, lower, lower_calc, upper, upper_calc))
  } else {
    cat("🎉 All confidence intervals match expected values.\n")
  }
  
  invisible(df_check)
}

# Run the function
verify_mcse_ci(df_summary, tolerance = 1e-3)






