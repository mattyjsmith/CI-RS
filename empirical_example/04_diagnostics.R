################################################################################
# Does the life table describe these patients? Observed other-cause deaths against the number expected from the life table
# (approximately: expected rate at the exit age/year times follow-up time), overall and by baseline subgroup, for each arm.
# Subgroups are small: read the pattern, not single numbers.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("empirical_example"))
rm(list = ls())
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

raw <- read_csv("data/data_50months.csv", show_col_types = FALSE)
dat <- raw %>%
  dplyr::mutate(
    other_death    = status50 != "alive" & status50 != "dead - prostatic ca",
    expected_other = mx * Time50 / 12,
    `CVD history`        = ifelse(hx == 1, "yes", "no"),
    `ECG`                = ifelse(ekg == "normal", "normal", "abnormal"),
    `Performance status` = ifelse(pf == "normal activity", "normal activity", "limited"),
    `Stage`              = paste("stage", Stage),
    `Age`                = as.character(cut(Age, c(0, 60, 75, 120), right = FALSE, labels = c("<60", "60-74", "75+")))
  )

all_arms <- dplyr::bind_rows(dat, dat %>% dplyr::mutate(Trt = "All four arms"))

overall <- all_arms %>% dplyr::group_by(Trt) %>%
  dplyr::summarise(n = dplyr::n(), observed_other = sum(other_death), expected_other = sum(expected_other),
                   observed_over_expected = observed_other / expected_other, .groups = "drop")
cat("Other-cause deaths: observed against expected from the life table, by arm:\n")
print(as.data.frame(overall %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

by_sub <- all_arms %>%
  tidyr::pivot_longer(c(`CVD history`, `ECG`, `Performance status`, `Stage`, `Age`), names_to = "subgroup", values_to = "level") %>%
  dplyr::group_by(Trt, subgroup, level) %>%
  dplyr::summarise(n = dplyr::n(), observed_other = sum(other_death), expected_other = sum(expected_other),
                   observed_over_expected = observed_other / expected_other, .groups = "drop")
for (arm in c("placebo", "All four arms")) {
  cat("\n", arm, ": observed / expected other-cause deaths by subgroup\n", sep = "")
  print(as.data.frame(by_sub %>% dplyr::filter(Trt == arm) %>% dplyr::select(-Trt) %>%
                        dplyr::mutate(across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)
}
write_csv(overall, "observed_vs_expected_overall.csv")
write_csv(by_sub,  "observed_vs_expected_by_subgroup.csv")
message("Saved observed_vs_expected_overall.csv and observed_vs_expected_by_subgroup.csv")
