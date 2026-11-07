################################################################################
# Tests of the corrected empirical pipeline, for each dose (a few minutes, point estimates only):
# 1. The relative survival machinery given the TRUE cause of death must reproduce the cause-specific method exactly.
# 2. The cause-specific g-formula should be close to the non-parametric (Aalen-Johansen) estimate in these randomised data.
# 3. Shows the relative survival estimate and how its weights compare with the observed causes of death.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("empirical_example"))
rm(list = ls())
suppressPackageStartupMessages({ library(dplyr); library(readr); library(survival); library(mexhaz); library(cmprsk) })
source("02_functions.R")

ok <- TRUE
for (tag in c("0.2mg", "1.0mg", "5.0mg")) {
  dat <- as.data.frame(read_csv(sprintf("analysis_data_%s.csv", tag), show_col_types = FALSE))
  young  <- gformula_est(dat, "young")
  rstrue <- gformula_est(dat, "rs_true")
  aj     <- aj_est(dat)
  rs     <- gformula_est(dat, "rs")

  cat("\n=====", tag, "vs placebo =====\n")
  cat("Estimates at 50 months (RD and log RR for prostate cancer death; RD and log RR for the competing event):\n")
  print(round(rbind("Aalen-Johansen" = aj, "Cause-specific g-formula" = young,
                    "Relative survival, TRUE cause of death" = rstrue, "Relative survival, estimated weights" = rs), 4))
  same <- all(abs(young - rstrue) < 1e-6); ok <- ok && same
  cat("Relative survival with the true cause of death equals the cause-specific method (to 1e-6):", same, "\n")
  cat("Largest difference between the cause-specific g-formula and Aalen-Johansen (RD):",
      round(max(abs(young[c("RD", "RD_comp")] - aj[c("RD", "RD_comp")])), 4), "\n")
  w <- gformula_fit(dat, "rs")$weights %>% dplyr::filter(died == 1)
  cat("Among the", nrow(w), "deaths: observed prostate cancer deaths", sum(w$eventType_int == 1),
      "| sum of cancer weights", round(sum(w$wE), 1), "| observed other-cause deaths", sum(w$eventType_int == 2),
      "| sum of other-cause weights", round(sum(w$wP), 1), "\n")
}
cat("\nAll equality checks passed:", ok, "\n")
