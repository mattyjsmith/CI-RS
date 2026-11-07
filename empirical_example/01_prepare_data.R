################################################################################
# Byar (1980) prostate cancer trial: analysis data sets for the empirical example.
# Each estrogen dose (0.2, 1.0 and 5.0 mg) is compared with the SAME placebo group (n = 127), administratively censored at
# 50 months (deaths after 50 months are alive at 50). One file per comparison: analysis_data_0.2mg.csv, _1.0mg.csv, _5.0mg.csv.
# Reads ../data_50months.csv (made by "Obtaining data.R": Byar data + US Human Mortality Database rates, mx).
# NOTE: the old dataRed.csv holds placebo vs the 0.2 mg arm only.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("empirical_example"))
rm(list = ls())
suppressPackageStartupMessages({ library(dplyr); library(readr) })

doses  <- c("0.2 mg estrogen", "1.0 mg estrogen", "5.0 mg estrogen")
tag_of <- function(dose) gsub(" ", "", sub(" estrogen$", "", dose))
raw    <- read_csv("data/data_50months.csv", show_col_types = FALSE)

prepare <- function(treated_arm) {
  raw %>%
    filter(Trt %in% c("placebo", treated_arm)) %>%
    mutate(
      rx            = as.integer(Trt == treated_arm),
      allCause      = as.integer(status50 != "alive"),
      eventType_int = case_when(status50 == "alive" ~ 0L,
                                status50 == "dead - prostatic ca" ~ 1L,
                                TRUE ~ 2L),        # 0 = alive at 50 months, 1 = prostate cancer death, 2 = other-cause death
      # Month in which the death (or the censoring at 50 months) is recorded: follow-up time rounded up to whole months.
      # Deaths at time 0 (less than a month) go into the first month.
      dtime_int     = pmax(1L, as.integer(ceiling(Time50))),
      # Follow-up time in YEARS for the excess hazard model; half a month for deaths at time 0
      Timey         = pmax(Time50, 0.5) / 12,
      vstatus       = allCause,
      hgBinary      = hg < 12,
      normalAct     = pf == "normal activity",
      ageCat        = as.character(Hmisc::cut2(Age, c(0, 60, 75, 100))),
      age_c         = Age - mean(Age),
      age_c2        = age_c^2
    ) %>%
    mutate(patno = row_number()) %>%
    as.data.frame()
}

needed <- c("rx", "allCause", "eventType_int", "dtime_int", "Timey", "vstatus", "mx", "hgBinary", "normalAct", "ageCat",
            "age_c", "age_c2", "hx", "Stage", "ekg", "sbp", "dbp", "hg", "sz", "sg", "ap", "bm", "wt", "Age")

for (d in doses) {
  dat <- prepare(d)
  stopifnot(all(needed %in% names(dat)), !anyNA(dat[, needed]))     # complete data: no patient is dropped by a model
  cat("\n=== ", d, " vs placebo ===\n", sep = "")
  print(dat %>% count(Trt, eventType_int, name = "n") %>%
          tidyr::pivot_wider(names_from = eventType_int, values_from = n, names_prefix = "type "))
  cat("Patients whose death is recorded at time 0 (put in month 1):", sum(dat$Time50 == 0), "\n")
  write_csv(dat, sprintf("analysis_data_%s.csv", tag_of(d)))
}
message('\nSaved analysis_data_0.2mg.csv, analysis_data_1.0mg.csv, analysis_data_5.0mg.csv')
