## =========================================================
## Byar (1980) prostate data + HMD life tables (AGE-OF-EVENT & YEAR-OF-EVENT)
## - Truncate follow-up at 50 months; deaths >50m -> alive (censored) at 50m
## =========================================================

## --- Setup ---
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("empirical_example"))
dir.create("data", showWarnings = FALSE)   # USHMD.txt (Human Mortality Database) must already be in empirical_example/data

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(stringr)
})

## --- 1) Import Byar1980 data ---
byar_url <- "https://raw.githubusercontent.com/IowaBiostat/data-sets/main/Byar1980/Byar1980.txt"
byar_raw <- read_tsv(byar_url, show_col_types = FALSE)
readr::write_csv(byar_raw, "data/rawdata.csv")

## Use the exact Status levels present in the file
status_levels <- c(
  "alive",
  "dead - cerebrovascular",
  "dead - heart or vascular",
  "dead - other ca",
  "dead - other specific non-ca",
  "dead - prostatic ca",
  "dead - pulmonary embolus",
  "dead - respiratory disease",
  "dead - unknown cause",
  "dead - unspecified non-ca"
)

## --- 2) Prepare Byar data: dates, flags, ages ---
byar <- byar_raw %>%
  mutate(
    Date        = as.Date(Date),                 # start-of-follow-up (calendar)
    Time_months = as.numeric(Time),
    
    # End-of-follow-up (event/censor) date (estimated)
    end_date_est = if (all(Time_months == floor(Time_months), na.rm = TRUE)) {
      Date %m+% months(Time_months)             # exact month arithmetic for integer months
    } else {
      Date + days(round(Time_months * 365.2425 / 12))  # proportional-day approx for fractional months
    },
    
    # Any event?
    Status      = as.character(Status),
    event_any   = Status != "alive",
    year_start  = year(Date),
    year_end_est = year(end_date_est),
    
    # Age at event/censor (years); last birthday integer for life-table match
    age_at_event  = as.numeric(Age) + Time_months / 12,
    age_death_int = floor(age_at_event),
    
    # Useful factors
    Stage       = factor(Stage),
    Trt         = factor(Trt),
    status_cause = factor(Status, levels = status_levels)
  )

## --- 3) Import HMD life table (local file) ---
hmd <- tryCatch(
  readr::read_table("data/USHMD.txt", skip = 1, comment = "#", na = c(".", "NA"), show_col_types = FALSE),
  error = function(e) readr::read_table("data/USHMD.txt", comment = "#", na = c(".", "NA"), show_col_types = FALSE)
)

## Standardize HMD columns
hmd <- hmd %>%
  mutate(
    Year    = as.integer(Year),
    Age_num = as.integer(parse_number(Age))   # e.g., "110+" -> 110
  )

## (Optional) Restrict to relevant years
hmd_range    <- range(hmd$Year, na.rm = TRUE)
years_needed <- range(byar$year_end_est, na.rm = TRUE)
message("HMD year range: ", paste(hmd_range, collapse = "–"),
        " | Byar end years observed: ", paste(years_needed, collapse = "–"))

hmd_sub <- hmd %>% filter(Year >= years_needed[1], Year <= years_needed[2])

## --- 4) Merge expected mortality (original, untruncated) ---
byar_joined_all <- byar %>%
  mutate(
    death_year_est = if_else(event_any, year(end_date_est), NA_integer_),
    hmd_year = coalesce(death_year_est, year_end_est),
    hmd_age  = as.integer(age_death_int)
  ) %>%
  left_join(
    hmd_sub %>% select(Year, Age_num, mx, qx, ex),
    by = c("hmd_year" = "Year", "hmd_age" = "Age_num")
  )

byar_joined_deaths <- byar_joined_all  # keep name if used downstream

## --- 5) Truncate follow-up at 50 months and recode status ---
cut_m <- 50

byar <- byar %>%
  mutate(
    # Truncated follow-up time (months)
    Time50 = pmin(Time_months, cut_m),
    
    # Updated status: deaths after 50 months become "alive" at 50 (censored)
    status50_chr = dplyr::case_when(
      Time_months > cut_m ~ "alive",
      TRUE                ~ Status
    ),
    status50 = factor(status50_chr, levels = status_levels),
    
    # Event indicator under the 50-month regime
    event_any50 = status50 != "alive",
    
    # Truncated end date (calendar) for joins
    end_date_50 = if (all(Time50 == floor(Time50), na.rm = TRUE)) {
      Date %m+% months(Time50)
    } else {
      Date + days(round(Time50 * 365.2425 / 12))
    },
    year_end_50 = lubridate::year(end_date_50),
    
    # Age at truncated event/censor (years; last birthday integer)
    age_at_event50  = as.numeric(Age) + Time50 / 12,
    age_death_int50 = floor(age_at_event50),
    
    # Death year under truncation (only if still an event by 50 months)
    death_year_50 = dplyr::if_else(event_any50, lubridate::year(end_date_50), NA_integer_)
  )

## --- 6) Merge expected mortality for truncated data (50 months) ---
byar_joined_50 <- byar %>%
  mutate(
    hmd_year_50 = dplyr::coalesce(death_year_50, year_end_50),
    hmd_age_50  = as.integer(age_death_int50)
  ) %>%
  left_join(
    hmd_sub %>% dplyr::select(Year, Age_num, mx, qx, ex),
    by = c("hmd_year_50" = "Year", "hmd_age_50" = "Age_num")
  )

## --- 7) Diagnostics ---
n_deaths_orig    <- byar %>% filter(event_any) %>% nrow()
n_unmatched_orig <- sum(is.na(byar_joined_deaths$mx))
message("Original: deaths total = ", n_deaths_orig,
        " | Unmatched after HMD join = ", n_unmatched_orig)

n_recode_alive <- byar %>% filter(Time_months > cut_m & event_any) %>% nrow()
n_events_by_50 <- byar %>% filter(event_any50) %>% nrow()
n_unmatched_50 <- sum(is.na(byar_joined_50$mx))
message("Truncation @50m: recoded deaths>50m to alive = ", n_recode_alive,
        " | deaths by 50m (any cause) = ", n_events_by_50,
        " | Unmatched after HMD join = ", n_unmatched_50)

# Peek at merged fields (truncated)
print(
  byar_joined_50 %>%
    select(Stage, Trt, Age, Time_months, Time50, status50, event_any50,
           end_date_50, age_at_event50, age_death_int50, mx, qx, ex) %>%
    head(10)
)

## --- 8) Save outputs ---
readr::write_csv(byar_joined_deaths, "data/data.csv")          # original follow-up
readr::write_csv(byar_joined_50,     "data/data_50months.csv") # truncated at 50 months