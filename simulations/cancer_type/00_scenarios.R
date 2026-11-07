################################################################################
# Cancer aggressiveness scenarios, defined once (sourced by 08, 10, 11 and the tests).
# lambda    : Weibull scale of the cancer death hazard (the shape stays 0.7). The cumulative hazard is lambda * t^0.7 * exp(linear predictor),
#             so a smaller lambda is a more indolent cancer. 0.25 is the main simulation.
# age_shift : years added to every patient's age (10 = the older population, in which other-cause death is common, so cancer
#             is a minority of deaths when it is indolent).
# Targets for the untreated 2-year cancer cumulative incidence: aggressive about 35%, medium about 12%, indolent about 3%.
# The lambda values for medium and indolent are STARTING POINTS: check them with 08_truth_cancertype.R (n_total <- 1000000) and
# adjust before the full runs. The "aggressive" scenario is the older_HR1 scenario of the OtherCause folder (same data and seeds).
################################################################################
scenarios <- data.frame(
  scenario  = c("aggressive", "medium", "indolent"),
  lambda    = c(0.25,         0.07,     0.017),
  age_shift = c(10,           10,       10),
  stringsAsFactors = FALSE
)
