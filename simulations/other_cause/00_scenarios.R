################################################################################
# The scenarios, defined once (sourced by 08, 10, 11 and the tests)
# HR_other  : hazard ratio of treatment on OTHER-CAUSE death (1 = no effect, as in the main simulation)
# age_shift : years added to every patient's age at diagnosis (0 = the main simulation; 10 = an older population,
#             with more background mortality). Age is centred on the shifted mean, so the cancer hazard and the
#             treatment assignment are unchanged: only the other-cause mortality differs.
# The baseline scenario (HR 1, age_shift 0) is the main simulation at n = 5000 (Paper folder), not repeated here.
################################################################################
scenarios <- data.frame(
  scenario  = c("base_HR1.3", "base_HR1.65", "older_HR1", "older_HR1.3", "older_HR1.65"),
  HR_other  = c(1.3,          1.65,          1,           1.3,           1.65),
  age_shift = c(0,            0,             10,          10,            10),
  stringsAsFactors = FALSE
)
