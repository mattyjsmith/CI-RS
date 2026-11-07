# 05_helpers_CL.R  (corrected)
# Cumulative incidence from predicted discrete-time hazards (pooled logistic regression models).
#
# inputData : one row per person per prediction time, with columns patno, dtime, hazardP (cancer death)
#             and hazardO (other-cause death). The hazards are those of the interval ENDING at dtime.
# timepts   : not used (kept so that old calls still work); the intervals are those in dtime.
# competing : FALSE = cumulative incidence of cancer death; TRUE = of other-cause death.
#
# Within an interval an other-cause death is taken to happen first, then a cancer death among the
# survivors. This matches how the cancer hazard is fitted (people with an other-cause death in the
# interval are left out of the cancer model).
# Returns the mean cumulative incidence over people at the end of each interval.
calculateCumInc <- function(inputData, timepts = NULL, competing = FALSE) {
  # Put rows in patno / time order (so any earlier re-ordering, e.g. by merge(), cannot mis-pair people)
  inputData  <- inputData[order(inputData$patno, inputData$dtime), ]
  n_time     <- length(unique(inputData$dtime))   # number of intervals (timepts is not needed)
  n_subjects <- length(unique(inputData$patno))
  stopifnot(nrow(inputData) == n_time * n_subjects)

  hP <- matrix(inputData$hazardP, nrow = n_time)   # one column per person, one row per interval
  hO <- matrix(inputData$hazardO, nrow = n_time)

  surv_end  <- apply((1 - hP) * (1 - hO), 2, cumprod)        # survival to the END of each interval
  surv_prev <- rbind(1, surv_end[-n_time, , drop = FALSE])    # survival to the START of each interval

  inc <- if (!competing) hP * (1 - hO) * surv_prev else hO * surv_prev
  rowMeans(apply(inc, 2, cumsum))
}
