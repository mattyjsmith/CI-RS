# 05_helpers.R
# Helper function to compute cumulative incidence from predicted hazards
# Works for cause-specific hazards from pooled logistic regression models

calculateCumInc <- function(inputData, timepts=cutTimes, competing = FALSE) {
  # inputData must contain:
  # - one row per person per time point
  # - columns: patno, dtime, hazardP (cancer), hazardO (other), s (joint survival prob)
  
  # Create a matrix to hold cumulative incidence at each time for each subject
  n_time <- length(timepts)
  n_subjects <- length(unique(inputData$patno))
  cumulativeIncidence <- matrix(NA, nrow = n_time, ncol = n_subjects)
  
  # Initial hazard at t = 0
  if (!competing) {
    cumulativeIncidence[1, ] <- with(inputData[inputData$dtime == 0, ],
                                     hazardP * (1 - hazardO))
  } else {
    cumulativeIncidence[1, ] <- inputData[inputData$dtime == 0, ]$hazardO
  }
  
  # Compute cumulative survival up to each time point
  survival_matrix <- t(aggregate(s ~ patno, data = inputData, FUN = cumprod)$s)
  
  for (i in 2:n_time) {
    idx <- inputData$dtime == (i - 1)
    subP <- inputData[idx, ]$hazardP
    subO <- inputData[idx, ]$hazardO
    if (!competing) {
      cumulativeIncidence[i, ] <- subP * (1 - subO) * survival_matrix[i - 1, ]
    } else {
      cumulativeIncidence[i, ] <- subO * survival_matrix[i - 1, ]
    }
  }
  
  # Return mean cumulative incidence across individuals at each time point
  rowMeans(apply(cumulativeIncidence, 2, cumsum))
}