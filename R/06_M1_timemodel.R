# GENERATED from the Paper folder's function; the only change is the time model of the logistic regressions (argument time_model).
# 06_bootstrap_function_M1.R
# Cause-specific pooled logistic regression (Young) with the parametric g-formula.
# Returns: total RD, total logRR, competing RD, competing logRR at 2 years.
# Corrections: see README_paper.md

boot_function_M1_tm <- function(data, indices, time_model = "linear"){
  
  dt <- data[indices,]
  dt$patno <- 1:nrow(dt)
  # Time model of the pooled logistic regressions: linear (as in the main simulation), interval-specific
  # intercepts (one per month), or a natural cubic spline with 4 df
  time_term <- switch(time_model,
                      linear = "dtime",
                      factor = "factor(dtime)",
                      spline = "ns(dtime, df = 4)")
  ends <- round(cutTimes + (cutTimes[2] - cutTimes[1]), 10)   # end time of each real (monthly) interval
  
  # Administrative censoring at 2 years since diagnosis
  dt$cause   <- ifelse(dt$Timey >= 2, 0, dt$cause)
  dt$vstatus <- ifelse(dt$Timey >= 2, 0, dt$vstatus)
  dt$Timey   <- ifelse(dt$Timey >= 2, 2, dt$Timey)
  
  dt$allCause  <- dt$vstatus != 0
  dt$eventType <- as.factor(dt$cause)
  
  dt$Tstart <- -0.01
  dt$dtime  <- ceiling(dt$Timey * 12)       # month of death or censoring
  
  # One row per person per month at risk
  longdata <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="allCause")
  
  longdata$CancerDeath <- longdata$allCause==1 & longdata$eventType==1
  longdata$otherDeath  <- longdata$allCause==1 & longdata$eventType==2
  
  # Cancer model: people with an other-cause death in the interval are not at risk of a cancer death
  # (other-cause death acts first within an interval, as assumed in calculateCumInc)
  longdata$CancerDeath[longdata$otherDeath==1] <- NA
  
  baseline <- longdata[longdata$dtime==0,]    # one zero-length dummy row per person: baseline covariates only
  longdata <- longdata[longdata$dtime>0,]     # the models are fitted on real person-intervals only
  n.obs    <- length(unique(longdata$patno))
  
  plrFitP <- glm(as.formula(paste("CancerDeath ~", time_term, "+ agecr + IsexH + trt + dep")), 
                 data = longdata, family=binomial())
  plrFitO <- glm(as.formula(paste("otherDeath ~", time_term, "+ agecr + IsexH + trt + dep")), 
                 data = longdata, family=binomial())
  
  treated <- baseline[rep(1:n.obs, each=length(ends)),]
  treated$dtime <- rep(ends, n.obs)
  treated$trt <- 1
  treated$hazardP <- predict(plrFitP, newdata = treated, type = 'response')
  treated$hazardO <- predict(plrFitO, newdata = treated, type = 'response')
  
  placebo <- baseline[rep(1:n.obs, each=length(ends)),]
  placebo$dtime <- rep(ends, n.obs)
  placebo$trt <- 0
  placebo$hazardP <- predict(plrFitP, newdata = placebo, type = 'response')
  placebo$hazardO <- predict(plrFitO, newdata = placebo, type = 'response')
  
  cumIncTreatedRS <- calculateCumInc(treated)
  cumIncPlaceboRS <- calculateCumInc(placebo)
  log.gcomptotrr <- log(cumIncTreatedRS[length(ends)] / cumIncPlaceboRS[length(ends)])
  gcomptotrd <- cumIncTreatedRS[length(ends)] - cumIncPlaceboRS[length(ends)]
  
  cumIncTreated5comp <- calculateCumInc(treated, competing = TRUE)
  cumIncPlacebo5comp <- calculateCumInc(placebo, competing = TRUE)
  log.gcompcomprr <- log(cumIncTreated5comp[length(ends)] / cumIncPlacebo5comp[length(ends)])
  gcompcomprd <- cumIncTreated5comp[length(ends)] - cumIncPlacebo5comp[length(ends)]
  
  return(c(gcomptotrd, log.gcomptotrr, gcompcomprd, log.gcompcomprr))
}
