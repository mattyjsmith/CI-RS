################################################################################
# Functions for the empirical example (corrected). The same corrections as in the simulation study:
#  * the models are fitted on real person-months only (the zero-length dtime = 0 row is used for baseline covariates only);
#  * hazards are predicted for the interval ENDING at each month (months 1..50) and the cumulative incidence is read at month 50;
#  * the other-cause model of the relative survival method is a fractional-outcome fit (weight 1 on every row);
#  * the cause-specific method leaves people with an other-cause death in a month out of that month's cancer model.
################################################################################

cutTimes <- 0:49              # split points; person-months are labelled by their END month dtime = 1..50
ends     <- cutTimes + 1      # 1..50

# Outcome models: the specification used in Analysis.R / AnalysisRS.R (cubic time + treatment-by-time for cancer death)
rhs_P <- "dtime + dtime2 + dtime3 + normalAct + ageCat + hx + hgBinary + rx + rx1 + rx2"   # cancer death
rhs_O <- "dtime + dtime2 + normalAct + ageCat + hx + hgBinary + rx"                         # other-cause death

# Excess hazard model: the specification used in AnalysisRS.R
ehm_formula <- Surv(Timey, vstatus) ~ rx + age_c + age_c2 + hgBinary + normalAct + hx + Stage + ekg + sbp + dbp +
                                      hg + sz + sg + ap + bm + wt

add_time_terms <- function(d) {
  d$dtime2 <- d$dtime^2
  d$dtime3 <- d$dtime2 * d$dtime
  d$rx1    <- d$dtime  * d$rx
  d$rx2    <- d$dtime2 * d$rx
  d
}

# Cumulative incidence from predicted discrete-time hazards (corrected helper).
# inputData: one row per person per month, columns patno, dtime, hazardP (cancer death), hazardO (other-cause death).
# Within a month an other-cause death is taken to happen first, then a cancer death among the survivors.
calculateCumInc <- function(inputData, competing = FALSE) {
  inputData  <- inputData[order(inputData$patno, inputData$dtime), ]
  n_time     <- length(unique(inputData$dtime))
  n_subjects <- length(unique(inputData$patno))
  stopifnot(nrow(inputData) == n_time * n_subjects)
  hP <- matrix(inputData$hazardP, nrow = n_time)
  hO <- matrix(inputData$hazardO, nrow = n_time)
  surv_end  <- apply((1 - hP) * (1 - hO), 2, cumprod)
  surv_prev <- rbind(1, surv_end[-n_time, , drop = FALSE])
  inc <- if (!competing) hP * (1 - hO) * surv_prev else hO * surv_prev
  rowMeans(apply(inc, 2, cumsum))
}

# Relative survival weights: probability that a death is a cancer death (wE) or an other-cause death (wP)
excess_hazard_weights <- function(dt, mx_mult = 1) {
  dt$mx <- dt$mx * mx_mult      # sensitivity analysis: life-table rates multiplied by mx_mult (1 = as observed)
  # Standardise the continuous covariates: the same model, but a much better-conditioned optimisation (fit and prediction
  # below both use the standardised data)
  for (v in c("age_c", "age_c2", "sbp", "dbp", "hg", "sz", "ap", "wt")) dt[[v]] <- as.numeric(scale(dt[[v]]))
  # mexhaz prints the optimiser's trace ("Iteration limit exceeded...") to the console whatever `verbose` is: capture and discard it
  fit <- NULL
  invisible(capture.output(fit <- suppressWarnings(suppressMessages(
    mexhaz(ehm_formula, data = dt, base = "weibull", expected = "mx", verbose = 0)))))
  eh <- vapply(seq_len(nrow(dt)), function(k) {
    predict(fit, time.pts = dt$Timey[k], data.val = dt[k, , drop = FALSE])$results$hazard
  }, numeric(1))
  den <- eh + dt$mx
  data.frame(excess = eh, wE = eh / den, wP = dt$mx / den)
}

# Parametric g-formula. method:
#   "young"   : cause-specific pooled logistic models (Young et al.)
#   "rs"      : relative survival: weights from the excess hazard model and the life table
#   "rs_true" : the relative survival machinery with the TRUE cause of death as the weights (must equal "young": a check)
gformula_fit <- function(dt, method = c("young", "rs", "rs_true"), mx_mult = 1) {
  method <- match.arg(method)
  dt <- as.data.frame(dt)
  dt$patno <- seq_len(nrow(dt))
  n.obs <- nrow(dt)

  if (method == "rs") {
    w <- excess_hazard_weights(dt, mx_mult)
    dt$wE <- w$wE; dt$wP <- w$wP
  } else {
    dt$wE <- as.numeric(dt$eventType_int == 1)
    dt$wP <- as.numeric(dt$eventType_int == 2)
  }

  # Person-month data
  dt$Tstart <- -0.01
  dt$dtime  <- dt$dtime_int
  long     <- survSplit(data = dt, cut = cutTimes, start = "Tstart", end = "dtime", event = "allCause")
  baseline <- long[long$dtime == 0, ]      # one zero-length dummy row per person: baseline covariates only
  long     <- long[long$dtime > 0, ]       # the models are fitted on real person-months only
  stopifnot(nrow(baseline) == n.obs, max(long$dtime) <= max(ends))
  long <- add_time_terms(long)
  died <- long$allCause == 1

  if (method == "young") {
    long$yP <- as.numeric(died & long$eventType_int == 1)
    long$yO <- as.numeric(died & long$eventType_int == 2)
    long$yP[long$yO == 1] <- NA            # other-cause death in the month: not at risk of a cancer death that month
    fitP <- glm(as.formula(paste("yP ~", rhs_P)), data = long, family = binomial())
    fitO <- glm(as.formula(paste("yO ~", rhs_O)), data = long, family = binomial())
  } else {
    long$wC <- ifelse(died, long$wE, 1)    # weight 1 on every row, except the death row (probability of a cancer death)
    long$yO <- ifelse(died, long$wP, 0)    # expected other-cause death (fractional outcome), weight 1 on every row
    fitP <- suppressWarnings(glm(as.formula(paste("allCause ~", rhs_P)), data = long, weights = wC, family = binomial()))
    fitO <- suppressWarnings(glm(as.formula(paste("yO ~", rhs_O)), data = long, family = quasibinomial()))
  }

  make_grid <- function(rx_value) {
    g <- baseline[rep(seq_len(n.obs), each = length(ends)), ]
    g$dtime <- rep(ends, n.obs)
    g$rx    <- rx_value
    g <- add_time_terms(g)
    g$hazardP <- predict(fitP, newdata = g, type = "response")
    g$hazardO <- predict(fitO, newdata = g, type = "response")
    g
  }
  treated <- make_grid(1L)
  placebo <- make_grid(0L)

  list(
    cifP_treated = calculateCumInc(treated),                    cifP_placebo = calculateCumInc(placebo),
    cifO_treated = calculateCumInc(treated, competing = TRUE),  cifO_placebo = calculateCumInc(placebo, competing = TRUE),
    weights = data.frame(patno = dt$patno, rx = dt$rx, eventType_int = dt$eventType_int, died = dt$allCause,
                         wE = dt$wE, wP = dt$wP, mx = dt$mx, Timey = dt$Timey)
  )
}

# The four estimands at 50 months: total effect on cancer death (RD, log RR) and effect on the competing event (RD, log RR)
gformula_est <- function(dt, method, mx_mult = 1) {
  f <- gformula_fit(dt, method, mx_mult)
  k <- length(ends)
  c(RD         = f$cifP_treated[k] - f$cifP_placebo[k],
    logRR      = log(f$cifP_treated[k] / f$cifP_placebo[k]),
    RD_comp    = f$cifO_treated[k] - f$cifO_placebo[k],
    logRR_comp = log(f$cifO_treated[k] / f$cifO_placebo[k]))
}

# Non-parametric benchmark: Aalen-Johansen cumulative incidence (randomised arms)
aj_fit <- function(dt) cmprsk::cuminc(ftime = pmax(dt$Time50, 0.5), fstatus = dt$eventType_int, group = dt$rx)
aj_est <- function(dt) {
  tp <- cmprsk::timepoints(aj_fit(dt), times = 50)$est
  g  <- function(arm, cause) tp[paste(arm, cause), 1]
  c(RD         = g(1, 1) - g(0, 1),
    logRR      = log(g(1, 1) / g(0, 1)),
    RD_comp    = g(1, 2) - g(0, 2),
    logRR_comp = log(g(1, 2) / g(0, 2)))
}
