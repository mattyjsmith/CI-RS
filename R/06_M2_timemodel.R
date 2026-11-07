# GENERATED from the Paper folder's function; the only change is the time model of the logistic regressions (argument time_model).
# CORRECTED (see README_paper.md): indexing, dummy row, other-cause model.
# 06_bootstrap_function.R

# This script defines the boot_function used for bootstrap resampling
# and estimation using a parametric g-formula approach with competing risks.

# Define the Bootstrap function
boot_function_M2_tm <- function(data, indices, time_model = "linear"){
  dt <- data[indices,]
  
  #p$tick()$print()  # update progress bar
  
  # Add patient number to the data
  dt$patno <- 1:nrow(dt)
  # Time model of the pooled logistic regressions: linear (as in the main simulation), interval-specific
  # intercepts (one per month), or a natural cubic spline with 4 df
  time_term <- switch(time_model,
                      linear = "dtime",
                      factor = "factor(dtime)",
                      spline = "ns(dtime, df = 4)")
  ends <- round(cutTimes + (cutTimes[2] - cutTimes[1]), 10)   # end time of each real (monthly) interval
  
  # Obtain the excess hazard Data
  # the excess hazard step is the slow one: it is done once per data set and passed in as a column
  if (is.null(dt$excesshazard)) {
  invisible(capture.output(ehmodel <- mexhaz(Surv(Timey, vstatus) ~ agecr + IsexH + dep + trt, 
                                             data = dt, base = "weibull", 
                                             expected = "exprate")))
  
  # Specify the excess hazard at the time of each event
  # for (k in 1:nrow(dt)) {
  #   pred.model1 <- predict(ehmodel, 
  #                          time.pts = dt$Timey[k], 
  #                          data.val = data.frame(agecr = dt$agecr[k],
  #                                                IsexH = dt$IsexH[k],
  #                                                dep = dt$dep[k],
  #                                                trt = dt$trt[k]))
  #   dt$excesshazard[k] <- pred.model1[["results"]][["hazard"]]
  # }
  
  dt$excesshazard <- vapply(
    seq_len(nrow(dt)),
    function(k) {
      predict(ehmodel,
              time.pts = dt$Timey[k],
              data.val = dt[k, , drop = FALSE]
      )$results$hazard
    },
    numeric(1)
  )
  }
  
  # Censor the patients at 2 years
  dt$cause <- ifelse(dt$Timey >= 2, 0, dt$cause)
  dt$vstatus <- ifelse(dt$Timey >= 2, 0, dt$vstatus)
  dt$Timey <- ifelse(dt$Timey >=2, 2, dt$Timey)
  
  # Calculate weights of cancer death
  dt$weightsE <- (dt$excesshazard)/(dt$excesshazard + dt$exprate)
  
  # Calculate weights of other cause of death
  dt$weightsP <- (dt$exprate)/(dt$excesshazard + dt$exprate)
  
  # Summmarise the weights for cancer death
  SumW <- sum(dt$weightsE[dt$vstatus==1]) ; SumW
  NY <- sum(dt$cause==1) ; NY       # Sanity check
  
  # Summmarise the weights for other cause of death
  SumP <- sum(dt$weightsP[dt$vstatus==1]) ; SumP
  ND <- sum(dt$cause==2) ; ND       # Sanity check
  
  
  #################
  # Data management
  #################
  
  # Create a binary variable to indicate all-cause death
  dt$allCause <- dt$vstatus != 0
  
  # Create a variable with 3 levels indicating the cause of death 
  dt$eventType <- as.factor(dt$cause)
  #as.factor(dt$cause)
  #table(dt$eventType)
  
  # Tabulate the type of event by the exposure level
  #table(dt$eventType, dt$cmb)
  
  # Data processing continues: Comorbidity (A=1) and no comorbidity (A=0)
  
  # Indicator for whether the patient was censored
  dt$eventCens <- dt$eventType == 0
  #table(dt$eventCens)
  
  # Indicator for whether the patient died due to cancer
  #data$eventProst <- data$eventType == 1
  
  
  # Convert data to long format 
  
  # Set the start time of the study 
  dt$Tstart = -0.01 #Needed for the long format
  
  # Check that your survival time variable is the same format as the cutTimes variable
  #summary(dt$finalsurvtime)
  dt$dtimex <- dt$Timey*12
  dt$dtime <- ceiling(dt$dtimex)
  
  # Split each record into multiple subrecords at each cut time.
  # Here, each patient contributes to each cut time that they are observed within. Once they
  # experience an all-cause event, they do not have any further records.
  # "dtime" is the time in months that the patient was followed-up for.
  longdata <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="allCause")  
  
  # We do the same for the censored data
  longdataCens <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="eventCens")
  
  # We put the information of censoring time into the data for the all-cause times.
  longdata$eventCens <- longdataCens$eventCens
  
  # Restrict input data set to records with dtime<K+1 for fitting pooled over time models
  # In other words, stop the study at K by removing those records that are after K
  baseline <- longdata[longdata$dtime==0,]   # one zero-length dummy row per person: baseline covariates only
  longdata <- longdata[longdata$dtime>0,]    # the models are fitted on real person-intervals only
  # To confirm the maximum of dtime should be less than K (since K(1) = 0)
  #summary(longdata$dtime)
  
  # Create 'baseline' - data collected at visit 0
  # Create the baseline data (i.e., confounders measured at baseline)
  
  # Number of subjects
  n.obs <- length(unique(longdata$patno))
  
  
  ############
  # Fit conditional pooled logistic regression (PLR) models later used for parametric g-formula and/or ipw
  # 
  # Note all models are pooled over all patients
  ############
  
  # 1) Weights should be 1 for all patient for all intervals, except for the last if they died (in which case their weight will be the actual weights)
  
  method1 <- longdata
  method1$weightsC <- ifelse(method1$allCause==1, method1$weightsE, 1)
  method1$weightsO <- ifelse(method1$allCause==1, method1$weightsP, 1)
  # CORRECTION: other-cause model = fractional outcome (expected other-cause death), weight 1 on every row, so the
  # hazard is (other-cause deaths)/(everyone at risk), as assumed by calculateCumInc
  method1$yO <- ifelse(method1$allCause==1, method1$weightsP, 0)
  
  # PLR for the death due to cancer
  suppressWarnings({
    invisible(capture.output(plrFitP1 <- glm(as.formula(paste("allCause ~", time_term, "+ agecr + IsexH + dep + trt")), 
                                             data = method1, family=binomial(), weight = weightsC)
    ))})
  
  # PLR for the death due to other causes (only need variables that would be in life table)
  suppressWarnings({
    invisible(capture.output(plrFitO1 <- glm(as.formula(paste("yO ~", time_term, "+ agecr + IsexH + dep + trt")),
                                             data = method1, family=quasibinomial())
    ))})
  
  
  # PLR for censoring
  #temp50<-longdata[longdata$dtime>50,]        # Select only those records for k>50
  #plrFitC <- glm(eventCens ~ agecr + IsexH + cmb, data = temp50,family=binomial()); summary(plrFitC)   # PLR model on those who have k>50 ???
  
  # Fit models for the numerators in the stabilized weights:
  #plrFitCnum <- glm(eventCens ~ cmb , data = temp50,family=binomial()); summary(plrFitCnum)
  #plrFitOnum <- glm(otherDeath ~ dtime + agecr + IsexH + cmb, data = longdata,family=binomial()); summary(plrFitOnum)
  
  
  
  # 2) outcome is Y = "weightsE * allCause", and remove weighting option
  
  # method2 <- longdata
  # method2$eventY <- method2$allCause * method2$weightsE
  # method2$eventD <- method2$allCause * method2$weightsP
  # 
  # # PLR for the death due to cancer
  # plrFitP2 <- glm(eventY ~ dtime + agecr + IsexH + cmb,
  #                 data = method2, family=binomial()); summary(plrFitP2)
  # 
  # # PLR for the death due to other causes
  # plrFitO2 <- glm(eventD ~ dtime + agecr + IsexH + cmb,
  #                 data = method2, family=binomial()); summary(plrFitO2)
  
  ############
  # Create expanded data sets to obtain parametric g-formula estimates 
  ############
  
  # Create expanded data so that hazards given A=1 can be computed for all subjects (even for those 
  # with A=0) -- note in an RCT this is not necessary, we could stratify models above and this step by treatment arm
  
  # Expand baseline so it contains a visit at each time point for every individual
  # where the baseline information has been carried forward at each time
  treated <- baseline[rep(1:n.obs,each=length(ends)),] #One row for each interval k+1 for k=0,...,K
  treated$dtime <- rep(ends,n.obs)
  treated$trt <-1      # Assign all patient to have had the treatment
  
  
  # "Predict" the conditional discrete hazards (cause specific for each cause of death/competing event 
  # conditioned for event and hazard of competing event) for each subject in each time interval
  treated$hazardP <- predict(plrFitP1, newdata = treated, type = 'response')       # Hazard from the outcome model
  treated$hazardO <- predict(plrFitO1, newdata = treated, type = 'response')       # Hazard from the competing event model
  treated$s <- (1-treated$hazardP) * (1-treated$hazardO)                      # Product of the hazards prior to 
  # sum(treated$hazardO < 0)           # sanity check  
  
  
  # Make analogous dataset for placebo
  placebo <- baseline[rep(1:n.obs,each=length(ends)),] #One row for each time
  placebo$dtime <- rep(ends,n.obs)
  placebo$trt <- 0 
  
  # Estimate conditional discrete hazard for each subject in each time interval
  placebo$hazardP <- predict(plrFitP1, newdata = placebo, type = 'response') 
  placebo$hazardO <- predict(plrFitO1, newdata = placebo, type = 'response')
  placebo$s <- (1-placebo$hazardP) * (1-placebo$hazardO)
  
  # Parametric g-formula estimate of risk without elimination of competing events
  cumIncTreatedRS <- calculateCumInc(treated)#a=1
  cumIncPlaceboRS <- calculateCumInc(placebo)#a=0
  
  # Parametric g-formula estimate of the total effect
  # Risk ratio
  log.gcomptotrr <- log(cumIncTreatedRS[length(ends)]/cumIncPlaceboRS[length(ends)]) #JGY need ci's on this, I think let's skip composite outcome this is all enough
  #logRR <- log(gcomptotrr)
  # print(gcomptotrr)       # RR = 1.264382 ... risk of cancer death among those with comorbidity is 1.27 times that of those without comorbidity
  
  # Risk difference
  gcomptotrd<-cumIncTreatedRS[length(ends)]-cumIncPlaceboRS[length(ends)] #JGY need ci's on this
  # print(gcomptotrd)       # RD = 0.08825672 ... risk of cancer death among those with comorbidity is 8.8% higher than those without comorbidity
  
  
  # Parametric g-formula estimate of the risk of the competing event itself
  cumIncTreated5comp <- calculateCumInc(treated, competing = TRUE)#a=1
  cumIncPlacebo5comp <- calculateCumInc(placebo, competing = TRUE)#a=0
  #(Total) effect on risk of competing event
  log.gcompcomprr <- log(cumIncTreated5comp[length(ends)]/cumIncPlacebo5comp[length(ends)]) 
  gcompcomprd <- cumIncTreated5comp[length(ends)]-cumIncPlacebo5comp[length(ends)] 
  #print(log.gcompcomprr)
  #print(gcompcomprd)
  
  return(c( gcomptotrd, log.gcomptotrr, gcompcomprd, log.gcompcomprr))
  
}