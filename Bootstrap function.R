# Define the function to obtain the statistics (i.e., risk difference or relative risk)
boot_function <- function(data, indices, estimand="RD"){
  dt <- data[indices,]
  
  p$tick()$print()  # update progress bar
  
  # Add patient ID
  dt$patno <- 1:nrow(dt)
  
  
  ########################################################
  # Estimate the weights for the probability of death-type
  ########################################################
  
  # Estimate the excess hazard
  invisible(capture.output(ehmodel <- mexhaz(Surv(finalsurvtime, vstatus) ~ agecr + IsexH + cmb + L1, 
                                             data = dt, base = "weibull", 
                                             expected = "exprate")))
  # Check that the coefficients are the same as those that were specified when simulating the data
  summary(ehmodel) 
  # The log-excess hazard coefficients should be close to the values used when simulating the data
  
  # Predict the excess hazard at each interval for each patient (10 seconds with 10,000 observations)
  for (j in 1:nrow(dt)) {
    pred.model1 <- predict(ehmodel, 
                           time.pts = dt$finalsurvtime[j], 
                           data.val = data.frame(agecr = dt$agecr[j],
                                                 IsexH = dt$IsexH[j],
                                                 cmb = dt$cmb[j],
                                                 L1 = dt$L1[j]))
    dt$excesshazard[j] <- pred.model1[["results"]][["hazard"]] 
  }
  # It is important that the prediction of the excess hazard is from a model that has a longer
  # follow-up time than the survival time of interest. In other words, the model estimates the
  # excess hazard over 11 years of follow-up, but the survival of interest is 5-years.
  
  # Censor the patients at 5 years (i.e., administrative censoring)
  dt$cause <- ifelse(dt$finalsurvtime >= 5, 0, dt$cause)
  dt$vstatus <- ifelse(dt$finalsurvtime >= 5, 0, dt$vstatus)
  dt$finalsurvtime <- ifelse(dt$finalsurvtime >=5, 5, dt$finalsurvtime)
  
  # Calculate weights for cancer-related death (w_{E})
  dt$weightsE <- (dt$excesshazard)/(dt$excesshazard + dt$exprate)
  
  # Calculate weights for other-cause death (w_{P})
  dt$weightsP <- (dt$exprate)/(dt$excesshazard + dt$exprate)
  
  # Summmarise the weights for cancer death
  SumW <- sum(dt$weightsE[dt$vstatus==1]) ; SumW
  NY <- sum(dt$cause==1) ; NY       # Sanity check
  
  # Summmarise the weights for other-cause death
  SumP <- sum(dt$weightsP[dt$vstatus==1]) ; SumP
  ND <- sum(dt$cause==2) ; ND       # Sanity check
  
  # Check the proportion of cancer deaths amongst all-cause death
  Prop.cancerdeaths <- SumW/sum(dt$vstatus); Prop.cancerdeaths
  # The proportion of cancer deaths should be between 30-90%   
  # Outside of this range might not provide reliable weights.
  
  
  
  #############################
  # Convert to long-format data
  #############################
  
  # Create a binary variable to indicate all-cause death
  dt$allCause <- dt$vstatus != 0
  
  # Create a variable with 3 levels indicating the cause of death (0 = cens, 1 = cancer, 2 = other)
  dt$eventType <- as.factor(dt$cause)
  
  # Indicator for whether the patient was censored
  dt$eventCens <- dt$eventType == 0
  
  # Set the start time of the study 
  data$Tstart = -0.01 # Needed for the long format
  
  # Check that your survival time variable is the same format as the cutTimes variable
  dt$dtimex <- dt$finalsurvtime*12
  dt$dtime <- ceiling(dt$dtimex)
  
  # Split each record into multiple sub-records at each cut time.
  # Here, each patient contributes to each cut time that they are observed within. Once they
  # experience an all-cause event, they do not have any further records.
  # "dtime" is the time in months that the patient was followed-up for.
  longdata <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="allCause")  
  
  # Do the same for the censored data
  longdataCens <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="eventCens")
  
  # Put the information of censoring time into the data for the all-cause times.
  longdata$eventCens <- longdataCens$eventCens
  
  # Restrict input data set to records with dtime<K+1 for fitting pooled over time models
  # In other words, stop the study at K by removing those records that are after K
  longdata <- longdata[longdata$dtime<length(cutTimes),]
  
  # Check the maximum of dtime should be less than K (since K(1) = 0)
  summary(longdata$dtime)
  
  # Create the baseline data (i.e., confounders measured at baseline)
  baseline <- longdata[longdata$dtime==0,]
  
  # Number of subjects
  n.obs <- length(unique(longdata$patno))
  
  
  
  
  ############################################################################################
  # Fit weighted conditional pooled logistic regression (WPLR) models later used for g-formula
  ############################################################################################
  
  # Create variables for the weights at each event time
  # The variable will be 1.0 when the event does not occur (i.e., no death)
  # The variable will be the weight when an event does occur (i.e., any death)
  method1 <- longdata
  
  # Weights for the cancer-related death
  method1$weightsC <- ifelse(method1$allCause==1, method1$weightsE, 1)
  
  # Weights for the other-cause death
  method1$weightsO <- ifelse(method1$allCause==1, method1$weightsP, 1)
  
  # WPLR model for the cancer-related death
  suppressWarnings({
    invisible(capture.output(plrFitP1 <- glm(allCause ~ dtime + agecr + IsexH + cmb + L1, 
                                             data = method1, family=binomial(), weight = weightsC)))})
  # Note that this model is weighted by the probability of cancer-related death (i.e., "weightsC")  
  
  # WPLR model for the other-cause death
  suppressWarnings({
    invisible(capture.output(plrFitO1 <- glm(allCause   ~ dtime + agecr + IsexH + cmb + L1, 
                                             data = method1, family=binomial(), weight = weightsO)))})
  # Note that this model is weighted by the probability of other-cause death (i.e., "weightsO")
  
  
  ####################################################################
  # Create expanded data sets to obtain parametric g-formula estimates 
  ####################################################################
  
  # Expand baseline so it contains a visit at each time point for every individual
  # where the baseline information has been carried forward at each time
  treated <- baseline[rep(1:n.obs,each=length(cutTimes)),] #One row for each interval k+1 for k=0,...,K
  treated$dtime <- rep(cutTimes,n.obs)
  treated$cmb <-1      # Assign all patient to have had the treatment
  
  # "Predict" the conditional discrete hazards (cause-specific for each cause of death/competing event 
  # conditioned for event and hazard of competing event) for each subject in each time interval
  treated$hazardP <- predict(plrFitP1, newdata = treated, type = 'response')       # Hazard from the outcome model
  treated$hazardO <- predict(plrFitO1, newdata = treated, type = 'response')       # Hazard from the competing event model
  treated$s <- (1-treated$hazardP) * (1-treated$hazardO)                      # Product of the hazards prior to time K
  # sum(treated$hazardO < 0)           # sanity check  
  
  
  # Make analogous dataset for placebo (do the same for when the exposure = 0)
  placebo <- baseline[rep(1:n.obs,each=length(cutTimes)),] #One row for each time
  placebo$dtime <- rep(cutTimes,n.obs)
  placebo$cmb <- 0 
  
  # Estimate conditional discrete hazard for each subject in each time interval
  placebo$hazardP <- predict(plrFitP1, newdata = placebo, type = 'response') 
  placebo$hazardO <- predict(plrFitO1, newdata = placebo, type = 'response')
  placebo$s <- (1-placebo$hazardP) * (1-placebo$hazardO)  
  
  
  ############################################################################################################
  # Calculate risks and corresponding treatment effects (RR/RD) by end of follow-up using parametric g-formula 
  ############################################################################################################
  
  # IMPORTANT!
  # Ensure you have the function "calculateCumInc" from the R script "Functions.R" from the GitHub repo: mattyjsmith/CI-RS
  
  # Total Causal Effect of X -> Y
  
  # Estimate the risk of cancer-related death at each event time for each level of the exposure (A=1 and A=0)
  cumIncTreatedRS <- calculateCumInc(treated)#a=1
  cumIncPlaceboRS <- calculateCumInc(placebo)#a=0
  
  if (estimand=="RR") {
    # Parametric g-formula estimate of the total effect
    # Risk ratio
    gcomptotrr<-cumIncTreatedRS[length(cutTimes)]/cumIncPlaceboRS[length(cutTimes)]  
    # RR = 1.174025 ... risk of cancer death among those with A=1 is 1.17 times that of those with A=0  
  } else {
    # Risk difference
    gcomptotrd<-cumIncTreatedRS[length(cutTimes)]-cumIncPlaceboRS[length(cutTimes)]    
    # RD = 0.06572498 ... risk of cancer death among those with A=1 is 6.6% higher than those with A=0
  }
  
}