#############################################################################################
# Title: On Causal Inference for the Relative Survival Setting
# Author: Matthew J. Smith
# Date: 12th February 2023
#############################################################################################


###################
# Acknowledgements:
###################

# 1) This code is based on the code written by Mats Stensrud et al in the paper:  
   # Jessica G. Young, Mats J. Stensrud, Eric J. Tchetgen Tchetgen, Miguel A. Hernán  (2020):
     # "A causal framework for classical statistical estimands in failure-time settings with 
     # competing events". https://doi.org/10.1002/sim.8471
   # Our code allows estimation when the cause of death (i.e., the outcome of interest) is 
     # not known.

# 2) The functions to simulate the cancer data set (i.e., patient characteristics and 
   # survival times were written and provided by Dr. Aurelien Belot. 

  
##########
# Outline:
##########

# This code shows how to estimate the total causal effect of an exposure on an
# outcome when the cause of death is unknown (i.e., relative survival setting).

# Exposure   : cmb (Binary: cmb = 0, cmb = 1) 
# Outcome    : Cancer-related death (death due to cancer)
# Confounders: age (15-99), sex (0 = male, 1 = female), L1 (0 = no, 1 = yes)



#############
# Set-up
#############

# Install necessary packages
  library(progress)
  library(dplyr)
  library(MASS)
  library(Epi)
  library(plyr)
  library(dplyr)
  library(mexhaz)
  library(utils)
  library(boot)

# Set working directory
  setwd("~/CI&RS")

# IMPORTANT!
  # Ensure you have the functions from the GitHub repo: mattyjsmith/CI-RS
      # 1. calculateCumInc
      # 2. cDataDesignOptim
      # 3. cdatasimulationT1WeibOptim
  # If you have the "Functions.R" script in your working directory then run the following:
  source("Functions.R")

# Ensure you have loaded the life table (from the GitHub repo):
  exprates <- readRDS(file="exprates.Rds")
    # The expected mortality rates are stratified by age, sex, year, and the exposure (cmb)


###############
# Preliminaries
###############
  
# Set the seed
  seed <- 1
  
# Set number of observations you want in your study
  n <- 10000
  
# Simulate the patient characteristics using the 'cDataDesignOptim' function
  # betaagecr: age (centred and rescaled) increases the chances of having the exposure by 2.5
  # betasex  : females have 0.8 times the odds of having the exposure compared to males
  # betaL1   : Those with L1 = 1 have 1.3 times the probability of having the exposure
  Simoptim <- cDataDesignOptim(n = n, seed=seed, cens.admin = 2015, 
                                ydiagmin=2005, ydiagmax=2010, 
                                betaagecr = 2.5, betasex = 0.8, betaL1 = 1.3)

  
# Simulate the survival time for each patient
  # Log-excess hazard scale (e.g., BetaSex=-0.07, so females have an excess hazard of 0.93
  # times the excess hazard of males ... a protective effect)
  Sigoptim <- cdatasimulationT1WeibOptim(mydat=Simoptim, seed=seed, myexpmort=exprates, 
                                lambda.weib=0.1, rho.weib=0.9, 
                                BetaAge = 0.5, Betasex = -0.07, BetaCmb = 0.3, BetaL1 = 0.2)
  
# State the number of time points (K) 
  # This is the number of cut times for the Pooled Logistic Regression
  cutTimes <- c(0:59)     
  # Here, K=60 because we are interested in survival at 5 years. 
  # The follow-up time (i.e., duee to administrative censoring) can be more than 5 years. 
  # The maximum follow-up time is 11 years (i.e., 2005 to end of 2015) 


# Rename your data set  
  data <- Sigoptim

# Add patient ID
  data$patno <- 1:n
  
  
########################################################
# Estimate the weights for the probability of death-type
########################################################

# Estimate the excess hazard
  invisible(capture.output(ehmodel <- mexhaz(Surv(finalsurvtime, vstatus) ~ agecr + IsexH + cmb + L1, 
                                             data = data, base = "weibull", 
                                             expected = "exprate")))
  # Check that the coefficients are the same as those that were specified when simulating the data
  summary(ehmodel) 
  # The log-excess hazard coefficients should be close to the values used when simulating the data
  
# Predict the excess hazard at each interval for each patient (10 seconds with 10,000 observations)
  for (k in 1:nrow(data)) {
    pred.model1 <- predict(ehmodel, 
                           time.pts = data$finalsurvtime[k], 
                           data.val = data.frame(agecr = data$agecr[k],
                                                 IsexH = data$IsexH[k],
                                                 cmb = data$cmb[k],
                                                 L1 = data$L1[k]))
    data$excesshazard[k] <- pred.model1[["results"]][["hazard"]] 
  }
  # It is important that the prediction of the excess hazard is from a model that has a longer
  # follow-up time than the survival time of interest. In other words, the model estimates the
  # excess hazard over 11 years of follow-up, but the survival of interest is 5-years.

# Censor the patients at 5 years (i.e., administrative censoring)
  data$cause <- ifelse(data$finalsurvtime >= 5, 0, data$cause)
  data$vstatus <- ifelse(data$finalsurvtime >= 5, 0, data$vstatus)
  data$finalsurvtime <- ifelse(data$finalsurvtime >=5, 5, data$finalsurvtime)
  
# Calculate weights for cancer-related death (w_{E})
  data$weightsE <- (data$excesshazard)/(data$excesshazard + data$exprate)
  
# Calculate weights for other-cause death (w_{P})
  data$weightsP <- (data$exprate)/(data$excesshazard + data$exprate)
  
# Summmarise the weights for cancer death
  SumW <- sum(data$weightsE[data$vstatus==1]) ; SumW
  NY <- sum(data$cause==1) ; NY       # Sanity check
  
# Summmarise the weights for other-cause death
  SumP <- sum(data$weightsP[data$vstatus==1]) ; SumP
  ND <- sum(data$cause==2) ; ND       # Sanity check
  
# Check the proportion of cancer deaths amongst all-cause death
  Prop.cancerdeaths <- SumW/sum(data$vstatus); Prop.cancerdeaths
  # The proportion of cancer deaths should be between 30-90%   
  # Outside of this range might not provide reliable weights.
  

  
#############################
# Convert to long-format data
#############################
  
# Create a binary variable to indicate all-cause death
  data$allCause <- data$vstatus != 0
  
# Create a variable with 3 levels indicating the cause of death (0 = cens, 1 = cancer, 2 = other)
  data$eventType <- as.factor(data$cause)
  
# Indicator for whether the patient was censored
  data$eventCens <- data$eventType == 0

# Set the start time of the study 
  data$Tstart = -0.01 # Needed for the long format
  
# Check that your survival time variable is the same format as the cutTimes variable
  data$dtimex <- data$finalsurvtime*12
  data$dtime <- ceiling(data$dtimex)
  
# Split each record into multiple sub-records at each cut time.
  # Here, each patient contributes to each cut time that they are observed within. Once they
  # experience an all-cause event, they do not have any further records.
  # "dtime" is the time in months that the patient was followed-up for.
  longdata <- survSplit(data=data, cut=cutTimes, start="Tstart", end="dtime", event="allCause")  
  
  # Do the same for the censored data
  longdataCens <- survSplit(data=data, cut=cutTimes, start="Tstart", end="dtime", event="eventCens")
  
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
  
# Parametric g-formula estimate of the total effect
  # Risk ratio
  gcomptotrr<-cumIncTreatedRS[length(cutTimes)]/cumIncPlaceboRS[length(cutTimes)]; print(gcomptotrr)       
  # RR = 1.174025 ... risk of cancer death among those with A=1 is 1.17 times that of those with A=0
  
  # Risk difference
  gcomptotrd<-cumIncTreatedRS[length(cutTimes)]-cumIncPlaceboRS[length(cutTimes)]; print(gcomptotrd)       
  # RD = 0.06572498 ... risk of cancer death among those with A=1 is 6.6% higher than those with A=0
  


############################################################################################################
# Construct confidence intervals - bootstrap the standard error
############################################################################################################
  
# To construct the confidence intervals you will need to get the Bootstrap function from the GitHub repo
  # Notice that the Bootstrap function is specifically designed for this analysis. If you want to adapt 
  # this analysis, you will need to adapt the Bootstrap function accordingly.
  source("Bootstrap function.R")
  
# Set the seed
  set.seed(1)
  
# Set number of bootstrap samples
  bootrep <- 300
  
# Specify the progress 
  p <- progress_estimated(bootrep + 1)
  
# Perform the bootstrap
  # "Estimand" option specifies whether you want Risk Difference (RD) or Relative Risk (RR). RD is the default.
  bsreps <- boot(data, boot_function, R=bootrep, estimand="RD")
  
# Mean of bootstrap replicates
  BSestimate <- mean(bsreps$t); BSestimate
  
# Confidence interval
  boot.ci1 <- boot.ci(bsreps, type=c('norm','basic','perc')); boot.ci1



##################################
# Create plots
##################################  
  
# Plot cumulative incidence
  plot(cutTimes,cumIncTreatedRS, type="s",ylim=c(0,1), ylab="Risk of cancer death", 
       xlab="Month",xlim=c(0,59), cex.lab=1.2, cex.main=1.2, lwd=1, lty=1, #xaxt="n",
       main ="") # Death due to Prostate Cancer: Estimand (5)
  lines(cutTimes, cumIncPlaceboRS, type="s",col=2,ylim=c(0,1),lwd=1,lty=1)
  #axis(1, at = seq(0, 59, by = 6))
  legend("topleft", c("A = 1", "A = 0"),
         col=c(1,2), lty=c(1,1),cex=0.9,pt.cex=1,lwd=2,bty="n")
