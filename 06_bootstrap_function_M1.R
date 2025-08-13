# 06_bootstrap_function.R

# This script defines the boot_function used for bootstrap resampling
# and estimation using a parametric g-formula approach with competing risks.

boot_function_M1 <- function(data, indices){

  dt <- data[indices,]

  #p$tick()$print()  # update progress bar
  
  # Add patient number to the data
  dt$patno <- 1:nrow(dt)
  
  # Censor the patients at 2 years
  dt$cause <- ifelse(dt$Timey >= 2, 0, dt$cause)
  dt$vstatus <- ifelse(dt$Timey >= 2, 0, dt$vstatus)
  dt$Timey <- ifelse(dt$Timey >=2, 2, dt$Timey)
  
  # Summmarise the weights for cancer death
  SumW <- table(dt$cause)[2]
  NY <- table(dt$cause)[2]
  
  # Summmarise the weights for other cause of death
  SumP <- table(dt$cause)[3]
  ND <- table(dt$cause)[3]
  
  # Summmarise the weights for cancer death stratified by deprivation levels
  SumWND <- table(dt$cause[dt$dep==0])[2]
  NYND <- table(dt$cause[dt$dep==0])[2]
  SumWD <- table(dt$cause[dt$dep==1])[2]
  NYD <- table(dt$cause[dt$dep==1])[2]
  
  # Summmarise the weights for other cause of death
  SumPND <- table(dt$cause[dt$dep==0])[3]
  NDND <- table(dt$cause[dt$dep==0])[3]
  SumPD <- table(dt$cause[dt$dep==1])[3]
  NDD <- table(dt$cause[dt$dep==1])[3]
  
  # Create a binary variable to indicate all-cause death
  dt$allCause <- dt$vstatus != 0
  dt$eventType <- as.factor(dt$cause)
  dt$eventCens <- dt$eventType == 0
  dt$eventCancer <- dt$eventType == 1
  
  dt$Tstart = -0.01
  dt$dtimex <- dt$Timey*12
  dt$dtime <- ceiling(dt$dtimex)
  
  longdata <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="allCause")  
  longdataCens <- survSplit(data=dt, cut=cutTimes, start="Tstart", end="dtime", event="eventCens")
  longdata$eventCens <- longdataCens$eventCens
  
  longdata$CancerDeath <- longdata$allCause==1 & longdata$eventType==1
  longdata$otherDeath <- longdata$allCause==1 & longdata$eventType==2
  
  longdata$CancerDeath[longdata$eventCens==1]<-NA
  longdata$otherDeath[longdata$eventCens==1]<-NA
  longdata$CancerDeath[longdata$otherDeath==1]<-NA
  
  longdata <- longdata[longdata$dtime<length(cutTimes),]
  baseline <- longdata[longdata$dtime==0,]
  n.obs <- length(unique(longdata$patno))
  
  plrFitP <- glm(CancerDeath ~ dtime + agecr + IsexH + trt + dep, 
                 data = longdata, family=binomial())
  plrFitO <- glm(otherDeath ~ dtime + agecr + IsexH + trt + dep, 
                 data = longdata, family=binomial())
  
  treated <- baseline[rep(1:n.obs, each=length(cutTimes)),]
  treated$dtime <- rep(cutTimes, n.obs)
  treated$trt <- 1
  treated$hazardP <- predict(plrFitP, newdata = treated, type = 'response')
  treated$hazardO <- predict(plrFitO, newdata = treated, type = 'response')
  treated$s <- (1-treated$hazardP) * (1-treated$hazardO)
  
  placebo <- baseline[rep(1:n.obs, each=length(cutTimes)),]
  placebo$dtime <- rep(cutTimes, n.obs)
  placebo$trt <- 0
  placebo$hazardP <- predict(plrFitP, newdata = placebo, type = 'response')
  placebo$hazardO <- predict(plrFitO, newdata = placebo, type = 'response')
  placebo$s <- (1-placebo$hazardP) * (1-placebo$hazardO)
  
  cumIncTreatedRS <- calculateCumInc(treated)
  cumIncPlaceboRS <- calculateCumInc(placebo)
  log.gcomptotrr <- log(cumIncTreatedRS[length(cutTimes)] / cumIncPlaceboRS[length(cutTimes)])
  gcomptotrd <- cumIncTreatedRS[length(cutTimes)] - cumIncPlaceboRS[length(cutTimes)]
  
  cumIncTreated5comp <- calculateCumInc(treated, competing = TRUE)
  cumIncPlacebo5comp <- calculateCumInc(placebo, competing = TRUE)
  log.gcompcomprr <- log(cumIncTreated5comp[length(cutTimes)] / cumIncPlacebo5comp[length(cutTimes)])
  gcompcomprd <- cumIncTreated5comp[length(cutTimes)] - cumIncPlacebo5comp[length(cutTimes)]
  
  return(c(gcomptotrd, log.gcomptotrr, gcompcomprd, log.gcompcomprr))
}