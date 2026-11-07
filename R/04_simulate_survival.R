# 04_simulate_survival.R
# CORRECTED VERSION (other-cause time over full follow-up; counterfactual worlds) (see README_paper.md): same scenarios and methods as the original,
# with the corrections found in October 2026 applied.
# Full version of cdatasimulationT1WeibOptim.v2 with support for:
# - Weibull or Generalised Weibull survival
# - Non-proportional hazards (NPH) by sex
# - Expected mortality via Lexis split
# - Competing risks (cancer vs other-cause death)

# Function to simulate the survival times
cdatasimulationT1WeibOptim.v2 <- function (mydat, dist.T1="Weibull", lambda.weib=NULL, rho.weib=NULL, 
                                           k1=NULL, r1=NULL, a1=NULL, myexpmort,
                                           NPH=FALSE, lambda.weib2=NULL, rho.weib2=NULL, 
                                           BetaAge=0, Betasex=0, Betatrt=0, Betadep=0, Betastage=0)
{
  
  #set.seed(seed)
  
  # **************************************************************************** #
  # mydat : Corresponds to the list returned by cDataDesignOptim to collect
  #          * tab2 = Design matrix (NEEDS TO BE ORDERED BY Ident)
  #          * and general information : n, NCluster, cens.admin, ydiagmin, ydiagmax
  # dist.T1                : Distribution of T1 "Weibull" or "GenWeibull"
  # lambda.weib, rho.weib  : scale & shape parameters for the Weibull baseline hazard
  # k1, r1, a1             : parameters for the Generalised Weibull Distribution
  # mu                     : mean of the multinormal random effect
  # Sigma                  : Variance-covariance matrix of the random Effect (so diagonal element= variances)
  # NPH                    : TRUE/FALSE Is there an NPH effect to simulate for IsexH?
  #        If true, lambda.weib2 & rho.weib2 needs to be given to define the baseline hazard for women
  #                while lambda.weib, rho.weib define the baseline hazard for men
  # lambda.weib2, rho.weib2 : scale & shape parameters for the Weibull baseline hazard (used only if NPH=TRUE)
  # myexpmort              : data.frame of the expected mortality hazard
  # lambda.weib2, rho.weib2 : scale & shape parameters for the Weibull baseline hazard (used only if NPH=TRUE)
  # BetaAge      : Age effect
  #                if BetaAge == 0    => no age effect on the cancer mortality hazard (i.e. T1)
  # Betasex     : Gender effect (Ref=women)
  #                If Betasex == 0    => No sex effect on the cancer mortality hazard (i.e. T1)
  # Betatrt     : Treatment effect (Ref=placebo)
  #                If Betatrt == 0    => No treatment effect on the cancer mortality hazard (i.e. T1)
  # Betadep     : Deprivation effect (Built at the cluster level): Higher level of depriv means higher deprivation
  #                If Betadepriv == 0 => No deprivation effect on the cancer mortality hazard (i.e. T1)
  # **************************************************************************** #
  
  # Survival Time from cancer T1  (Weibull)
  # lambda.weib=0.25; rho.weib=0.7;
  # tim <- seq(0,10,by=0.01)
  # haz <- lambda.weib*rho.weib*tim^(rho.weib-1)
  # plot(tim, haz, type="l")
  
  # Survival Time from cancer T1  (Generalised Weibull)
  # k=2; r=0.5; a=0.3;
  # t=seq(0,10,by=0.01)
  # haz<-k*(r^k)*(t)^(k-1)/(1+((r*t)^k)/a)
  # plot(t, haz, type="l")
  
  # EXAMPLE
  # tata = cdatasimulationT1WeibOptim(mydat=toto, lambda.weib=0.4, rho.weib=0.4,
  #                                   Sigma=matrix(c(0.3,0.5*sqrt(0.3)*sqrt(0.5),0.5*sqrt(0.3)*sqrt(0.5),0.5),2,2),
  #                                   BetaAge = 0, Betasex=0, BetaTTT=0, Betadepriv=0, myexpmort=tauxatt)
  
  # **************************************************************************** #
  
  
  # Warning messages to make sure you have included the necessary arguments
  if (dist.T1=="Weibull" & (is.null(lambda.weib) | is.null(rho.weib))) {stop("Arguments lambda.weib and rho.weib should be given for Weibull simulation!")}

  ########
  # Step 1
  ########
  
  #print("STEP 1 ")
  #print(date())
  
  # The list you created in cDataDesignOptim forms the cancer dataset you will use
  tab2 <- mydat$tab2
  
  
  # Collect Information useful from the Design
  # Number of obs
  n <- mydat$n
  # Year of administrative censoring 
  cens.admin <- mydat$cens.admin
  # Earliest entry year in study
  ydiagmin <- mydat$ydiagmin
  # Latest entry year in study
  ydiagmax <- mydat$ydiagmax
  # Maximum potential time of follow up
  maxptime <- mydat$maxptime
  
  
  # The linear model 
  tab2$exp.xbeta <- exp(BetaAge * tab2$agecr + Betasex * tab2$IsexH +      # Age and sex
                          Betatrt * tab2$trt + Betadep * tab2$dep +
                          Betastage * tab2$stage
  ) 
  # exp(-xbeta)
  tab2$exp.Mxbeta   <- (tab2$exp.xbeta)^(-1)
  
  # Time from diagnosis to censoring
  tab2$potentime   <- tab2$cens.admin - tab2$year
  
  ########
  # Step 2
  ########
  
  #print("STEP 2 ")
  #print(date())
  
  # Survival Time from cancer T1  (Weibull)
  # Random number of n values

  tab2$Trt0 <- 0
  tab2$Trt1 <- 1
  
  # Counterfactual world: As if everyone was treated 
  temp.ut1  <- runif(nrow(tab2))
  
  tab2$exp.xbetaT1 <- exp(BetaAge * tab2$agecr + Betasex * tab2$IsexH +      # Age and sex
                                     Betatrt * tab2$Trt1 + Betadep * tab2$dep + Betastage * tab2$stage)
  tab2$TimeT1Cause1 <- exp((1/rho.weib) * (log(-1/lambda.weib*log(temp.ut1)) - log(tab2$exp.xbetaT1))) # cancer

  
  # Counterfactual world: As if everyone was not treated 
  temp.ut0  <- runif(nrow(tab2))
  
  tab2$exp.xbetaT0 <- exp(BetaAge * tab2$agecr + Betasex * tab2$IsexH +      # Age and sex
                                     Betatrt * tab2$Trt0 + Betadep * tab2$dep + Betastage * tab2$stage)
  
  tab2$TimeT0Cause1 <- exp((1/rho.weib) * (log(-1/lambda.weib*log(temp.ut0)) - log(tab2$exp.xbetaT0))) # cancer


  tab2$T1 <- ifelse(tab2$trt==1, tab2$TimeT1Cause1, tab2$TimeT0Cause1) # cancer
  
  # Survival Time from expected mortality T2
  # Minimum of time from cancer or time to censoring
  tab2$T1cens <- pmin(tab2$T1, tab2$potentime)
  # Flag those who died due to cancer (i.e., cancer time less than censoring time)
  tab2$cause1 <- ifelse(tab2$T1<tab2$potentime,1,0)
  # Age at censoring
  # CORRECTION: the other-cause time T2 is simulated over the FULL potential follow-up (to administrative
  # censoring), not only up to the observed cancer/censoring time. T2 is a person-level latent time shared by
  # both counterfactual worlds, so it must exist after the observed cancer death too (e.g. when the other
  # treatment arm would have lived longer).
  tab2$ageout <- tab2$age+tab2$potentime
  # Year at censoring
  tab2$yearout <- tab2$year+tab2$potentime
  
  # Create a Lexis object to split time along the two timescales (agediag and ydiag)
  # Object that represents follow-up in multiple states on multiple time scales (i.e., death from cancer, or censored, time to death or censor)
  trylx <- Lexis(entry       = list(yearscale=year, agescale=age), 
                 exit        = list(yearscale=yearout, agescale=ageout),
                 exit.status = cause1, 
                 data=tab2, tol=min(tab2$T1)-0.00001)
  
  # Split time along two time-axes
  # Divide each row of the Lexis object into disjoint follow-up intervals according to the supplied break points
  # Year-scale
  trylxsplit <- splitLexis( trylx, breaks = seq(ydiagmin,ydiagmin+maxptime,1), time.scale="yearscale" )
  # Age-scale
  trylxsplit <- splitLexis( trylxsplit, breaks = seq(ceiling(min(tab2$age)),ceiling(max(tab2$age)+maxptime),1), time.scale="agescale" )
  
  ###  tab3 <- rename(tab3, c("age.entier"="AGEX" ,"year.entier"="year","sex"="sex"))
  ###  tab4 <- join(tab3, tauxatt, by = c("AGEX","year","sex"), match="first")
  
  # Truncate the variable to an integer
  trylxsplit$AGERT <-trunc(trylxsplit$agescale)
  trylxsplit$YEARRT <- trunc(trylxsplit$yearscale)
  trylxsplit$SEXRT <- trylxsplit$sex
  trylxsplit$DEPRT <- trylxsplit$dep
  
  # trylxsplit <- rename(trylxsplit, c("age.entier"="AGEX" ,"year.entier"="year","sex"="sex"))
  
  # Join the data sets together
  tab5 <- join(data.frame(trylxsplit), myexpmort, by = c("AGERT","YEARRT","SEXRT", "DEPRT"), match="first")
  
  # Calculation of T2bis (Used to calculate the "final" T2)
  # Calculate the n values for the expected mortality rate with rate equal to 'exprate' in the life table data
  tab5$T2bis <- rexp(nrow(tab5), tab5$exprate)
  
  
  ########
  # Step 3
  ########
  
  #print("STEP 3 ")
  #print(date())
  
  # 
  # Extract those where T2bis is less than lex.dur for the variables of interest
  tempT2 <- tab5[tab5$T2bis<tab5$lex.dur, c("Ident", "yearscale", "agescale", "age", "T2bis")]
  # Order by the variables of interest 
  tempT2 <- tempT2[order(tempT2$Ident, tempT2$yearscale, tempT2$agescale, tempT2$T2bis),]
  # Merge the cancer data to the life table data
  tabtemp <- join(tab2, tempT2, by=c("Ident"), match="first")
  # Expected survival time is the age difference (if T2bis is false), otherwise use potential follow up time
  tabtemp$T2 <- ifelse(is.na(tabtemp$T2bis)==F, tabtemp$agescale+tabtemp$T2bis-tabtemp$age, tabtemp$potentime)
  # Flag those with other causes of death if T2bis observed
  tabtemp$cause2 <- ifelse(is.na(tabtemp$T2bis)==F, 1, 0)
  
  # Final Expected Survival Time: T2
  
  # Final Observed Survival Time: finalsurvtime
  
  #  Temps       <- data.frame(T1 = T1, T2 = tabtemp$T2)
  # Final survival time is the minimum of T1 (time to cancer death) or T2 (time to other cause death), whichever occurred first
  Timey  <- pmin(tabtemp$T1, tabtemp$T2)
  
  # Join the data sets of: cancer data, time to other cause death, flag of other cause death, and actual final survival time
  datasimulation <- data.frame(tab2, T2 = tabtemp$T2, cause2=tabtemp$cause2, Timey)
  

  datasimulation$TimeT1Cause2 <- datasimulation$T2 # other cause
  datasimulation$TimeT1 <- pmin(datasimulation$TimeT1Cause1, datasimulation$TimeT1Cause2)
  datasimulation$CauseT1 <- ( 0
                    + 1 * (datasimulation$TimeT1 == datasimulation$TimeT1Cause1 & datasimulation$TimeT1 < datasimulation$potentime )
                    + 2 * (datasimulation$TimeT1 == datasimulation$TimeT1Cause2 & datasimulation$TimeT1 < datasimulation$potentime ))
  
  datasimulation$TimeT0Cause2 <- datasimulation$T2 # other cause
  datasimulation$TimeT0 <- pmin(datasimulation$TimeT0Cause1, datasimulation$TimeT0Cause2)
  datasimulation$CauseT0 <- ( 0
                    + 1 * (datasimulation$TimeT0 == datasimulation$TimeT0Cause1 & datasimulation$TimeT0 < datasimulation$potentime )
                    + 2 * (datasimulation$TimeT0 == datasimulation$TimeT0Cause2 & datasimulation$TimeT0 < datasimulation$potentime ))
  
  
  ########
  # Step 4
  ########
  
  
  #print("STEP 4 ")
  #print(date())
  
  # Vital Status ("vstatus")  : 1 ==> death, 0 ==> alive
  
  # Create a flag variable that holds the vital status if final survival time was the actual survival time
  datasimulation[, "vstatus"] <- ifelse(datasimulation[, "Timey"] < datasimulation[, "potentime"], 1 , 0)
  
  # Cause of death : 1 ==> Cancer, 2 ==> Other causes
  
  # Create a flag variable for the cause (cancer or otherwise)
  datasimulation$cause <-( 0
                           + 1 * (datasimulation$Timey == datasimulation$T1 & datasimulation$vstat == 1)
                           + 2 * (datasimulation$Timey == datasimulation$T2 & datasimulation$vstat == 1) )
  
  
  # Update of the yearout and ageout to apply for finalsurvtime and not T1 only
  datasimulation$yearout <- datasimulation$year + datasimulation$Timey
  datasimulation$ageout <- datasimulation$age + datasimulation$Timey
  
  # ----- datasimulation : final data frame, with the life-table rates at the end of follow-up ----- #
  datasimulation$truncyearout <- trunc(datasimulation$yearout)
  datasimulation$truncageout <- trunc(datasimulation$ageout)
  datasimulation <- merge(datasimulation, myexpmort,
                          by.x=c("truncyearout","truncageout","sex", "dep"),
                          by.y=c("YEARRT","AGERT","SEXRT","DEPRT"))
  
  ########
  # Step 5
  ########
  
  #print("STEP 5 ")
  #print(date())
  
  # Return the data set that holds the variables
  return(datasimulation)
  
}
