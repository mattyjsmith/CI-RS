# User defined functions

# Function to simulate the patient characteristics
cDataDesignOptim <- function (n = 500, seed=1234, cens.admin = 2012, ydiagmin=1998, ydiagmax=2000,
                              NCluster = 10, ClusterSize = NULL,
                              betaagecr = 2.5, betasex = 0.8, betaL1 = 1.3) {
  set.seed(seed)
  
  # **************************************************************************** #
  # n					   : Sample size
  # cens.admin : Year of administrative censoring (if == 1995 : Maximum potential time of F-Up maxptime = 15. If 1991, maxptime = 11)
  #              Needs to be greater than 1990
  # ydiagmin, ydiagmax : Min and max of the period of diagnosis (parameter linked to cens.admin)
  # NCluster : Number of cluster per sample
  # ClusterSize : Vector of number of observations per cluster
  #              If NULL, then balanced cluster size given n and NCluster
  
  # EXAMPLE
  # toto <- cDataDesignOptim(n = 10000, NCluster = 100, cens.admin = 1995, ydiagmin=1980, ydiagmax=1981)
  # **************************************************************************** #
  
  # ------------------------------ Covariables --------------------------------- #
  
  ########
  # Step 1
  ########  
  print("STEP 1 ")
  print(date())
  
  # Maximum potential time of F-Up
  maxptime <- cens.admin - ydiagmin
  
  # ID patient
  Ident <- c(1:n)
  
  # Sex : Men ( sex = 1 ) - Women ( sex = 2 )
  sex <- ifelse(runif(n) < 0.5, 2, 1)
  IsexH <- ifelse(sex==2,1,0)
  
  # Treatment
  TTT <- ifelse(runif(n) < 0.5, 0, 1)
  
  # Year of diagnosis
  year <- runif(n, min = ydiagmin, max = ydiagmax)
  
  # Age at diagnosis (continuous)
  age   <- c(runif(n*0.25, min = 30, max = 65), runif(n*0.35, min = 65, max = 75), runif(n*0.40, min = 75, max = 85))
  agecr <- (age - mean(age))/10
  
  # Deprivation (continuous)
  depriv <- rnorm(NCluster,0,sd=1.5) # Value of deprivation is from a N~(0,1.5) taken from n number of clusters
  
  # Ethnicity
  ethnic <- rbinom(n,1,prob=0.90)   # Roughly 90% of observed ethnicity is white
  
  # Diagnostic route
  route <- rbinom(n,size=1,prob=0.30)   # 1 = emergency (about 30% A&E), 0 = other route
  
  # Create stage variable
  stage <- rbinom(n,size=1,prob=0.55) # 1 = late stage (about 55% late stage)
  
  # Unmeasured blood tests
  L1 <- rbinom(n,size=1,prob=0.55) # 1 = late stage (about 55% late stage)
  
  # Comorbidity
  comorbid <- rbinom(n,1,prob=0.40)  # Roughly 40% of patients have comorbidity
  
  # Create exposure variable
  z <- 1 + betaagecr*agecr + betasex*IsexH + betaL1*L1 # Linear combination
  pr <- 1/(1+exp(-z))         # Inv-logit function
  cmb <- rbinom(n, 1, prob = pr)  # Roughly 40% of patients have comorbidity
  
  
  # Create a data frame to hold the ID of the cluster and the deprivation value assigned to it
  DF.depriv<- data.frame(idcluster=seq(1:NCluster), depriv=depriv)
  
  
  # Assign patients to clusters
  if (is.null(ClusterSize)) {
    cluster <- sample(rep(1:NCluster,each=n/NCluster)) 
  } else {
    if (sum(ClusterSize)!=n)
      stop("Total Number of obs different from the total sum of clusters' size !")
    cluster <- sample(rep(1:length(ClusterSize),times=ClusterSize))
  }
  
  #hist(cluster)
  
  
  ########
  # Step 2
  ########
  
  print("STEP 2 ")
  print(date())
  
  # tab 2 : dataframe with 11 columns
  #        Ident         : Id patient
  #        sex           : Sex
  #        year          : Year of diagnosis
  #        year.entier   : Truncated Year of diagnosis
  #        age           : Age at diagnosis
  #        age.entier    : Truncated Age at diagnosis
  #        agec          : Age centered
  #        exp.xbeta     : Exponential of the linear predictor
  #        cens.admin    : Year of administrative censoring
  #        maxptime      : Maximum potential time of F-Up
  #        depriv        : Deprivation level (higher=more deprived)
  
  # Create a data frame to hold the variables you have created
  tab2              <- data.frame(Ident, sex, IsexH, TTT, year, age, cluster, ethnic, cmb, route, stage, L1)
  
  # Merge the cluster data frame to the cancer data frame
  tab2              <- merge(tab2, DF.depriv, by.x="cluster", by.y="idcluster")
  
  # Create centered and rescaled versions of the age variable
  tab2$agec         <- tab2$age - mean(tab2$age)
  tab2$agecr        <- tab2$agec/10
  tab2$agecr2       <- tab2$agecr^2
  tab2$agecr2plus   <- tab2$agecr2*(tab2$agecr>0)
  tab2$maxptime     <- rep(maxptime, n)
  tab2$cens.admin   <- rep(cens.admin, n)
  
  # Create categorical deprivation variable, and dummy variables  
  library(gtools)
  tab2$EDI <- as.numeric(quantcut(tab2$depriv, q=5))              # Create quintiles of deprivation
  tab2$edi2 = tab2$edi3 = tab2$edi4 = tab2$edi5 = 0   # Create dummy variable
  tab2$edi2[tab2$EDI == 2] <- 1                       # Assign patient with EDI of 2 to be dummy variable 2
  tab2$edi3[tab2$EDI == 3] <- 1                       # As above
  tab2$edi4[tab2$EDI == 4] <- 1                       # As above
  tab2$edi5[tab2$EDI == 5] <- 1                       # As above
  
  # Interaction between stage and comorbidity
  tab2$IntStageCmb  <- tab2$stage*tab2$cmb
  
  # Sort by ID of the patients
  tab2 <- tab2[order(tab2$Ident),]
  
  
  ########
  # Step 3
  ########
  
  print("STEP 3 ")
  print(date())
  
  
  # Return the data frame you have created
  return(list(tab2 = tab2, n=n, NCluster=NCluster, cens.admin = cens.admin, ydiagmin=ydiagmin, ydiagmax=ydiagmax, maxptime=maxptime))
  
}

# Function to simulate the survival time and vital status
cdatasimulationT1WeibOptim <- function (mydat, dist.T1="Weibull", lambda.weib=NULL, rho.weib=NULL, k1=NULL, r1=NULL, a1=NULL,
                                        mu=c(0,0), Sigma=matrix(c(0,0,0,0),2,2),
                                        NPH=FALSE, lambda.weib2=NULL, rho.weib2=NULL, myexpmort,
                                        BetaAge=0, Betasex=0, BetaTTT=0, Betadepriv=0,
                                        BetaEdi2=0, BetaEdi3=0, BetaEdi4=0, BetaEdi5=0,
                                        BetaEthnic=0, BetaCmb=0, BetaRoute=0, BetaStage=0, BetaL1,
                                        BetaAge2=0, BetaAgeSqPlus=0,
                                        BetaStageCmb=0, seed=1234) {
  
  set.seed(seed)
  
  # **************************************************************************** #
  # mydat : Corresponds to the list returned by cDataDesignOptim to collect
  #          * tab2 = Design matrix (NEEDS TO BE ORDERED BY Ident)
  #          * and general information : n, NCluster, cens.admin, ydiagmin, ydiagmax
  # dist.T1                : Distribution of T1 "Weibull" or "GenWeibull"
  # lambda.weib, rho.weib  : scale & shape parameters for the Weibull baseline hazard
  # k1, r1, a1             : parameters for the Generalised Weibull Distribution
  # mu                     : mean of the multinormal random effect
  # Sigma                  : Variance-covariance matrix of the random Effect
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
  # BetaTTT     : Treatment effect (Ref=placebo)
  #                If BetaTTT == 0    => No treatment effect on the cancer mortality hazard (i.e. T1)
  # Betadepriv  : Deprivation effect (Built at the cluster level): Higher level of depriv means higher deprivation
  #                If Betadepriv == 0 => No deprivation effect on the cancer mortality hazard (i.e. T1)
  # BetaEthnic  : Ethnicity effect (ref White)
  #                If BetaEthnic == 0 => No ethnicity effect on the cancer mortality hazard (i.e. T1)
  # BetaCmb     : Comorbidity effect (ref 'no comorbidity')
  #                If BetaCmb    == 0 => No comorbidity effect on the cancer mortality hazard (i.e. T1)
  # BetaRoute   : Route effect (ref 'non-emergency route')
  #                If BetaRoute  == 0 => No route effect on the cancer mortality hazard (i.e. T1)
  # BetaStage   : Stage effect (ref 'early stage')
  #                If BetaStage  == 0 => No stage effect on the cancer mortality hazard (i.e. T1)
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
  if (dist.T1=="GenWeibull" & (is.null(k1) | is.null(r1) | is.null(a1))) {stop("Arguments k1, r1 and a1 should be given for Genearlised Weibull simulation!")}
  if (NPH==TRUE & (is.null(lambda.weib2) | is.null(rho.weib2))) {stop("Arguments lambda.weib2 and rho.weib2 should be given for NPH simulation!")}
  
  ########
  # Step 1
  ########
  
  print("STEP 1 ")
  print(date())
  
  # The list you created in cDataDesignOptim forms the cancer dataset you will use
  tab2 <- mydat$tab2
  
  
  # Collect Information useful from the Design
  # Number of obs
  n <- mydat$n
  # Number of clusters
  NCluster <- mydat$NCluster
  # Year of administrative censoring 
  cens.admin <- mydat$cens.admin
  # Earliest entry year in study
  ydiagmin <- mydat$ydiagmin
  # Latest entry year in study
  ydiagmax <- mydat$ydiagmax
  # Maximum potential time of follow up
  maxptime <- mydat$maxptime
  
  
  # The linear model 
  tab2$exp.xbeta <- exp(BetaAge * tab2$agecr + Betasex * tab2$IsexH +       # Age and sex
                          BetaTTT * tab2$TTT + Betadepriv * tab2$depriv +   # Treatment
                          BetaEdi2 * tab2$edi2 + BetaEdi3 * tab2$edi3 +     # Deprivation
                          BetaEdi4 * tab2$edi4 + BetaEdi5 * tab2$edi5 +     # Deprivation (continued)
                          BetaEthnic * tab2$ethnic + BetaCmb * tab2$cmb +   # Ethnicity and Comorbidity
                          BetaRoute * tab2$route + BetaStage * tab2$stage + # Route and stage
                          BetaL1 * tab2$L1 +                                # Unmeasured variable
                          BetaAge2 * tab2$agecr2 +                          # Non-linear age
                          BetaAgeSqPlus * tab2$agecr2plus +
                          BetaStageCmb * tab2$IntStageCmb                   # Interaction of Stage and comorbidity
  ) 
  
  tab2$exp.Mxbeta   <- exp(-BetaAge * tab2$agecr - Betasex * tab2$IsexH -     
                             BetaTTT * tab2$TTT - Betadepriv * tab2$depriv -
                             BetaEdi2 * tab2$edi2 - BetaEdi3 * tab2$edi3 -     # Deprivation
                             BetaEdi4 * tab2$edi4 - BetaEdi5 * tab2$edi5 -    # Deprivation (continued)
                             BetaEthnic * tab2$ethnic - BetaCmb * tab2$cmb -
                             BetaRoute * tab2$route - BetaStage * tab2$stage -
                             BetaL1 * tab2$L1 -                                # Unmeasured variable
                             BetaAge2 * tab2$agecr2 -
                             BetaAgeSqPlus * tab2$agecr2plus - 
                             BetaStageCmb * tab2$IntStageCmb)
  
  # Time from diagnosis to censoring
  tab2$potentime   <- tab2$cens.admin - tab2$year
  
  # Contribution to the random effect for each cluster
  wi.unique <- mvrnorm(NCluster, mu=mu, Sigma=Sigma)
  
  # Variance and SE?
  tab2$wi <- wi.unique[tab2$cluster,1]
  tab2$vi <- wi.unique[tab2$cluster,2]
  
  ########
  # Step 2
  ########
  
  print("STEP 2 ")
  print(date())
  
  # Survival Time from cancer T1  (Weibull)
  # Random number of n values
  temp.ut1  <- runif(n)
  # If NPH is included in the model
  if (NPH==FALSE){
    # If Wiebull or GenWeibull is specified
    if (dist.T1=="Weibull"){ tab2$T1  <- exp((1/rho.weib) * (log(-1/lambda.weib*log(temp.ut1)) - log(tab2$exp.xbeta) - tab2$wi - tab2$vi*tab2$TTT )) }
    else if (dist.T1=="GenWeibull"){tab2$T1  <-  1/r1*(a1*((1/(1-temp.ut1))^(1/(a1*(tab2$exp.xbeta*exp(tab2$wi + tab2$vi*tab2$TTT))))-1))^(1/k1)}
  }
  else if (NPH==TRUE)
  { tab2$T1  <- (tab2$IsexH)*exp((1/rho.weib) * (log(-1/lambda.weib*log(temp.ut1)) - log(tab2$exp.xbeta) - tab2$wi - tab2$vi*tab2$TTT ))+
    (1-tab2$IsexH)*exp((1/rho.weib2) * (log(-1/lambda.weib2*log(temp.ut1)) - log(tab2$exp.xbeta) - tab2$wi - tab2$vi*tab2$TTT ))         }
  
  # ???
  tab2$T1 <- ifelse(tab2$T1<=.Machine$double.eps^0.5, tab2$T1+10*.Machine$double.eps^0.5, tab2$T1)
  tab2$T1 <- ifelse(is.na(tab2$T1), tab2$maxptime, tab2$T1)
  
  # Survival Time from expected mortality T2
  # Minimum of time from cancer or time to censoring
  tab2$T1cens <- pmin(tab2$T1, tab2$potentime)
  # Flag those who died due to cancer (i.e., cancer time less than censoring time)
  tab2$cause1 <- ifelse(tab2$T1<tab2$potentime,1,0)
  # Age at censoring
  tab2$ageout <- tab2$age+tab2$T1cens
  # Year at censoring
  tab2$yearout <- tab2$year+tab2$T1cens
  
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
  trylxsplit$CMBRT <- trylxsplit$cmb
  
  
  # trylxsplit <- rename(trylxsplit, c("age.entier"="AGEX" ,"year.entier"="year","sex"="sex"))
  
  # Join the data sets together
  tab5 <- join(data.frame(trylxsplit), myexpmort, by = c("AGERT","YEARRT","SEXRT","CMBRT"), match="first")
  
  # Calculation of T2bis (Used to calculate the "final" T2)
  # Calculate the n values for the expected mortality rate with rate equal to 'exprate' in the life table data
  tab5$T2bis <- rexp(nrow(tab5), tab5$exprate)
  
  
  
  
  ########
  # Step 3
  ########
  
  print("STEP 3 ")
  print(date())
  
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
  finalsurvtime  <- pmin(tabtemp$T1, tabtemp$T2)
  
  # Join the data sets of: cancer data, time to other cause death, flag of other cause death, and actual final survival time
  datasimulation <- data.frame(tab2, T2 = tabtemp$T2, cause2=tabtemp$cause2, finalsurvtime)
  
  
  ########
  # Step 4
  ########
  
  
  print("STEP 4 ")
  print(date())
  
  # Vital Status ("vstatus")  : 1 ==> death, 0 ==> alive
  
  # Create a flag variable that holds the vital status if final survival time was the actual survival time
  datasimulation[, "vstatus"] <- ifelse(datasimulation[, "finalsurvtime"] < datasimulation[, "potentime"], 1 , 0)
  
  # Cause of death : 1 ==> Cancer, 2 ==> Other causes
  
  # Create a flag variable for the cause (cancer or otherwise)
  datasimulation$cause <-( 0
                           + 1 * (datasimulation$finalsurvtime == datasimulation$T1 & datasimulation$vstat == 1)
                           + 2 * (datasimulation$finalsurvtime == datasimulation$T2 & datasimulation$vstat == 1) )
  
  # Update of the yearout and ageout to apply for finalsurvtime	and not T1 only
  datasimulation$yearout <- datasimulation$year + datasimulation$finalsurvtime
  datasimulation$ageout <- datasimulation$age + datasimulation$finalsurvtime
  
  # ----- datasimulation : Dataframe finale ----- #
  datasimulation$truncyearout <- trunc(datasimulation$yearout)
  datasimulation$truncageout <- trunc(datasimulation$ageout)
  datasimulation <- merge(datasimulation, myexpmort,
                          by.x=c("truncyearout","truncageout","sex","cmb"),
                          by.y=c("YEARRT","AGERT","SEXRT","CMBRT"))
  
  
  ########
  # Step 5
  ########
  
  print("STEP 5 ")
  print(date())
  
  # Return the data set that holds the variables
  return(datasimulation)
  
}

# Utility function for parametric g-formula estimators 
  # Jessica G. Young, Mats J. Stensrud, Eric J. Tchetgen Tchetgen, Miguel A. Hernán  (2020):
  # "A causal framework for classical statistical estimands in failure-time settings with competing events"
  # https://doi.org/10.1002/sim.8471
calculateCumInc <- function(inputData, timepts=cutTimes, competing=FALSE){
  
  # Create a matrix to hold the data set
  cumulativeIncidence <- matrix(NA, ncol = length(unique(inputData$patno)), nrow = length(timepts))
  
  # Insert the event probabilities at the first time interval
  # Needs to account for "temporal order" Dk+1 before Yk+1
  if(!competing) cumulativeIncidence[1,] <- inputData[inputData$dtime==0,]$hazardP*(1-inputData[inputData$dtime==0,]$hazardO)
  else cumulativeIncidence[1,] <- inputData[inputData$dtime==0,]$hazardO
  
  # Create a matrix compatible with 'cumulativeIncidence' with survival probabilities at each time 
  survivalProb <- t(aggregate(s~patno, data = inputData, FUN = cumprod)$s) #split the long data into subsets per subject
  for(i in 2:length(timepts)){
    subInputDataP <- inputData[inputData$dtime==(i-1),]$hazardP #OBS: dtime starts at 0 - hazard of cancer
    subInputDataO <- inputData[inputData$dtime==(i-1),]$hazardO #OBS: dtime starts at 0 - hazard of other cause
    if(!competing) cumulativeIncidence[i,] <- subInputDataP * (1-subInputDataO) * survivalProb[(i-1),] # OBS: survivalProb entry i is time point i-1
    else cumulativeIncidence[i,] <- subInputDataO * survivalProb[(i-1),]
  }
  meanCumulativeIncidence <- rowMeans(apply(cumulativeIncidence, MARGIN = 2,cumsum)) 
  return(meanCumulativeIncidence)
}