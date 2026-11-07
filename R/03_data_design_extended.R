# 03_data_design.R
# Function to simulate baseline data and treatment assignment

# Function to simulate the study design
cDataDesignOptim <- function (n = 500, cens.admin = 2015, ydiagmin=2005, ydiagmax=2010,
                              betaagecr = 0.5, betasex = 0.8, betadep = 0.5, betastage = 0.0, age_shift = 0)
{
  
  #set.seed(seed)
  
  # **************************************************************************** #
  # n					   : Sample size
  # cens.admin : Year of administrative censoring (if == 1995 : Maximum potential time of F-Up maxptime = 15. If 1991, maxptime = 11)
  #              Needs to be greater than 1990
  # ydiagmin, ydiagmax : Min and max of the period of diagnosis (parameter linked to cens.admin)
  
  # EXAMPLE
  # toto <- cDataDesignOptim(n = 10000, cens.admin = 1995, ydiagmin=1980, ydiagmax=1981)
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
  IsexH <- ifelse(sex==1,0,1)  # Males 0, Females 1
  
  # Year of diagnosis
  year <- runif(n, min = ydiagmin, max = ydiagmax)
  
  # Age at diagnosis (continuous)
  age   <- c(runif(n*0.25, min = 30, max = 65), runif(n*0.35, min = 65, max = 75), runif(n*0.40, min = 75, max = 85)) + age_shift   # age_shift > 0: an older population
  
  # Create dep
  dep <- rbinom(n,size=1,prob=0.5) # 1 = deprived
  
  # Create stage
  stage <- rbinom(n,size=1,prob=0.40) # 1 = deprived
  
  
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
  #        dep        : Deprivation level (higher=more deprived)
  
  # Create a data frame to hold the variables you have created
  tab2              <- data.frame(Ident, sex, IsexH, year, age, dep, stage)
  
  # Create centered and rescaled versions of the age variable
  tab2$agec         <- tab2$age - 70 - age_shift   # centred on the shifted age distribution
  tab2$agecr        <- tab2$agec/10
  tab2$maxptime     <- rep(maxptime, n)
  tab2$cens.admin   <- rep(cens.admin, n)
  
  # Sort by ID of the patients
  tab2 <- tab2[order(tab2$Ident),]
  
  # Treatment
  tab2$trt <- rbinom(n,1, prob=plogis(0 + betaagecr*tab2$agecr + betasex*tab2$IsexH + betadep*tab2$dep + betastage*tab2$stage))
  
  ########
  # Step 3
  ########
  
  print("STEP 3 ")
  print(date())
  
  
  # Return the data frame you have created
  return(list(tab2 = tab2, n=n, cens.admin = cens.admin, ydiagmin=ydiagmin, ydiagmax=ydiagmax, maxptime=maxptime))
  
}
