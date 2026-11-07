# 02_load_lifetable.R
# Load and adjust life table data for deprivation

# Load life table
expratesASD <- as.data.frame(readxl::read_excel(here::here("data", "LT.xlsx")))

# Keep only least (1) and most (5) deprived groups
expratesASD <- expratesASD[expratesASD$DEPRT %in% c(1, 5), ]

# Recode deprivation: 0 = least deprived, 1 = most deprived
expratesASD$DEPRT <- ifelse(expratesASD$DEPRT == 1, 0, 1)

# Define factor for modifying deprivation effect
factor <- 1  # You can increase this >1 to amplify deprivation differences

# Compute average mortality across deprivation groups
expratesASD$expratesMean <- NA
for (i in unique(expratesASD$YEARRT)) {
  for (j in unique(expratesASD$SEXRT)) {
    for (k in unique(expratesASD$AGERT)) {
      idx <- expratesASD$YEARRT == i & expratesASD$SEXRT == j & expratesASD$AGERT == k
      expratesASD$expratesMean[idx] <- mean(expratesASD$exprate[idx], na.rm = TRUE)
    }
  }
}

# Adjust exprate for increased deprivation effect
expratesASD$expratesOther <- abs(expratesASD$exprate - expratesASD$expratesMean) + expratesASD$expratesMean
expratesASD$exprateAdj <- ifelse(
  expratesASD$DEPRT == 1,
  expratesASD$exprate * factor,
  expratesASD$exprate + ((1 - factor) * expratesASD$expratesOther)
)

# Overwrite exprate with adjusted version
expratesASD$exprate <- expratesASD$exprateAdj