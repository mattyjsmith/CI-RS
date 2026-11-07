# 01_setup_packages.R
# Load and install required packages

required_packages <- c(
  "progress", "dplyr", "MASS", "Epi", "plyr", "readxl", "mexhaz",
  "utils", "boot", "gtools", "cobalt", "twang", "ggplot2"
)

# install.packages(c("progress", "dplyr", "MASS", "Epi", "plyr", "readxl", "mexhaz", "utils", "boot", "gtools", "cobalt", "twang", "ggplot2"))


# Function to install missing packages
load_or_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Apply to all required packages
invisible(lapply(required_packages, load_or_install))