#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - This code estimates size-scaling exponents from data on the rate if it was not 
#   provided in the original paper
# 
# A. Load libraries
#
# B. Fit log-log model by species
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#== A. LOAD LIBRARIES & READ DATA ==================================================
rm(list = ls())

# *** Not finished!

# Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "readxl",
          "RCurl",
          "magrittr")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
pkg_info(pkgs)

# package   version
# 1   dplyr   0.8.0.1
# 2   RCurl 1.95-4.10
# 3  readxl     1.3.1
# 4 tidylog     0.1.0
# 5   tidyr     0.8.3

# Will crate a csv that one can read once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/metabolism_data.xlsx"))

dat <- read_excel("data/consumption_data.xlsx")

glimpse(dat)

cols = c(1, 2, 3, 14, 15, 16, 17)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)


#== B. FIT LOG-LOG MODEL BY SPECIES ================================================

# Check which species have "data" as source in exponent data set and redo them here

#==** Coregonus albula =============================================================
c_albula <- dat %>% 
  filter(species == "Coregonus albula")

sort(unique(c_albula$temp_c))

# 4C
d4 <- c_albula %>% filter(temp_c == 4)
d4
plot(log(d4$metabolic_rate) ~ log(d4$mass_g))
summary(lm(log(d4$metabolic_rate) ~ log(d4$mass_g)))
