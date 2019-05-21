#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore consumption and metabolic scaling data
# 
# A. Load libraries & read data
#
# B. Explore data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES & READ DATA =============================================
rm(list = ls())

#==** Load packages ====
# Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "readxl",
          "RCurl",
          "ggplot2",
          "viridis",
          "plyr",
          "RColorBrewer")

# Install packages
#install.packages(pkgs)

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
pkg_info(pkgs)

# package   version
# 1        dplyr   0.8.0.1
# 2      ggplot2     3.0.0
# 3         plyr     1.8.4
# 4 RColorBrewer     1.1-2
# 5        RCurl 1.95-4.10
# 6       readxl     1.3.1
# 7      tidylog     0.1.0
# 8        tidyr     0.8.3
# 9      viridis     0.5.1

#==** Read data ====
# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

con <- read_excel("data/consumption_data.xlsx")
met <- read_excel("data/metabolism_data.xlsx")

glimpse(dat)

cols = c(1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)


#======== B. EXPLORE DATA ===========================================================