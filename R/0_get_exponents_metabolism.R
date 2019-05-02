#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - This code estimates size-scaling exponents from data on the rate if it was not provided in the original paper
# 
# A. Load libraries
#
# B. Fit log-log model by species
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES & READ DATA ========
rm(list = ls())

#====**** Load packages ========
## Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "readxl",
          "RCurl")

## Install packages
#install.packages(pkgs)

## Load all packages
lapply(pkgs, library, character.only = TRUE)

## Print package version
# x <- devtools::session_info(pkgs = pkgs)
# x <- as.data.frame(x$packages)
# x <- dplyr::filter(x, package %in% pkgs) %>% 
# dplyr::select(-`*`, -date, -source) %>% 
# dplyr::arrange(package)
# x

# package   version
# 1   dplyr   0.8.0.1
# 2   RCurl 1.95-4.10
# 3  readxl     1.3.1
# 4 tidylog     0.1.0
# 5   tidyr     0.8.3

#====**** Read data ========
# Will crate a csv that one can read once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/metabolism_data.xlsx"))

dat <- read_excel("data/metabolism_data.xlsx", col_types = "guess")

cols = c(1, 2, 3, 14, 15, 16, 17)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)


#======== B. FIT LOG-LOG MODEL BY SPECIES ========
#====**** Coregonus albula ========
c_albula <- dat %>% 
  filter(species == "Coregonus albula")

sort(unique(c_albula$temp_c))

# 4C
d4 <- c_albula %>% filter(temp_c == 4)
d4
plot(log(d4$metabolic_rate) ~ log(d4$mass_g))
summary(lm(log(d4$metabolic_rate) ~ log(d4$mass_g)))

# 8C
d8 <- c_albula %>% filter(temp_c == 8)
d8
plot(log(d8$metabolic_rate) ~ log(d8$mass_g))
summary(lm(log(d8$metabolic_rate) ~ log(d8$mass_g)))

# 15C
d15 <- c_albula %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))

#====**** Coregonus fontanae ========
c_fontanae <- dat %>% 
  filter(species == "Coregonus fontanae")

sort(unique(c_fontanae$temp_c))

# 4C
d4 <- c_fontanae %>% filter(temp_c == 4)
d4
plot(log(d4$metabolic_rate) ~ log(d4$mass_g))
summary(lm(log(d4$metabolic_rate) ~ log(d4$mass_g)))

# 8C
d8 <- c_fontanae %>% filter(temp_c == 8)
d8
plot(log(d8$metabolic_rate) ~ log(d8$mass_g))
summary(lm(log(d8$metabolic_rate) ~ log(d8$mass_g)))

# 15C
d15 <- c_fontanae %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))

#====**** Abramis brama ========
a_brama <- dat %>% 
  filter(species == "Abramis brama")

sort(unique(a_brama$temp_c))

# 5C
d5 <- a_brama %>% filter(temp_c == 5)
d5
plot(log(d5$metabolic_rate) ~ log(d5$mass_g))
summary(lm(log(d5$metabolic_rate) ~ log(d5$mass_g)))

# 10C
d10 <- a_brama %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g)))

# 15C
d15 <- a_brama %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))

# 20C
d20 <- a_brama %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g)))

# 23C
d23 <- a_brama %>% filter(temp_c == 23)
d23
plot(log(d23$metabolic_rate) ~ log(d23$mass_g))
summary(lm(log(d23$metabolic_rate) ~ log(d23$mass_g)))

#====**** Rutilus rutilus ========
# Not sure yet if I should use this paper or stick to the same as for Cmax...

#====**** Salvelinus confluentus ========



