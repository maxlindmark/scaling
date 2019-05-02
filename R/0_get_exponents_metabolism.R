#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - This code estimates size-scaling exponents from data on the rate if it was not provided in the original paper
# 
# A. Load libraries
#
# B. 
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

min(subset(dat, Species == "Rainbow trout")$M)
max(subset(dat, Species == "Rainbow trout")$M)

# 5C
d5 <- dat %>% filter(Species == "Rainbow trout" & Temp == 5)
head(d5)
plot(log(d5$Consumption) ~ log(d5$M))
summary(lm(log(d5$Consumption) ~ log(d5$M)))
summary(lm(log(d5$Consumption) ~ log(d5$M)))$coefficients[2] 
summary(lm(log(d5$Consumption) ~ log(d5$M)))$coefficients[1]

# 10C
d10 <- dat %>% filter(Species == "Rainbow trout" & Temp > 9 & Temp < 10.3)
head(d10)
plot(log(d10$Consumption) ~ log(d10$M))
summary(lm(log(d10$Consumption) ~ log(d10$M)))
summary(lm(log(d10$Consumption) ~ log(d10$M)))$coefficients[2] 
summary(lm(log(d10$Consumption) ~ log(d10$M)))$coefficients[1]

# 15C
d15 <- dat %>% filter(Species == "Rainbow trout" & Temp == 15)
head(d15)
plot(log(d15$Consumption) ~ log(d15$M))
summary(lm(log(d15$Consumption) ~ log(d15$M)))
summary(lm(log(d15$Consumption) ~ log(d15$M)))$coefficients[2] 
summary(lm(log(d15$Consumption) ~ log(d15$M)))$coefficients[1]

# 21C
d21 <- dat %>% filter(Species == "Rainbow trout" & Temp > 19 & Temp < 23)
head(d21, 20)
plot(log(d21$Consumption) ~ log(d21$M))
summary(lm(log(d21$Consumption) ~ log(d21$M)))
summary(lm(log(d21$Consumption) ~ log(d21$M)))$coefficients[2] 
summary(lm(log(d21$Consumption) ~ log(d21$M)))$coefficients[1]

# 24C
d24 <- dat %>% filter(Species == "Rainbow trout" & Temp == 24.3)
head(d24, 20)
plot(log(d24$Consumption) ~ log(d24$M))
summary(lm(log(d24$Consumption) ~ log(d24$M)))
summary(lm(log(d24$Consumption) ~ log(d24$M)))$coefficients[2] 
summary(lm(log(d24$Consumption) ~ log(d24$M)))$coefficients[1]

d2 <- data.frame(b = c(0.6524, 0.5772, 0.74361, 0.74619, 0.70172), 
                 temp = c(24, 21, 15, 10, 5))

ggplot(d2, aes(temp, b)) + 
  geom_point(size = 5) +
  theme_bw(base_size = 18) + 
  stat_smooth(method = "lm") +
  NULL







