#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.03: Max Lindmark
#
# - Explore body growth rate data
# 
# A. Load libraries & read data
#
# B. Fit LMM
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES & READ DATA ========
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
          "RColorBrewer",
          "lmerTest",
          "magrittr")

# Install packages
#install.packages(pkgs)

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# x <- devtools::session_info(pkgs = pkgs)
# x <- as.data.frame(x$packages)
# x <- dplyr::filter(x, package %in% pkgs) %>% 
# dplyr::select(-`*`, -date, -source) %>% 
# dplyr::arrange(package)
# x

# package   version
# 1         dplyr   0.8.0.1
# 2       ggplot2     3.0.0
# 3      lmerTest     3.0-1
# 4      magrittr       1.5
# 5          plyr     1.8.4
# 6  RColorBrewer     1.1-2
# 7         RCurl 1.95-4.10
# 8        readxl     1.3.1
# 9       tidylog     0.1.0
# 10        tidyr     0.8.3
# 11      viridis     0.5.1

#==** Read data ====
# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

dat <- read_excel("data/growth_data.xlsx")

glimpse(dat)

cols = c(1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)


#======== B. FILTER DATA AND FIT LINEAR MIXED MODELS ========
#==** Filter & reshape data ====
# Mass for analysis
dat$mass <- dat$geom_mean_mass_g

# Replace NA with -9...
dat$mass[is.na(dat$mass)] <- -9

dat$mass <- ifelse(dat$mass == -9,
                   dat$mid_mass_sizecl_g,
                   dat$mass)

dat$mass[is.na(dat$mass)] <- -9

dat <- dat %>% filter(mass > 0)

# Normalize mass with respect to max mass
dat$mass_norm <- dat$mass / dat$w_max_published_g

# Filter n>1 species
s_dat <- dat %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Calculate mean optimum temperature within species
s_dat$mean_opt_temp_c <- ave(s_dat$opt_temp_c, s_dat$common_name)

# Center each size class' optimum temperature by mean optimum temperature for that species
s_dat$opt_temp_c_ct <- s_dat$opt_temp_c - s_dat$mean_opt_temp_c

# Define more readable palette
nb.cols <- length(unique(s_dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat, aes(log10(mass_norm), opt_temp_c_ct, 
                  color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = mycolors) +  
  NULL

#==** Fit models ====
# Starting already with the mixed effects model, on conceptual grounds (pseudoreplication)
s_dat$log10_mass_n <- log10(s_dat$mass_norm)

me1 <- lmer(opt_temp_c_ct ~ log10_mass_n + (log10_mass_n | species), s_dat)
summary(me1)

# Even the simplest random effects model results in singularity. Not much I can do about that then!
test <- lmer(opt_temp_c_ct ~ 1 + (1|species), s_dat)

# https://stats.stackexchange.com/questions/393901/singular-fit-with-simplest-random-structure-in-lmer-lme4-is-a-bayesian-approa

# Fitting fixed (intecepts) effects only
m1 <- lm(opt_temp_c_ct ~ log10_mass_n + species, s_dat)
summary(m1)

s_dat$pred <- predict(m1)

ggplot(s_dat, aes(log10(mass_norm), opt_temp_c_ct, 
                  color = common_name)) + 
  geom_point(size = 5) +
  geom_line(data = s_dat, aes(log10(mass_norm), pred), 
            size =3, alpha = 0.5) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = mycolors) +  
  NULL

# Will not bother with interaction here because I hardly have the amount of data needed for that to be meaningful...
