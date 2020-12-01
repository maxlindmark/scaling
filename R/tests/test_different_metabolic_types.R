#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of metabolic and maximum consumption
# rate as a function of temperature and mass
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit models 
#    - After this section you can run the others separately
#
# D. Model validation (convergence, fit, residuals)
#    - This part is LONG and contains lots of big plots. Skip if you don't want that
# 
# E. Plot predictions
# 
# F. Additional calculations on the posterior
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(RColorBrewer)
library(ggmcmc)
library(RCurl)
library(readxl)
library(magrittr)
library(viridis)
library(patchwork)
library(bayesplot)
library(scales)
library(bayesplot)
library(tidylog)
library(brms)

# B. READ IN DATA ==================================================================
# Read in your data file
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

# Filter data points at below optimum temperatures
met <- met %>% filter(above_peak_temp == "N")

# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species_ab))

# Mean-center predictor variables
met$log_mass_ct <- met$log_mass - mean(met$log_mass)
met$temp_arr_ct <- met$temp_arr - mean(met$temp_arr)

# Use mass-specific values
met$y_spec <- met$y / met$mass_g

# Masses for prediction
mass_pred_met <-  seq(from = min(met$log_mass_ct), 
                      to = max(met$log_mass_ct),
                      length.out = 100)

# Prepare data for JAGS
met_data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
met_data = list(
  y = log(met$y_spec), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_ct,
  temp = met$temp_arr_ct,
  mass_pred = mass_pred_met)


# Fit model to see if there is a difference in slope or intercept for different metabolic types brms

met$met_type <- "standard"
met$met_type <- ifelse(met$type == "Standard", "standard", "rout_rest")

unique(met$met_type)

m <- brm(log(y_spec) ~ log_mass_ct*met_type + temp_arr_ct -1 + (1 + met_type | species_n), data = met,
         inits = "random", chains = 4, iter = 2000, family = gaussian(),
         control = list(max_treedepth = 15, adapt_delta = 0.99))

summary(m)

posterior <- as.array(m)

p1 <- mcmc_areas(
  posterior,
  pars = c("b_met_typerout_rest", "b_met_typestandard"),
  prob = 0.95,
  point_est = "mean"
) + ggtitle("Difference in intercept")

p2 <- mcmc_areas(
  posterior,
  pars = c("b_log_mass_ct:met_typestandard"),
  prob = 0.95,
  point_est = "mean"
) + ggtitle("Difference in slope")

p1/p2
