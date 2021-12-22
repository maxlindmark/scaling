#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model to metabolism using a log-linear model 
# inspired by the MTE, and perform model selection using WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Model selection
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(rjags)
library(RColorBrewer)
library(ggmcmc)
library(RCurl)
library(readxl)
library(magrittr)
library(viridis)
library(patchwork)
library(bayesplot)

# sessionInfo()
# other attached packages:
# [1] bayesplot_1.7.2    patchwork_1.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_2.0.1     readxl_1.3.1       RCurl_1.98-1.2    
# [8] ggmcmc_1.4.1       ggplot2_3.3.2      tidyr_1.1.2        dplyr_1.0.2        RColorBrewer_1.1-2 rjags_4-10         coda_0.19-4   


# B. READ IN DATA ==================================================================
# Read in your data file
dat <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

str(dat)

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_peak_temp == "N")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Mean-center predictor variables
dat$log_mass_ct <- dat$log_mass - mean(dat$log_mass)
dat$temp_arr_ct <- dat$temp_arr - mean(dat$temp_arr)

# Use mass-specific values
#dat$y_spec <- dat$y / dat$mass_g

# Add dummy variable for metabolic type (0 for routine/resting, 1 for standard)
dat$met_type <- ifelse(dat$type == "Standard", 1, 0)
dat$met_s <- ifelse(dat$type == "Standard", 1, 0)
dat$met_r <- ifelse(!dat$type == "Standard", 1, 0)

plot(dat$met_s ~ dat$met_r)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Arrange by species name so there's no confusion between the index and the name
dat <- dat %>% arrange(species_n)

# Data in list-format for JAGS
data = list(
  y = log(dat$y), 
  n_obs = length(dat$y), 
  species_n = dat$species_n,
  mass = dat$log_mass_ct,
  temp = dat$temp_arr_ct,
  met_r = dat$met_r,
  met_s = dat$met_s,
  spec_s = distinct(dat, species, .keep_all = TRUE)$met_s, # New variable to go in the 1:34 for loop 
  spec_r = distinct(dat, species, .keep_all = TRUE)$met_r # New variable to go in the 1:34 for loop 
)


# C. MODEL SELECTION ===============================================================
# Here we fit models with different hierarchical structures
# Specifically, we consider:

# M1  - all coefficients vary by species
# M2  - intercept, mass, temperature vary by species
# M3a - intercept and mass vary by species
# M3b - intercept and temperature vary by species
# M4  - intercept varies by species
# M5  - no interaction, all coefficients vary by species
# M6a - no interaction, intercept and mass vary by species
# M6b - no interaction, intercept and temperature vary by species
# M7  - no interaction, intercept varies by species

# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

#**** M1 ===========================================================================
# M1  - all coefficients vary by species

model1 = "JAGS_models/log_linear/metabolism/m1.txt"

# Manually set initial values, because otherwise all the chains get the same
# NOTE I don't do it for all parameters...
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    mu_b3 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    sigma_b3 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples, see JAGS 4.3 user manual
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    mu_b3 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    sigma_b3 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    mu_b3 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    sigma_b3 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm1 = jags.model(model1,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

# Just to check the specified initial values are used!
# https://stats.stackexchange.com/questions/231787/rjags-does-not-seem-to-use-initial-values-specified
# jm1 = jags.model(model, data = data, n.adapt = 0, n.chains = 3, inits = inits)
# jm1$state()
# coda.samples(jm1, c('mu_b0'), n.iter = 1)
# jm1$state()[[1]]$mu_b0

update(jm1, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj1 = jags.samples(jm1, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd1 <- -2*sum(log(summary(zj1$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC1 <- sum((summary(zj1$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m1 <- lppd1 + 2*pd.WAIC1


#**** M2 ===========================================================================
# M2  - intercept, mass, temperature vary by species

model2 = "JAGS_models/log_linear/metabolism/m2.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    b3 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm2 = jags.model(model2,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm2, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj2 = jags.samples(jm2, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd2 <- -2*sum(log(summary(zj2$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC2 <- sum((summary(zj2$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m2 <- lppd2 + 2*pd.WAIC2


#**** M3a ==========================================================================
# M3a - intercept and mass vary by species

model3a = "JAGS_models/log_linear/metabolism/m3a.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    mu_b1 = 0.1,
    b2 = 0.1,
    b3 = 0.11,
    sigma = 0.1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    mu_b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    mu_b1 = 2,
    b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b1 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm3a = jags.model(model3a,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm3a, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj3a = jags.samples(jm3a, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd3a <- -2*sum(log(summary(zj3a$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC3a <- sum((summary(zj3a$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m3a <- lppd3a + 2*pd.WAIC3a


#**** M3b ==========================================================================
# M3b - intercept and temperature vary by species

model3b = "JAGS_models/log_linear/metabolism/m3b.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    b1 = 0.1,
    mu_b2 = 0.1,
    b3 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    b1 = 1,
    mu_b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    b1 = 2,
    mu_b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm3b = jags.model(model3b,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm3b, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj3b = jags.samples(jm3b, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd3b <- -2*sum(log(summary(zj3b$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC3b <- sum((summary(zj3b$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m3b <- lppd3b + 2*pd.WAIC3b


#**** M4 ===========================================================================
# M4  - intercept varies by species

model4 = "JAGS_models/log_linear/metabolism/m4.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    b1 = 2,
    b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm4 = jags.model(model4,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm4, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj4 = jags.samples(jm4, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd4 <- -2*sum(log(summary(zj4$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC4 <- sum((summary(zj4$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m4 <- lppd4 + 2*pd.WAIC4


#**** M5 ===========================================================================
# M5  - no interaction, all coefficients vary by species

model5 = "JAGS_models/log_linear/metabolism/m5.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm5 = jags.model(model5,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm5, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj5 = jags.samples(jm5, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd5 <- -2*sum(log(summary(zj5$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC5 <- sum((summary(zj5$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m5 <- lppd5 + 2*pd.WAIC5


#**** M6a ===========================================================================
# M6a - no interaction, intercept and mass vary by species

model6a = "JAGS_models/log_linear/metabolism/m6a.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    mu_b1 = 0.1,
    b2 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b1 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    mu_b1 = 1,
    b2 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    mu_b1 = 2,
    b2 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b1 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm6a = jags.model(model6a,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm6a, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj6a = jags.samples(jm6a, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd6a <- -2*sum(log(summary(zj6a$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC6a <- sum((summary(zj6a$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m6a <- lppd6a + 2*pd.WAIC6a


#**** M6b ===========================================================================
# M6b - no interaction, intercept and temperature vary by species
model6b = "JAGS_models/log_linear/metabolism/m6b.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    b1 = 0.1,
    mu_b2 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    b1 = 1,
    mu_b2 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    b1 = 2,
    mu_b2 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm6b = jags.model(model6b,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm6b, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj6b = jags.samples(jm6b, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd6b <- -2*sum(log(summary(zj6b$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC6b <- sum((summary(zj6b$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m6b <- lppd6b + 2*pd.WAIC6b


#**** M7 ===========================================================================
# M7  - no interaction, intercept varies by species
model7 = "JAGS_models/log_linear/metabolism/m7.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    b1 = 1,
    b2 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    b1 = 2,
    b2 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm7 = jags.model(model7,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm7, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj7 = jags.samples(jm7, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd7 <- -2*sum(log(summary(zj7$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC7 <- sum((summary(zj7$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m7 <- lppd7 + 2*pd.WAIC7


#** COMPARE WAIC ===================================================================
# WAIC
waic_m1
waic_m2
waic_m3a
waic_m3b
waic_m4
waic_m5
waic_m6a
waic_m6b
waic_m7

# WAIC suggests model 1, closely followed by model 2
# > waic_m1
# [1] 274.5855
# > waic_m2
# [1] 275.7643
# > waic_m3a
# [1] 580.3847
# > waic_m3b
# [1] 659.8612
# > waic_m4
# [1] 923.1542
# > waic_m5
# [1] 280.6449
# > waic_m6a
# [1] 622.775
# > waic_m6b
# [1] 661.1844
# > waic_m7
# [1] 956.0734

# Calculate delta WAIC
waic_m1 - waic_m1
waic_m2 - waic_m1
waic_m3a - waic_m1
waic_m3b - waic_m1
waic_m4 - waic_m1
waic_m5 - waic_m1
waic_m6a - waic_m1
waic_m6b - waic_m1
waic_m7 - waic_m1

# > waic_m1 - waic_m1
# [1] 0
# > waic_m2 - waic_m1
# [1] 1.178716
# > waic_m3a - waic_m1
# [1] 305.7992
# > waic_m3b - waic_m1
# [1] 385.2756
# > waic_m4 - waic_m1
# [1] 648.5686
# > waic_m5 - waic_m1
# [1] 6.059374
# > waic_m6a - waic_m1
# [1] 348.1895
# > waic_m6b - waic_m1
# [1] 386.5988
# > waic_m7 - waic_m1
# [1] 681.4878
