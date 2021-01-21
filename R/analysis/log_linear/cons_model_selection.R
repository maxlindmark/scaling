#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model to maximum consumption using a log-linear model 
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
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

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
dat$y_spec <- dat$y / dat$mass_g

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
data = list(
  y = log(dat$y_spec), 
  n_obs = length(dat$y), 
  species_n = dat$species_n,
  mass = dat$log_mass_ct,
  temp = dat$temp_arr_ct
)


# C. MODEL SELECTION ===============================================================
# Here we fit models with different hierarcial structures
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

model1 = "JAGS_models/log_linear/growth_consumption/m1.txt"

# Manually set initial values, because otherwise all the chains get the same
# NOTE I don't do it for all parameters...
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    mu_b3 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    sigma_b3 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples, see JAGS 4.3 user manual
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    mu_b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    sigma_b3 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    mu_b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
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

update(jm1, n.iter = n.iter) 

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

model2 = "JAGS_models/log_linear/growth_consumption/m2.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    b3 = 1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm2 = jags.model(model2,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm2, n.iter = n.iter) 

# Monitor the likelihood to calculate WAIC
zj2 = jags.samples(jm2, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = n.iter, 
                  thin = thin,
                  inits = inits)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd2 <- -2*sum(log(summary(zj2$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC2 <- sum((summary(zj2$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m2 <- lppd2 + 2*pd.WAIC2


# Since this model is similar in WAIC to the model without interaction, I want to check the estimate of it
test <- coda.samples(jm2,
                     variable.names = c("b3"),
                     n.iter = n.iter,
                     thin = thin)

summary(test)


#**** M3a ==========================================================================
# M3a - intercept and mass vary by species

model3a = "JAGS_models/log_linear/growth_consumption/m3a.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    b2 = 0.1,
    b3 = 1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm3a = jags.model(model3a,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

burn.in = 10000 # Length of burn-in

update(jm3a, n.iter = n.iter) 

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

model3b = "JAGS_models/log_linear/growth_consumption/m3b.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    mu_b2 = 0.1,
    b3 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    mu_b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    b1 = 2,
    mu_b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm3b = jags.model(model3b,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

burn.in = 10000 # Length of burn-in

update(jm3b, n.iter = n.iter) 

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

model4 = "JAGS_models/log_linear/growth_consumption/m4.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    b1 = 2,
    b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm4 = jags.model(model4,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

burn.in = 10000 # Length of burn-in

update(jm4, n.iter = n.iter) 

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

model5 = "JAGS_models/log_linear/growth_consumption/m5.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm5 = jags.model(model5,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

burn.in = 10000 # Length of burn-in

update(jm5, n.iter = n.iter) 

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

model6a = "JAGS_models/log_linear/growth_consumption/m6a.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    b2 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    b2 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    b2 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm6a = jags.model(model6a,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

burn.in = 10000 # Length of burn-in

update(jm6a, n.iter = n.iter) 

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
model6b = "JAGS_models/log_linear/growth_consumption/m6b.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    mu_b2 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    mu_b2 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    b1 = 2,
    mu_b2 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm6b = jags.model(model6b,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

burn.in = 10000 # Length of burn-in

update(jm6b, n.iter = n.iter) 

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
model7 = "JAGS_models/log_linear/growth_consumption/m7.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    sigma = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    b1 = 2,
    b2 = 2,
    sigma = 2,
    sigma_b0 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm7 = jags.model(model7,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

burn.in = 10000 # Length of burn-in

update(jm7, n.iter = n.iter) 

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

# WAIC suggests model 5 is best fitting with model 2 as a close runner up
# > # WAIC
# > waic_m1
# [1] 564.5484
# > waic_m2
# [1] 563.3664
# > waic_m3a
# [1] 708.3846
# > waic_m3b
# [1] 630.3754
# > waic_m4
# [1] 750.1536
# > waic_m5
# [1] 560.2499
# > waic_m6a
# [1] 726.3189
# > waic_m6b
# [1] 634.3972
# > waic_m7
# [1] 774.189

# Calculate delta WAIC
waic_m1 - waic_m5
waic_m2 - waic_m5
waic_m3a - waic_m5
waic_m3b - waic_m5
waic_m4 - waic_m5
waic_m5 - waic_m5
waic_m6a - waic_m5
waic_m6b - waic_m5
waic_m7 - waic_m5

# > waic_m1 - waic_m5
# [1] 4.298438
# > waic_m2 - waic_m5
# [1] 3.116463
# > waic_m3a - waic_m5
# [1] 148.1346
# > waic_m3b - waic_m5
# [1] 70.12543
# > waic_m4 - waic_m5
# [1] 189.9037
# > waic_m5 - waic_m5
# [1] 0
# > waic_m6a - waic_m5
# [1] 166.069
# > waic_m6b - waic_m5
# [1] 74.14731
# > waic_m7 - waic_m5
# [1] 213.9391
