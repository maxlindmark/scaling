#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.06: Max Lindmark
#
# - Code to fit non-linear (polynomial) models to maximum consumption rate, 
#   and to do WAIC model selection, using only species with optimum data
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

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       
# readxl_1.3.1      [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1
# tidyr_1.0.0        dplyr_0.8.3       [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file(s)
con <- read.csv("data/con_analysis.csv")

# Find species with optimum data
opt <- con %>% filter(above_optimum == "Y") %>% droplevels()

opt_spec <- unique(opt$species)

# There is a lot of variation in rates between species (see exploratory script. 
# Here we will fit models where y is relative to max for that species
con <- con %>% 
  dplyr::group_by(species) %>% 
  #dplyr::mutate(y_norm = y/max(y)) %>% 
  dplyr::mutate(y_norm = y/mean(y)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(temp_norm_ct = temp_norm - mean(temp_norm)) %>% 
  dplyr::filter(species %in% opt_spec) %>% 
  dplyr::mutate(species_n = as.numeric(as.factor(species))) %>% 
  droplevels()

unique(con$species)

# Temp-range used for prediction
temp_pred = seq(from = min(con$temp_norm_ct), 
                to = max(con$temp_norm_ct),
                length.out = 50)

# Data list for consumption model
data = list(
  # y = con$y_norm,
  # y = con$y,
  y = log(con$y_norm), 
  #y = log(con$y), 
  n_obs = length(con$y_norm), 
  mass = con$log_mass_norm_ct,
  temp = con$temp_norm_ct,
  temp_pred = temp_pred,
  species_n = con$species_n,
  small = 2.5,
  medium = 0,
  large = 2.5
)


# C. MODEL SELECTION ===============================================================
# Here we fit models with different hierarcial structures
# Specifically, we consider:

# -- Without interaction
# M1 - b0 
# M2 - b1
# M3 - b2 + b3 # If b3 random then b2 - both temperature coefficients 

# -- With interaction
# M4 - b0 
# M5 - b1
# M6 - b2 + b3 # If b3 random then b2 - both temperature coefficients 
# M7 - b4


#**** M1 ===========================================================================
model = "R/analysis/non-linear model/models/m1.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m1 <- WAIC


#**** M2 ===========================================================================
model = "R/analysis/non-linear model/models/m2.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m2 <- WAIC


#**** M3 ===========================================================================
model = "R/analysis/non-linear model/models/m3.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m3 <- WAIC


#**** M4 ===========================================================================
model = "R/analysis/non-linear model/models/m4.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m4 <- WAIC


#**** M5 ===========================================================================
model = "R/analysis/non-linear model/models/m5.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m5 <- WAIC


#**** M6 ===========================================================================
model = "R/analysis/non-linear model/models/m6.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m6 <- WAIC


#**** M7 ===========================================================================
model = "R/analysis/non-linear model/models/m7.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj = jags.samples(jm, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd <- -2*sum(log(summary(zj$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m7 <- WAIC


#** COMPARE WAIC ===================================================================
# WAIC
waic_m1
waic_m2
waic_m3
waic_m4
waic_m5
waic_m6
waic_m7

# > waic_m1     # Selected
# [1] 633.0699
# > waic_m2
# [1] 841.8353
# > waic_m3
# [1] 759.1478
# > waic_m4
# [1] 635.7621
# > waic_m5
# [1] 843.4833
# > waic_m6
# [1] 759.5018
# > waic_m7
# [1] 838.4903



