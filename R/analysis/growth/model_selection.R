#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical models to growth rate rate using a log-linear model 
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

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file(s)
dat <- read.csv("data/growth_analysis.csv")

str(dat)

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_optimum == "N")

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Mass-range used for prediction
mass_pred = seq(from = min(dat$log_mass_norm_ct), 
                to = max(dat$log_mass_norm_ct),
                length.out = 100)

# Filter only positive growth rates
dat <- dat %>% filter(., G > 0)

# Data in list-format for JAGS
data = list(
  y = log(dat$G), 
  n_obs = length(dat$G), 
  species_n = dat$species_n,
  mass = dat$log_mass_norm_ct,
  temp = dat$temp_norm_arr_ct
)


# C. MODEL SELECTION ===============================================================
# Here we fit models with different hierarcial structures
# Specifically, we consider:

# M1 - all coefficients vary by species
# M2 - intercept, mass, temperature vary by species
# M3a - intercept and mass vary by species
# M3b - intercept and temperature vary by species
# M4 - intercept vary by species
# M5 no interaction, full random
# M6 no interaction, intercept random
# M7 no interaction, mass random
# M8 no interaction, temperature random

#**** M1 ===========================================================================
model = "R/analysis/growth/models/m1.txt"

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
model = "R/analysis/growth/models/m2.txt"

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


#**** M3a ===========================================================================
model = "R/analysis/growth/models/m3a.txt"

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
waic_m3a <- WAIC


#**** M3b ===========================================================================
model = "R/analysis/growth/models/m3b.txt"

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
waic_m3b <- WAIC


#**** M4 ===========================================================================
model = "R/analysis/growth/models/m4.txt"

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
model = "R/analysis/growth/models/m5.txt"

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
model = "R/analysis/growth/models/m6.txt"

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
model = "R/analysis/growth/models/m7.txt"

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


#**** M8 ===========================================================================
model = "R/analysis/growth/models/m8.txt"

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
waic_m8 <- WAIC

#** COMPARE WAIC ===================================================================
# WAIC
waic_m1
waic_m2
waic_m3a
waic_m3b
waic_m4
waic_m5
waic_m6
waic_m7
waic_m8

# WAIC suggests model 1 is most parsimonious

#-- With interaction
# > waic_m1
# [1] 32.51117

# > waic_m2
# [1] 37.45695

# > waic_m3a
# [1] 72.0526

# > waic_m3b
# [1] 78.72478

# > waic_m4
# [1] 100.5671

#-- No interaction
# > waic_m5
# [1] 34.71809
# 
# > waic_m6
# [1] 98.97725
# 
# > waic_m7
# [1] 313.4243
# 
# > waic_m8
# [1] 337.2472

# M1 - all coefficients vary by species
# M2 - intercept, mass, temperature vary by species
# M3a - intercept and mass vary by species
# M3b - intercept and temperature vary by species
# M4 - intercept vary by species

# M5 no interaction, full random
# M6 no interaction, intercept random
# M7 no interaction, mass random
# M8 no interaction, temperature random

