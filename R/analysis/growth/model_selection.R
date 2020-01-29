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
# mass_pred = seq(from = min(dat$log_mass_norm_ct), 
#                 to = max(dat$log_mass_norm_ct),
#                 length.out = 100)

mass_pred = seq(from = min(dat$log_mass_ct), 
                to = max(dat$log_mass_ct),
                length.out = 100)

# Filter only positive growth rates
dat <- dat %>% filter(., G > 0)

# Data in list-format for JAGS
data = list(
  y = log(dat$G), 
  n_obs = length(dat$G), 
  species_n = dat$species_n,
  #mass = dat$log_mass_norm_ct,
  mass = dat$log_mass_ct,
  #temp = dat$temp_norm_arr_ct
  temp = dat$temp_arr - mean(dat$temp_arr)
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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 31.45414

# cs = coda.samples(jm,
#                   variable.names = c(
#                     "mu_b3"), 
#                   n.iter = 10000, 
#                   thin = 5)
# summary(cs)
# js = jags.samples(jm, 
#                   variable.names = c("mu_b3"), 
#                   n.iter = 10000, 
#                   thin = 5)
# 
# ecdf(js$mu_b3)(0) 


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 32.07077


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 31.41895


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 31.28701


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 31.88844


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 32.25407


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 31.70706


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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 25.00297

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
pd.WAIC <- sum((summary(zj$log_pd, sd)$stat)^2) # Penalty (variance across samples)

# WAIC = model fit + 2*penalty
WAIC <- lppd + 2*pd.WAIC

c(pd.WAIC, WAIC)
waic_m8 <- WAIC

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 18.96065


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
# [1] 35.47166
# > waic_m2
# [1] 39.36716
# > waic_m3a
# [1] 70.76162
# > waic_m3b
# [1] 79.2966
# > waic_m4
# [1] 101.2919
# > waic_m5
# [1] 36.33327
# > waic_m6
# [1] 98.58492
# > waic_m7
# [1] 206.925
# > waic_m8
# [1] 223.3752

# delta waic
waic_m1 - waic_m1
waic_m2 - waic_m1
waic_m3a - waic_m1
waic_m3b - waic_m1
waic_m4 - waic_m1
waic_m5 - waic_m1
waic_m6 - waic_m1
waic_m7 - waic_m1
waic_m8 - waic_m1


# M1 - all coefficients vary by species
# M2 - intercept, mass, temperature vary by species
# M3a - intercept and mass vary by species
# M3b - intercept and temperature vary by species
# M4 - intercept vary by species

# M5 no interaction, full random
# M6 no interaction, intercept random
# M7 no interaction, mass random
# M8 no interaction, temperature random

#** COMPARE S.E. OF WAIC ===========================================================
#**** First, example  ==============================================================
# Based on Statistical Rethikning version 1
# Instead of summing the fit and the variance for all observations, as I do with sum
# above to calculate WAIC, I will use the individual observatuins to calculate standard error

# In stat rethinking, waic is calculated with sligthly different code (but same equation)
#lppd <- sapply( 1:n_cases , function(i) log_sum_exp(logprob[i,]) - log(n_samples) )
#pWAIC <- sapply( 1:n_cases , function(i) var(logprob[i,]) )
#-2*( sum(lppd) - sum(pWAIC) )

# I do the summing directly. I can make sure it's the same:
# WAIC for model 8 is:
waic_m8
# [1] 223.3752

# Using Stat-rethinking method, we get (i.e. same variables but not summing:
# (lppd_ind = lppd and pd.WAIC_ind = pWAIC in statistical rethinking)
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
-2*( sum(lppd_ind) - sum(pd.WAIC_ind) )
# [1] 223.3752

# So, now that we know it's the same, we can use his equation to calculate standard
# error of the WAIC
# from stat rethink:
# waic_vec <- -2*( lppd - pWAIC )
# sqrt( n_cases*var(waic_vec) )

# Calculate sample size
n_cases <- nrow(data.frame(data))

waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
# [1] 18.9656

# Stat rethink argues for using the standard error of the WAIC.
# Here I will look at the interval of the waic, not the difference,
# but using the same approach:
# (2.6 corresponding to 99% interval)
waic <- -2*(sum(lppd_ind) - sum(pd.WAIC_ind))
waic_se <- sqrt(n_cases*var(waic_vec))

waic + c(-1, 1) * waic_se * 2.6


#**** Repeat for M1 and M5  ========================================================
# Refit models
# M1
model = "R/analysis/growth/models/m1.txt"

jm1 = jags.model(model,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm1, n.iter = burn.in) 

zj1 = jags.samples(jm1, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = 10000, 
                   thin = 1)

# M5
model = "R/analysis/growth/models/m5.txt"

jm5 = jags.model(model,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm5, n.iter = burn.in) 

zj5 = jags.samples(jm5, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = 10000, 
                   thin = 1)


# Calculate WAIC_se
lppd_ind_1 <- log(summary(zj1$pd, mean)$stat)
pd.WAIC_ind_1 <- (summary(zj1$log_pd, sd)$stat)^2

lppd_ind_5 <- log(summary(zj5$pd, mean)$stat)
pd.WAIC_ind_5 <- (summary(zj5$log_pd, sd)$stat)^2

n_cases <- nrow(data.frame(data))

waic_vec_1 <- -2*(lppd_ind_1 - pd.WAIC_ind_1)
waic_vec_5 <- -2*(lppd_ind_5 - pd.WAIC_ind_5)

sqrt(n_cases*var(waic_vec))
sqrt(n_cases*var(waic_vec))

waic_1 <- -2*(sum(lppd_ind_1) - sum(pd.WAIC_ind_1))
waic_1_se <- sqrt(n_cases*var(waic_vec_1))

waic_5 <- -2*(sum(lppd_ind_5) - sum(pd.WAIC_ind_5))
waic_5_se <- sqrt(n_cases*var(waic_vec_5))

# Compare WAIC and se for two models
# (1.96 corresponding to 95% interval)

waic_1 + (c(-1, 1) * waic_1_se * 1.96)
waic_5 + (c(-1, 1) * waic_5_se * 1.96)




