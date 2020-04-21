#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchicals model to maximum consumption using a log-linear model 
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
dat <- read.csv("data/con_analysis.csv")
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

# Data in list-format for JAGS
data = list(
  y = log(dat$y), 
  n_obs = length(dat$y), 
  species_n = dat$species_n,
  mass = dat$log_mass_ct,
  #mass = dat$log_mass_norm_ct,
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
model = "R/analysis/log-linear model/models/m1.txt"

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
# [1] 103.3267


#**** M2 ===========================================================================
model = "R/analysis/log-linear model/models/m2.txt"

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
#[1] 102.3396

# TEST: Checking estimates:
cs <- coda.samples(jm,
                   variable.names = c("b3"),
                   n.iter = 10000, 
                   thin = 5)

summary(cs)

js = jags.samples(jm, 
                  variable.names = c("b3"), 
                  n.iter = 10000, 
                  thin = 5)

ecdf(js$b3)(0) 

cs %>% mcmc_dens() 


#**** M3a ==========================================================================
model = "R/analysis/log-linear model/models/m3a.txt"

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
#[1] 101.667


#**** M3b ==========================================================================
model = "R/analysis/log-linear model/models/m3b.txt"

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
# [1] 88.79128

#**** M4 ===========================================================================
model = "R/analysis/log-linear model/models/m4.txt"

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
#[1] 91.02917


#**** M5 ===========================================================================
model = "R/analysis/log-linear model/models/m5.txt"

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
# [1] 102.6251

#**** M6 ===========================================================================
model = "R/analysis/log-linear model/models/m6.txt"

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
#[1] 93.48178


#**** M7 ===========================================================================
model = "R/analysis/log-linear model/models/m7.txt"

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
#[1] 59.10719


#**** M8 ===========================================================================
model = "R/analysis/log-linear model/models/m8.txt"

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

# Standard error of WAIC
n_cases <- nrow(data.frame(data))
lppd_ind <- log(summary(zj$pd, mean)$stat)
pd.WAIC_ind <- (summary(zj$log_pd, sd)$stat)^2
waic_vec <- -2*(lppd_ind - pd.WAIC_ind)
sqrt(n_cases*var(waic_vec))
#[1] 48.22182


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

# > waic_m1
# [1] 529.8837
# > waic_m2
# [1] 528.5802
# > waic_m3a
# [1] 659.2791
# > waic_m3b
# [1] 587.4765
# > waic_m4
# [1] 695.2167
# > waic_m5
# [1] 526.6564
# > waic_m6
# [1] 711.8353
# > waic_m7
# [1] 1125.48
# > waic_m8
# [1] 1349.298

# Calculate delta WAIC
waic_m1 - waic_m5
waic_m2 - waic_m5
waic_m3a - waic_m5
waic_m3b - waic_m5
waic_m4 - waic_m5
waic_m5 - waic_m5
waic_m6 - waic_m5
waic_m7 - waic_m5
waic_m8 - waic_m5

# > waic_m1 - waic_m5
# [1] 3.520652
# > waic_m2 - waic_m5
# [1] 2.038109
# > waic_m3a - waic_m5
# [1] 132.4823
# > waic_m3b - waic_m5
# [1] 47.96525
# > waic_m4 - waic_m5
# [1] 157.4132
# > waic_m5 - waic_m5
# [1] 0
# > waic_m6 - waic_m5
# [1] 185.4572
# > waic_m7 - waic_m5
# [1] 681.9205
# > waic_m8 - waic_m5
# [1] 672.0204

# WAIC suggests model 5 is best fitting
