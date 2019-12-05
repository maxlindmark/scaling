#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of metabolic rate as a function of 
# temperature with different group-effects and compare DIC and WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Model selection
#
# D. Model validation
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
dat <- read.csv("data/met_analysis.csv")

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_optimum == "N")

str(dat)
summary(dat)
glimpse(dat)

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Plot response vs explanatory - color temp
ggplot(dat, aes(log_mass_norm_ct, log(y), color = temp_norm_arr_ct)) + 
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 11) +
  scale_color_viridis(discrete = FALSE, option = "plasma") +
  theme(aspect.ratio = 1) +
  labs(x = "ln(standardized mass)",
       y = "ln(metabolic rate)",
       color = "Standardized\nArrhenius\ntemperature") +
  NULL
#ggsave("figs/supplement/metabolism_mass_scatter.pdf", plot = last_plot(), scale = 1, width = 14, height = 14, units = "cm", dpi = 300)

# Plot response vs explanatory - color species
ggplot(dat, aes(log_mass_norm_ct, log(y), color = species)) + 
  geom_point(size = 2, alpha = 0.7) +
  #  stat_smooth(method = "lm", size = 2, alpha = 0.4, color = "black") +
  theme_classic(base_size = 11) +
  scale_color_viridis(discrete = TRUE, option = "plasma") +
  guides(color = FALSE) +
  theme(aspect.ratio = 1) +
  labs(x = "ln(standardized mass)",
       y = "ln(metabolic rate)",
       color = "Standardized\nArrhenius\ntemperature") +
  NULL

# Plot response vs explanatory - color mass
ggplot(dat, aes(temp_norm_arr_ct, log(y), color = log_mass_norm_ct)) + 
  geom_point(size = 2, alpha = 0.7) +
  #  stat_smooth(method = "lm", size = 2, alpha = 0.4, color = "black") +
  theme_classic(base_size = 11) +
  scale_color_viridis(discrete = FALSE, option = "plasma") +
  theme(aspect.ratio = 1) +
  labs(x = "Standardized Arrhenius temperature",
       y = "ln(metabolic rate)",
       color = "ln(standardized mass)") +
  NULL
#ggsave("figs/supplement/metabolism_temp_scatter.pdf", plot = last_plot(), scale = 1, width = 14, height = 14, units = "cm", dpi = 300)

# Plot no. data points per species
ggplot(dat, aes(species_ab)) + 
  geom_bar(size = 5) +
  theme_classic(base_size = 15) +
  scale_fill_viridis(discrete = TRUE, option = "plasma") +  
  theme(axis.text.x = element_text(angle = 90)) +
  NULL

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Mass-range used for prediction
mass_pred = seq(from = min(dat$log_mass_norm_ct), 
                to = max(dat$log_mass_norm_ct),
                length.out = 100)

# Data in list-format for JAGS
data = list(
  y = log(dat$y), 
  n_obs = length(dat$y), 
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

#**** M1 ===========================================================================
cat(
  "model{
  
  for(i in 1:n_obs){
     y[i] ~ dnorm(mu[i], tau)
     mu[i] <- 
       b0[species_n[i]] +                # varying intercept
       b1[species_n[i]]*mass[i] +        # varying mass-exponent 
       b2[species_n[i]]*temp[i] +        # varying activation energy
       b3[species_n[i]]*mass[i]*temp[i]  # varying M*T interaction
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # global mean
  mu_b1 ~ dnorm(-0.25, 0.5)  # global mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # global activation energy
  mu_b3 ~ dnorm(0, 0.5)      # global interaction
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  sigma_b3 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  tau_b3 <- 1/sigma_b3^2
  
}", fill = TRUE, file = "R/analysis/model_selection/m1_metabolism.txt")

model = "R/analysis/model_selection/m1_metabolism.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Extract DIC
dic_m1 <- dic.samples(jm, n.iter = 2500, type = "pD")
dic_m1

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


#**** M2 =============================================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1[species_n[i]]*mass[i] +   # varying mass-exponent
      b2[species_n[i]]*temp[i] +   # varying activation energy
      b3*mass[i]*temp[i]           # non-varying M*T interaction
  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  #-- Priors	
  b3 ~ dnorm(0, 0.5)         # global interaction
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # varying mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying activation energy
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/model_selection/m2_metabolism.txt")

model = "R/analysis/model_selection/m2_metabolism.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

# Extract DIC
dic_m2 = dic.samples(jm, n.iter = 2500, type = "pD")
dic_m2

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


#**** M3a =============================================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +          # varying intercept 
      b1[species_n[i]]*mass[i] +  # varying mass-exponent
      b2*temp[i] +                # non-varying activation energy
      b3*mass[i]*temp[i]          # non-varying M*T interaction
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  
  #-- Priors	
  b2 ~ dnorm(-0.6, 0.5)
  b3 ~ dnorm(0, 0.5)
  mu_b0 ~ dnorm(0, 0.5)
  mu_b1 ~ dnorm(-0.25, 0.5)
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  
  }", fill = TRUE, file = "R/analysis/model_selection/m3a_metabolism.txt")

model = "R/analysis/model_selection/m3a_metabolism.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

dic_m3a = dic.samples(jm, n.iter = 2500, type = "pD")
dic_m3a

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


#**** M3b =============================================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +          # varying intercept 
      b1*mass[i] +                # non-varying mass-exponent
      b2[species_n[i]]*temp[i] +  # varying activation energy
      b3*mass[i]*temp[i]          # non-varying M*T interaction
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  #-- Priors	
  b1 ~ dnorm(-0.25, 0.5)
  b3 ~ dnorm(0, 0.5)
  mu_b0 ~ dnorm(0, 0.5)
  mu_b2 ~ dnorm(-0.6, 0.5)
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/model_selection/m3b_metabolism.txt")

model = "R/analysis/model_selection/m3b_metabolism.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

dic_m3b = dic.samples(jm, n.iter = 2500, type = "pD")
dic_m3b

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


#**** M4 =============================================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
      mu[i] <- 
      b0[species_n[i]] + # varying intercept
      b1*mass[i] +       # non-varying mass exponent
      b2*temp[i] +       # non-varying activation energy
      b3*mass[i]*temp[i] # non-varying M*T interaction
  # Add log likelihood computation for each observation!
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
  b0[j] ~ dnorm(mu_b0, tau_b0)
  }
  
  #-- Priors	
  b1 ~ dnorm(-0.25, 0.5)
  b2 ~ dnorm(-0.6, 0.5)
  b3 ~ dnorm(0, 0.5)
  mu_b0 ~ dnorm(0, 0.5)
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  
  }", fill = TRUE, file = "R/analysis/model_selection/m4_metabolism.txt")

model = "R/analysis/model_selection/m4_metabolism.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

dic_m4 = dic.samples(jm, n.iter = 2500, type = "pD")
dic_m4

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


#**** M5 =============================================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1[species_n[i]]*mass[i] +   # varying mass-exponent
      b2[species_n[i]]*temp[i]     # varying activation energy
  
  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # varying mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying activation energy
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/model_selection/m5_consumption.txt")

model = "R/analysis/model_selection/m5_consumption.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

dic_m5 = dic.samples(jm, n.iter = 2500, type = "pD")
dic_m5

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


#** COMPARE DIC & WAIC =============================================================
# WAIC
waic_m1
waic_m2
waic_m3a
waic_m3b
waic_m4
waic_m5

# > waic_m1
# [1] 531.4119
# > waic_m2
# [1] 529.0335
# > waic_m3a
# [1] 674.0834
# > waic_m3b
# [1] 587.8616
# > waic_m4
# [1] 710.8867

# DIC
dic_m1
dic_m2
dic_m3a
dic_m3b
dic_m4

# Both waic and DIC suggest model 2 is most parsimonious


# D. MODEL VALIDATION ==============================================================
model = "R/analysis/model_selection/m2_metabolism.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 

#** Sample from the posterior ======================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm,
                   variable.names = c("b0", "b1", "b2", 
                                      "mu_b0", "mu_b1", "mu_b2", 
                                      "sigma_b0", "sigma_b1", "sigma_b2",
                                      "sigma", "b3"), # Population level
                   n.iter = samples, 
                   thin = n.thin)

summary(cs) # Get the mean estimate and SE and 95% CIs

cs_df <- data.frame(summary(cs)[1])
cs_df$Parameter <- row.names(cs_df)

# Get global temperature-coefficient
cs_df %>% filter(Parameter == "mu_b2") %>% select(statistics.Mean)

# Get global mass-coefficient
cs_df %>% filter(Parameter == "mu_b1") %>% select(statistics.Mean)


#** Evaluate convergence ===========================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

# Plot posterior densities of species intercepts
cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",  "b0[14]",
                          "b0[15]", "b0[16]", "b0[17]", "b0[18]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free_y") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",  "b0[14]",
                          "b0[15]", "b0[16]", "b0[17]", "b0[18]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 4, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL

# Plot posterior densities of species mass-effects
cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",  "b1[14]",
                          "b1[15]", "b1[16]", "b1[17]", "b1[18]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2) +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",  "b1[14]",
                          "b1[15]", "b1[16]", "b1[17]", "b1[18]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 4, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL

# Plot posterior densities of temperature-effects
cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]", "b2[14]",
                          "b2[15]", "b2[16]", "b2[17]", "b2[18]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2) +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]", "b2[14]",
                          "b2[15]", "b2[16]", "b2[17]", "b2[18]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 4, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL

# Plot posterior densities of random means, interaction and standard deviations
cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "b3",
                          "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 3, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "b3",
                          "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 4, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL

# Time series of running means
# cs_df %>% 
#   ggs_running(.) + 
#   facet_wrap(~ Parameter, ncol = 4, scales = "free") +
#   geom_line(size = 1.1, alpha = 0.8) +
#   theme_classic(base_size = 11) + 
#   scale_color_brewer(palette = rev("Dark2")) + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   theme(axis.text.x = element_text(size = 6)) +
#   NULL
#ggsave("figs/supplement/growth_running_mean.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)

# Rhat
cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.002) +
  theme_classic(base_size = 11) +
  geom_point(size = 2) +
  theme(aspect.ratio = 1)+
  NULL

#ggsave("figs/supplement/growth_rhat.pdf", plot = last_plot(), scale = 1, width = 14, height = 14, units = "cm", dpi = 300)


