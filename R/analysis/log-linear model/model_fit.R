#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of maximum consumption and metabolic rate as a function of 
# temperature with group-effects and assess model fit
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit models
# 
# D. Evaluate model fit
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
library(tidybayes)

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file(s)
met <- read.csv("data/met_analysis.csv")
con <- read.csv("data/con_analysis.csv")

# Filter data points at below optimum temperatures
met <- met %>% filter(above_optimum == "N")
con <- con %>% filter(above_optimum == "N")

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(met$species, 1, 1)
sp2 <- gsub( ".*\\s", "", met$species )
met$species_ab <- paste(sp1, sp2, sep = ".")

sp1 <- substring(con$species, 1, 1)
sp2 <- gsub( ".*\\s", "", con$species )
con$species_ab <- paste(sp1, sp2, sep = ".")

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species))
con$species_n <- as.numeric(as.factor(con$species))

# Mass-range used for prediction
mass_pred_met = seq(from = min(met$log_mass_norm_ct), 
                    to = max(met$log_mass_norm_ct),
                    length.out = 100)

mass_pred_con = seq(from = min(con$log_mass_norm_ct), 
                    to = max(con$log_mass_norm_ct),
                    length.out = 100)


# Data in list-format for JAGS
met_data = list(
  y = log(met$y), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_norm_ct,
  temp = met$temp_norm_arr_ct,
  mass_pred = mass_pred_met
)

con_data = list(
  y = log(con$y), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_norm_ct,
  temp = con$temp_norm_arr_ct,
  mass_pred = mass_pred_con
)


# C. FIT MODELS ====================================================================
# Refit chosen models from the model selection part
# Need to modify them to track new variables for assessing model fit
#**** Metabolism (M2) ==============================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
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
  
  # Predictions
  for(k in 1:length(mass_pred)){
      
    pred_warm[k] <- mu_b0 + mu_b1*mass_pred[k] + mu_b2*-1.5 + b3*mass_pred[k]*-1.5
    pred_cold[k] <- mu_b0 + mu_b1*mass_pred[k] + mu_b2*1.5 + b3*mass_pred[k]*1.5
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b3 ~ dnorm(0, 0.5)         # non-varying interaction
  mu_b0 ~ dnorm(0, 0.5)      # global intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # global mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # global activation energy
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/m2_metabolism_pred.txt")

met_model = "R/analysis/m2_metabolism_pred.txt"

jm_met = jags.model(met_model,
                    data = met_data, 
                    n.adapt = 5000, 
                    n.chains = 3)


#**** Consumption (M5) =============================================================
# Predictions
cat(
  "model{
  
    for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
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

    for(k in 1:length(mass_pred)){
  
        pred_warm[k] <- mu_b0 + mu_b1*mass_pred[k] + mu_b2*-1.5
        pred_cold[k] <- mu_b0 + mu_b1*mass_pred[k] + mu_b2*1.5
  
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # global intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # global mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # global activation energy
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/m5_consumption_pred.txt")

model_con = "R/analysis/m5_consumption_pred.txt"

jm_con = jags.model(model_con,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3)


# D. EVALUATE MODEL FIT ============================================================
# Plot mean of simulated data vs mean of observed data
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# First convert your matrix 
cs_fit_met = coda.samples(jm_met,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

cs_fit_con = coda.samples(jm_con,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

# Convert to data frames
cs_fit_df_met <- data.frame(as.matrix(cs_fit_met))
cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df_met)))

# Metabolism
p3 <- ggplot(cs_fit_df_met, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_met$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df_met$p_mean), digits = 3))) +
  labs(x = "Mean simulated metabolism", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  ggtitle("Metabolism") +
  NULL

p4 <- ggplot(cs_fit_df_met, aes(cv_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_met$cv_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3, 
           label = paste("P =", round(mean(cs_fit_df_met$p_cv), digits = 3))) +
  labs(x = "cv simulated metabolism", y = "") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) +
  NULL

# Consumption
p5 <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_con$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "C", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df_con$p_mean), digits = 3))) +
  labs(x = "Mean simulated consumption", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  ggtitle("Consumption") +
  NULL

p6 <- ggplot(cs_fit_df_con, aes(cv_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_con$cv_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "D", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -2, vjust = 3.3, 
           label = paste("P =", round(mean(cs_fit_df_con$p_cv), digits = 3))) +
  labs(x = "cv simulated consumption", y = "") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) +
  NULL

(p3 + p4) / (p5 + p6)
#ggsave("figures/supp/cv_mean_fit.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)
