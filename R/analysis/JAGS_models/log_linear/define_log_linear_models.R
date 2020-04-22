#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define and save models for rjags for model selection using WAIC
# 
# A. Specify models
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. SPECIFY MODELS ================================================================
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

#**** M1 ===========================================================================
# All coefficients vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
     y[i] ~ dnorm(mu[i], tau)
     mu[i] <- 
       b0[species_n[i]] +                # varying intercept
       b1[species_n[i]]*mass[i] +        # varying mass-exponent 
       b2[species_n[i]]*temp[i] +        # varying temperature coefficient
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
  mu_b2 ~ dnorm(-0.6, 0.5)   # global temperature coefficient
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
  
}", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m1.txt")


#**** M2 ===========================================================================
# Intercept, mass, temperature vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1[species_n[i]]*mass[i] +   # varying mass-exponent
      b2[species_n[i]]*temp[i] +   # varying temperature coefficient
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
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying temperature coefficient
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m2.txt")


#**** M3a =========================================================================
# Intercept and mass vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +          # varying intercept 
      b1[species_n[i]]*mass[i] +  # varying mass-exponent
      b2*temp[i] +                # non-varying temperature coefficient
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
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m3a.txt")


#**** M3b ==========================================================================
# Intercept and temperature vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +          # varying intercept 
      b1*mass[i] +                # non-varying mass-exponent
      b2[species_n[i]]*temp[i] +  # varying temperature coefficient
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
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m3b.txt")


#**** M4 ===========================================================================
# Intercept varies by species
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
      mu[i] <- 
      b0[species_n[i]] + # varying intercept
      b1*mass[i] +       # non-varying mass exponent
      b2*temp[i] +       # non-varying temperature coefficient
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
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m4.txt")


#**** M5 ===========================================================================
# No interaction, full random
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1[species_n[i]]*mass[i] +   # varying mass-exponent
      b2[species_n[i]]*temp[i]     # varying temperature coefficient
  
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
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying temperature coefficient
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m5.txt")


#**** M6a ==========================================================================
# No interaction, intercept and mass random
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1[species_n[i]]*mass[i] +   # varying mass-exponent
      b2*temp[i]                   # non-varying temperature coefficient
  
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
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # varying mass-exponent
  b2 ~ dnorm(-0.6, 0.5)      # non-varying temperature coefficient
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m6a.txt")


#**** M6b ==========================================================================
# No interaction, intercept and temperature random
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1*mass[i] +                 # non-varying mass-exponent
      b2[species_n[i]]*temp[i]     # varying temperature coefficient
  
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
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  b1 ~ dnorm(-0.25, 0.5)     # non-varying mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying temperature coefficient
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m6b.txt")


#**** M7 ===========================================================================
# No interaction, intercept random
cat(
  "model{
  
  for(i in 1:n_obs){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1*mass[i] +                 # non-varying mass effect
      b2*temp[i]                   # non-varying temperature coefficient
  
  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    }
  
  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  b1 ~ dnorm(-0.25, 0.5)
  b2 ~ dnorm(-0.6, 0.5)
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/log_linear/m7.txt")
