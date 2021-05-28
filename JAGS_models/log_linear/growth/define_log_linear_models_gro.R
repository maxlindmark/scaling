#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define and save models for rjags for model selection using WAIC
# - These models are used for growth and consumption data. For metabolism
# - data we in addition have a fixed effect of metabolic type (routine or standard),
# - so those models are written in another script.
# 
# A. Specify models for model selection
# 
# B. Add predictions and fit to the selected model code
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. SPECIFY MODELS ================================================================
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

#**** M1 ===========================================================================
# All coefficients vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3[species_n[i]]*mass_intra[i] + b4[species_n[i]]*temp_intra[i] + b5[species_n[i]]*mass_intra[i]*temp_intra[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b3[j] ~ dnorm(mu_b3, tau_b3)
    b4[j] ~ dnorm(mu_b4, tau_b4)
    b5[j] ~ dnorm(mu_b5, tau_b5)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(-0.25, 1)
    mu_b4 ~ dnorm(-0.6, 1)   
    mu_b5 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
    sigma_b4 ~ dunif(0, 10)
    sigma_b5 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
    tau_b4 <- 1/(sigma_b4*sigma_b4)
    tau_b5 <- 1/(sigma_b5*sigma_b5)
  
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m1.txt")


#**** M2 ===========================================================================
# Intercept, mass, temperature vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3[species_n[i]]*mass_intra[i] + b4[species_n[i]]*temp_intra[i] + b5*mass_intra[i]*temp_intra[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b3[j] ~ dnorm(mu_b3, tau_b3)
    b4[j] ~ dnorm(mu_b4, tau_b4)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(-0.25, 1)
    mu_b4 ~ dnorm(-0.6, 1)   
    b5 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
    sigma_b4 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
    tau_b4 <- 1/(sigma_b4*sigma_b4)
    
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m2.txt")


#**** M3a =========================================================================
# Intercept and mass vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3[species_n[i]]*mass_intra[i] + b4*temp_intra[i] + b5*mass_intra[i]*temp_intra[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(-0.25, 1)
    b4 ~ dnorm(-0.6, 1)   
    b5 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
    
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
    
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m3a.txt")


#**** M3b ==========================================================================
# Intercept and temperature vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*mass_intra[i] + b4[species_n[i]]*temp_intra[i] + b5*mass_intra[i]*temp_intra[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b4[j] ~ dnorm(mu_b4, tau_b4)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    b3 ~ dnorm(-0.25, 1)
    mu_b4 ~ dnorm(-0.6, 1)   
    b5 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b4 ~ dunif(0, 10)
    
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b4 <- 1/(sigma_b4*sigma_b4)
    
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m3b.txt")


#**** M4 ===========================================================================
# Intercept varies by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*mass_intra[i] + b4*temp_intra[i] + b5*mass_intra[i]*temp_intra[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    b3 ~ dnorm(-0.25, 1)
    b4 ~ dnorm(-0.6, 1)   
    b5 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m4.txt")


#**** M5 ===========================================================================
# No interaction, full random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3[species_n[i]]*mass_intra[i] + b4[species_n[i]]*temp_intra[i]
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b3[j] ~ dnorm(mu_b3, tau_b3)
    b4[j] ~ dnorm(mu_b4, tau_b4)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(-0.25, 1)
    mu_b4 ~ dnorm(-0.6, 1)   
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
    sigma_b4 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
    tau_b4 <- 1/(sigma_b4*sigma_b4)
  
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m5.txt")


#**** M6a ==========================================================================
# No interaction, intercept and mass random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3[species_n[i]]*mass_intra[i] + b4*temp_intra[i]
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(-0.25, 1)
    b4 ~ dnorm(-0.6, 1)   
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
    
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
    
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m6a.txt")


#**** M6b ==========================================================================
# No interaction, intercept and temperature random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*mass_intra[i] + b4[species_n[i]]*temp_intra[i]
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b4[j] ~ dnorm(mu_b4, tau_b4)
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    b3 ~ dnorm(-0.25, 1)
    mu_b4 ~ dnorm(-0.6, 1)   
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b4 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b4 <- 1/(sigma_b4*sigma_b4)
  
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m6b.txt")


#**** M7 ===========================================================================
# No interaction, intercept random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)            
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*mass_intra[i] + b4*temp_intra[i]
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    
  }
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    b3 ~ dnorm(-0.25, 1)
    b4 ~ dnorm(-0.6, 1)   
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    
}", fill = TRUE, file = "JAGS_models/log_linear/growth/m7.txt")


# B. Add predictions to selected models ============================================
#**** M1 ===========================================================================
# All coefficients vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood    
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3[species_n[i]]*mass_intra[i] + b4[species_n[i]]*temp_intra[i] + b5[species_n[i]]*mass_intra[i]*temp_intra[i]  
  
  # Simulate for comparison with data (evalute fit)
    y_sim[i] ~ dnorm(mu[i], tau)
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b3[j] ~ dnorm(mu_b3, tau_b3)
    b4[j] ~ dnorm(mu_b4, tau_b4)
    b5[j] ~ dnorm(mu_b5, tau_b5)
  }
  
  # Predictions
  for(k in 1:length(mass_pred)){
    pred[k] <- mu_b0 + b1*mass_pred[k] + b2*temp_pred + mu_b3*mass_intra_pred[k] + mu_b4*temp_intra_pred + mu_b5*temp_intra_pred*mass_intra_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
    
    cv_y <- sd(y[])/mean(y[])
    cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
    p_cv <- step(cv_y_sim - cv_y)
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    b1 ~ dnorm(-0.25, 1)  
    b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(-0.25, 1)
    mu_b4 ~ dnorm(-0.6, 1)   
    mu_b5 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
    sigma_b4 ~ dunif(0, 10)
    sigma_b5 ~ dunif(0, 10)

  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
    tau_b4 <- 1/(sigma_b4*sigma_b4)
    tau_b5 <- 1/(sigma_b5*sigma_b5)
  
}", fill = TRUE, file = "JAGS_models/log_linear/selected_models/m1_pred_fit_gro.txt")
