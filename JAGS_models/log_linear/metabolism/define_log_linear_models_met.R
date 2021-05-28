#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define and save models for rjags for model selection using WAIC 
# - These models are used for metabolism data
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
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i] + b3[species_n[i]]*mass[i]*temp[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    mu_b1 ~ dnorm(0.75, 1)  
    mu_b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
    sigma_b2 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
    tau_b2 <- 1/(sigma_b2*sigma_b2)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
  
}", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m1.txt")


#**** M2 ===========================================================================
# Intercept, mass, temperature vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i] + b3*mass[i]*temp[i]
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    mu_b1 ~ dnorm(0.75, 1)  
    mu_b2 ~ dnorm(-0.6, 1)   
    b3 ~ dnorm(0, 1)         
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
    sigma_b2 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
    tau_b2 <- 1/(sigma_b2*sigma_b2)
  
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m2.txt")


#**** M3a =========================================================================
# Intercept and mass vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
    
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1[species_n[i]]*mass[i] + b2*temp[i] + b3*mass[i]*temp[i]  
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    mu_b1 ~ dnorm(0.75, 1)
    b2 ~ dnorm(-0.6, 1)
    b3 ~ dnorm(0, 1)
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
  
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m3a.txt")


#**** M3b ==========================================================================
# Intercept and temperature vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1*mass[i] + b2[species_n[i]]*temp[i] + b3*mass[i]*temp[i]
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    b1 ~ dnorm(0.75, 1)
    mu_b2 ~ dnorm(-0.6, 1)
    b3 ~ dnorm(0, 1)
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b2 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b2 <- 1/(sigma_b2*sigma_b2)
  
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m3b.txt")


#**** M4 ===========================================================================
# Intercept varies by species
cat(
  "model{
  
  for(i in 1:n_obs){
    
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1*mass[i] + b2*temp[i] + b3*mass[i]*temp[i] 
      
  # Add log likelihood computation for each observation!
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    b1 ~ dnorm(0.75, 1)
    b2 ~ dnorm(-0.6, 1)
    b3 ~ dnorm(0, 1)
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
  
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m4.txt")


#**** M5 ===========================================================================
# No interaction, full random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood  
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i] 
    
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    mu_b1 ~ dnorm(0.75, 1)          
    mu_b2 ~ dnorm(-0.6, 1)           
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
    sigma_b2 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
    tau_b2 <- 1/(sigma_b2*sigma_b2)
  
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m5.txt")


#**** M6a ==========================================================================
# No interaction, intercept and mass random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1[species_n[i]]*mass[i] + b2*temp[i] 
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    mu_b1 ~ dnorm(0.75, 1)          
    b2 ~ dnorm(-0.6, 1)              
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
    
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m6a.txt")


#**** M6b ==========================================================================
# No interaction, intercept and temperature random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1*mass[i] + b2[species_n[i]]*temp[i] 
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    b0_type ~ dnorm(0, 0.04)
    b1 ~ dnorm(0.75, 1)            
    mu_b2 ~ dnorm(-0.6, 1)          
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b2 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b2 <- 1/(sigma_b2*sigma_b2)
  
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m6b.txt")


#**** M7 ===========================================================================
# No interaction, intercept random
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1*mass[i] + b2*temp[i] 
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
  }
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    b1 ~ dnorm(0.75, 1)
    b2 ~ dnorm(-0.6, 1)
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    
  }", fill = TRUE, file = "JAGS_models/log_linear/metabolism/m7.txt")


# B. Add predictions to selected models ============================================
#**** M1 ===========================================================================
# All coefficients vary by species
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0_r[species_n[i]]*met_r[i] + b0_s[species_n[i]]*met_s[i] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i] + b3[species_n[i]]*mass[i]*temp[i] # this works nicely
  
  # Simulate for comparison with data (evalute fit)
    y_sim[i] ~ dnorm(mu[i], tau)
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
    
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0s[j] ~ dnorm(mu_b0_s, tau_b0_s)
    b0_s[j] <- spec_s[j] * b0s[j]
    b0r[j] ~ dnorm(mu_b0_r, tau_b0_r)
    b0_r[j] <- spec_r[j] * b0r[j]
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Predictions (assuming met_type[i] = 0; i.e. resting or routine metabolic rate)
  for(k in 1:length(mass_pred)){
    pred[k] <- mu_b0_r + mu_b1*mass_pred[k] + mu_b2*temp_pred + mu_b3*temp_pred*mass_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
    
    cv_y <- sd(y[])/mean(y[])
    cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
    p_cv <- step(cv_y_sim - cv_y)
  
  # Priors	
    mu_b0_s ~ dnorm(-2, 0.04)
    mu_b0_r ~ dnorm(-1, 0.04)
    mu_b1 ~ dnorm(0.75, 1)  
    mu_b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(0, 1)      
    sigma ~ dunif(0, 10) 
    sigma_b0_s ~ dunif(0, 10)
    sigma_b0_r ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
    sigma_b2 ~ dunif(0, 10)
    sigma_b3 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0_s <- 1/(sigma_b0_s*sigma_b0_s)
    tau_b0_r <- 1/(sigma_b0_r*sigma_b0_r)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
    tau_b2 <- 1/(sigma_b2*sigma_b2)
    tau_b3 <- 1/(sigma_b3*sigma_b3)
  
}", fill = TRUE, file = "JAGS_models/log_linear/selected_models/m1_pred_fit_met.txt")

