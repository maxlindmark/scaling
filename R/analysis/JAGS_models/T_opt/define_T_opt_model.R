#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define and save models for rjags for model selection using WAIC
# 
# A. Specify models for model selection
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. SPECIFY MODELS ================================================================
# Here we fit models with different hierarcial structures
# Specifically, we consider:

# M1  - intercept and slope vary by species
# M2  - intercept varies by species

#**** M1 Random intercept & slope ==================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[species_n[i]] + b1[species_n[i]]*mass[i] # varying intercept & slope

  # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)

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
  
  # Predictions
  for(k in 1:length(mass_pred)){
    pred[k] <- mu_b0 + mu_b1*mass_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)    
    mu_b1 ~ dnorm(0, 0.04)    
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
    sigma_b1 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_b1 <- 1/(sigma_b1*sigma_b1)
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/T_opt/m1_T_opt_pred_fit.txt")


#**** M2 Random intercept ==========================================================
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood  
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[species_n[i]] + b1*mass[i] # varying intercept fixed slope

  # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)

  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
  }
  
  # Predictions
  for(k in 1:length(mass_pred)){
    pred[k] <- mu_b0 + b1*mass_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)    
    b1 ~ dnorm(0, 0.04)       
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10)
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/T_opt/m2_T_opt_pred_fit.txt")
