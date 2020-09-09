#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define quadratic model for consumption rate data that includes peak temperature
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#**** M1 - b0 varies by species, second degree polynomial ==========================
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood  
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
    
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
  for(k in 1:length(temp_pred)){
    pred_large[k] <- mu_b0 + b1*mass_pred_l + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- mu_b0 + b1*mass_pred_s + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)
    b1 ~ dnorm(0, 0.04)
    b2 ~ dnorm(0, 0.04)
    b3 ~ dnorm(0, 0.04)
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10) 
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    
  }", fill = TRUE, file = "R/analysis/JAGS_models/polynomial/polynomial_model1.txt")


#**** M1b - b0 varies by species, second degree polynomial, sigma increase with mass ====
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood  
    y[i] ~ dnorm(mu[i], tau[i])
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
    sigma[i] <- alpha.sigma + b1.sigma*mass[i]
    
    tau[i] <- 1/sigma[i]^2

  # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau[i])
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau[i])
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau[i]))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
  }
    
  # Predictions
  for(k in 1:length(temp_pred)){
    pred_large[k] <- mu_b0 + b1*mass_pred_l + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- mu_b0 + b1*mass_pred_s + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)
    b1 ~ dnorm(0, 0.04)
    b2 ~ dnorm(0, 0.04)
    b3 ~ dnorm(0, 0.04)
    sigma_b0 ~ dunif(0, 10)
    alpha.sigma ~ dgamma(0.001, 0.001) # intercept in model of sigma
    b1.sigma ~ dnorm(0, 0.01)

  # Derived quantiles
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    
  }", fill = TRUE, file = "R/analysis/JAGS_models/polynomial/polynomial_model1b.txt")


#**** M1c - b0 varies by species, second degree polynomial, sigma increase with temp ====
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood  
    y[i] ~ dnorm(mu[i], tau[i])
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
    sigma[i] <- alpha.sigma + b1.sigma*temp[i]
    
    tau[i] <- 1/sigma[i]^2

  # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau[i])
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau[i])
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau[i]))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
  }
    
  # Predictions
  for(k in 1:length(temp_pred)){
    pred_large[k] <- mu_b0 + b1*mass_pred_l + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- mu_b0 + b1*mass_pred_s + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)
    b1 ~ dnorm(0, 0.04)
    b2 ~ dnorm(0, 0.04)
    b3 ~ dnorm(0, 0.04)
    sigma_b0 ~ dunif(0, 10)
    alpha.sigma ~ dgamma(0.001, 0.001) # intercept in model of sigma
    b1.sigma ~ dnorm(0, 0.01)

  # Derived quantiles
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    
  }", fill = TRUE, file = "R/analysis/JAGS_models/polynomial/polynomial_model1c.txt")


#**** M1d - b0 varies by species, second degree polynomial, sigma increase with mass and temp ====
cat(
  "model{
  
  for(i in 1:n_obs){
  
  # Likelihood  
    y[i] ~ dnorm(mu[i], tau[i])
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
    sigma[i] <- alpha.sigma + b1.sigma*temp[i] + b2.sigma*mass[i]
    
    tau[i] <- 1/sigma[i]^2

  # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau[i])
  
  # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau[i])
  
  # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau[i]))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
  }
    
  # Predictions
  for(k in 1:length(temp_pred)){
    pred_large[k] <- mu_b0 + b1*mass_pred_l + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- mu_b0 + b1*mass_pred_s + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)
    b1 ~ dnorm(0, 0.04)
    b2 ~ dnorm(0, 0.04)
    b3 ~ dnorm(0, 0.04)
    sigma_b0 ~ dunif(0, 10)
    alpha.sigma ~ dgamma(0.001, 0.001) # intercept in model of sigma
    b1.sigma ~ dnorm(0, 0.01)
    b2.sigma ~ dnorm(0, 0.01)

  # Derived quantiles
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    
  }", fill = TRUE, file = "R/analysis/JAGS_models/polynomial/polynomial_model1d.txt")


#**** M2 - b0 varies by species, third degree polynomial ===========================
cat(
  "model{
  
  for(i in 1:n_obs){
    
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i] + b4*temp[i]*temp[i]*temp[i]
  
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
  for(k in 1:length(temp_pred)){
    pred_large[k] <- mu_b0 + b1*mass_pred_l + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*temp_pred[k]*temp_pred[k]
    pred_small[k] <- mu_b0 + b1*mass_pred_s + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*temp_pred[k]*temp_pred[k]
  } 

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04)
    b1 ~ dnorm(0, 0.04)
    b2 ~ dnorm(0, 0.04)
    b3 ~ dnorm(0, 0.04)
    b4 ~ dnorm(0, 0.04)
    sigma ~ dunif(0, 10) 
    sigma_b0 ~ dunif(0, 10) 
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
  
  }", fill = TRUE, file = "R/analysis/JAGS_models/polynomial/polynomial_model2.txt")

