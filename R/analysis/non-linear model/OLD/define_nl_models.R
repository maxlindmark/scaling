#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define and save models for rjags for model selection using WAIC
# 
# A. Specify models non-linear models
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. SPECIFY MODELS ================================================================
# Here we fit models with different hierarcial structures
# Specifically, we consider:

# -- Without interaction
# M1 - b0 
# M2 - b1
# M3 - b2 + b3 # If b3 random then b2 - both temperature coefficients 

# -- With interaction
# M4 - b0 
# M5 - b1
# M6 - b2 + b3 # If b3 random then b2 - both temperature coefficients 
# M7 - b4


#** TEST ===========================================================================
# cat(
#   "model{
#   
#   for(i in 1:n_obs){
#     # Simulate for comparison with data
#     y_sim[i] ~ dnorm(mu[i], tau)
#     
#     y[i] ~ dnorm(mu[i], tau)
#     
#     mu[i] <- b0 + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
#   
#     # Add log likelihood computation for each observation
#     pd[i] <- dnorm(y[i], mu[i], tau)
#   
#     # Calculates the log PPD
#     log_pd[i] <- log(dnorm(y[i], mu[i], tau))
#   }
#   
#   # Predictions
#   for(k in 1:length(temp_pred)){
#       
#     pred_large[k] <- b0 + b1*4 + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
#     pred_medium[k] <- b0 + b1*0 + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
#     pred_small[k] <- b0 + b1*-4 + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
#     
#   } 
# 
#   # Model fit
#   mean_y <- mean(y[])
#   mean_y_sim <- mean(y_sim[])
#   p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
#   
#   cv_y <- sd(y[])/mean(y[])
#   cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
#   p_cv <- step(cv_y_sim - cv_y)
# 
#   #-- Priors	
#   b0 ~ dnorm(0, 5)
#   b1 ~ dnorm(0, 5)
#   b2 ~ dnorm(0, 5)
#   b3 ~ dnorm(0, 5)
#   sigma ~ dunif(0, 10) 
#   tau <- 1/sigma^2
#   
#   }", fill = TRUE, file = "R/analysis/non-linear model/models/test.txt")


#** No Mass-Temperature Interaction ================================================
#**** M1 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
  
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
      
    pred_large[k] <- mu_b0 + b1*large + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_medium[k] <- mu_b0 + b1*medium + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- mu_b0 + b1*small + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b1 ~ dnorm(0, 5)
  b2 ~ dnorm(0, 5)
  b3 ~ dnorm(0, 5)
  mu_b0 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2

  }", fill = TRUE, file = "R/analysis/non-linear model/models/m1.txt")


#**** M2 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[species_n[i]]*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred_large[k] <- b0 + mu_b1*large + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_medium[k] <- b0 + mu_b1*medium + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- b0 + mu_b1*small + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b0 ~ dnorm(0, 5)
  b2 ~ dnorm(0, 5)
  b3 ~ dnorm(0, 5)
  mu_b1 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b1 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b1 <- 1/sigma_b1^2
  
  }", fill = TRUE, file = "R/analysis/non-linear model/models/m2.txt")


#**** M3 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1*mass[i] + b2[species_n[i]]*temp[i] + b3[species_n[i]]*temp[i]*temp[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b2[j] ~ dnorm(mu_b2, tau_b2)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred_large[k] <- b0 + b1*large + mu_b2*temp_pred[k] + mu_b3*temp_pred[k]*temp_pred[k]
    pred_medium[k] <- b0 + b1*medium + mu_b2*temp_pred[k] + mu_b3*temp_pred[k]*temp_pred[k]
    pred_small[k] <- b0 + b1*small + mu_b2*temp_pred[k] + mu_b3*temp_pred[k]*temp_pred[k]
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b0 ~ dnorm(0, 5)
  b1 ~ dnorm(0, 5)
  mu_b2 ~ dnorm(0, 5)
  mu_b3 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b2 ~ dunif(0, 10)
  sigma_b3 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b2 <- 1/sigma_b2^2
  tau_b3 <- 1/sigma_b3^2
  
  }", fill = TRUE, file = "R/analysis/non-linear model/models/m3.txt")


#** Mass-Temperature Interaction ===================================================
#**** M4 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0[species_n[i]] + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i] + b4*temp[i]*mass[i]
  
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
      
    pred_large[k] <- mu_b0 + b1*large + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*4
    pred_medium[k] <- mu_b0 + b1*medium + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*0
    pred_small[k] <- mu_b0 + b1*small + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*-4
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b1 ~ dnorm(0, 5)
  b2 ~ dnorm(0, 5)
  b3 ~ dnorm(0, 5)
  b4 ~ dnorm(0, 5)
  mu_b0 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2

  }", fill = TRUE, file = "R/analysis/non-linear model/models/m4.txt")


#**** M5 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[species_n[i]]*mass[i] + b2*temp[i] + b3*temp[i]*temp[i] + b4*temp[i]*mass[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred_large[k] <- b0 + mu_b1*large + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*large
    pred_medium[k] <- b0 + mu_b1*medium + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*medium
    pred_small[k] <- b0 + mu_b1*small + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*small
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b0 ~ dnorm(0, 5)
  b2 ~ dnorm(0, 5)
  b3 ~ dnorm(0, 5)
  b4 ~ dnorm(0, 5)
  mu_b1 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b1 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b1 <- 1/sigma_b1^2
  
  }", fill = TRUE, file = "R/analysis/non-linear model/models/m5.txt")


#**** M6 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1*mass[i] + b2[species_n[i]]*temp[i] + b3[species_n[i]]*temp[i]*temp[i] + b4*temp[i]*mass[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b2[j] ~ dnorm(mu_b2, tau_b2)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred_large[k] <- b0 + b1*large + mu_b2*temp_pred[k] + mu_b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*large
    pred_medium[k] <- b0 + b1*medium + mu_b2*temp_pred[k] + mu_b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*medium
    pred_small[k] <- b0 + b1*small + mu_b2*temp_pred[k] + mu_b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*small
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b0 ~ dnorm(0, 5)
  b1 ~ dnorm(0, 5)
  mu_b2 ~ dnorm(0, 5)
  mu_b3 ~ dnorm(0, 5)
  b4 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b2 ~ dunif(0, 10)
  sigma_b3 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b2 <- 1/sigma_b2^2
  tau_b3 <- 1/sigma_b3^2
  
  }", fill = TRUE, file = "R/analysis/non-linear model/models/m6.txt")


#**** M7 =========================================================================== 
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i] + b4[species_n[i]]*temp[i]*mass[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b4[j] ~ dnorm(mu_b4, tau_b4)
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred_large[k] <- b0 + b1*large + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + mu_b4*temp_pred[k]*large
    pred_medium[k] <- b0 + b1*medium + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + mu_b4*temp_pred[k]*medium
    pred_small[k] <- b0 + b1*small + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + mu_b4*temp_pred[k]*small
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b0 ~ dnorm(0, 5)
  b1 ~ dnorm(0, 5)
  b2 ~ dnorm(0, 5)
  b3 ~ dnorm(0, 5)
  mu_b4 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  sigma_b4 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b4 <- 1/sigma_b4^2
  
  }", fill = TRUE, file = "R/analysis/non-linear model/models/m7.txt")
