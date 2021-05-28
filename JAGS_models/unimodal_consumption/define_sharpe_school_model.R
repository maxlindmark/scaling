#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Define and save models for rjags
# 
# A. Specify models for model selection
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. SPECIFY MODELS ================================================================
#**** Sharpe-Schoolfield with species-varying activation energy and b0

cat(
  "model{

    tref <- -10
    bk <- 8.62e-05

  for(i in 1:n_obs){

  # Likelihood
    y[i] ~ dnorm(mu[i], tau[i])
    mu[i] <- (b0[species_n[i]] * exp(E[species_n[i]]*((1/(bk*(tref + 273.15))) - 1/(bk*(temp[i] + 273.15))))) / (1 + exp(Eh*((1/(bk*(Th + 273.15))) - (1/(bk*(temp[i] + 273.15))))))

    sigma[i] <- alpha.sigma + b1.sigma*temp[i]
    tau[i] <- 1/sigma[i]^2

  # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau[i])
  }

  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    E[j] ~ dnorm(mu_E, tau_E)
  }

  # Predictions
  for(k in 1:length(temp_pred)){
   pred[k] <- (mu_b0 * exp(mu_E*((1/(bk*(tref + 273.15))) - 1/(bk*(temp_pred[k] + 273.15))))) / (1 + exp(Eh*((1/(bk*(Th + 273.15))) - (1/(bk*(temp_pred[k] + 273.15))))))
  }

  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below

  # Priors
    mu_b0 ~ dnorm(1, 0.04)
    mu_E ~ dnorm(0.5, 4)
    #sigma ~ dunif(0, 3)
    sigma_b0 ~ dunif(0, 3)
    sigma_b ~ dunif(0, 3)
    sigma_E ~ dunif(0, 3)
    Eh ~ dnorm(2, 0.25)
    Th ~ dnorm(5, 0.25)
    alpha.sigma ~ dgamma(0.001, 0.001) # intercept in model of sigma
    b1.sigma ~ dnorm(0, 0.01)

  #-- Derived quantiles
    #tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
    tau_E <- 1/(sigma_E*sigma_E)

  }", fill = TRUE, file = "JAGS_models/unimodal_consumption/sharpe_school.txt")
