#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of metabolic and maximum consumption
# rate as a function of temperature and mass
# 
#   Here we use a student-t distribution for the likelihood and compare residuals
#   with the normal model
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit models 
#    - After this section you can run the others separately
#
# D. Model validation (convergence, fit, residuals)
#    - This part is LONG and contains lots of big plots. Skip if you don't want that
# 
# E. Plot predictions
# 
# F. Additional calculations on the posterior
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
library(MCMCvis)
library(scales)

# > sessionInfo()
# other attached packages:
# [1] scales_1.1.0       MCMCvis_0.14.0     bayesplot_1.7.1    patchwork_0.0.1   
# [5] viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [9] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1     
# [13] tidyr_1.0.0       dplyr_0.8.3        RColorBrewer_1.1-2 rjags_4-10        
# [17] coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

# Filter data points at below optimum temperatures
met <- met %>% filter(above_peak_temp == "N")

# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species_ab))

# Mean-center predictor variables
met$log_mass_ct <- met$log_mass - mean(met$log_mass)
met$temp_arr_ct <- met$temp_arr - mean(met$temp_arr)

# Use mass-specific values
met$y_spec <- met$y / met$mass_g

# Masses for prediction
mass_pred_met <-  seq(from = min(met$log_mass_ct), 
                      to = max(met$log_mass_ct),
                      length.out = 100)

# Temperature for prediction
temp_pred_met <- 0 # This means we use mean temperature as it is centered 

# Prepare data for JAGS
met_data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
met_data = list(
  y = log(met$y_spec), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_ct,
  temp = met$temp_arr_ct,
  mass_pred = mass_pred_met,
  temp_pred = temp_pred_met)


# C. FIT MODEL =====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

# Metabolic rate ===================================================================
# Select model with lowest WAIC (see met_model_selection.R)
met_model = "model{
  
  for(i in 1:n_obs){
    
    # Likelihood
    y[i] ~ dt(mu[i], tau, k)            
    mu[i] <- b0[species_n[i]] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i] + b3*mass[i]*temp[i]  
    
    # Simulate for comparison with data (evalute fit)
    y_sim[i] ~ dt(mu[i], tau, k)
    
    # Add log likelihood computation for each observation
    pd[i] <- dt(y[i], mu[i], tau, k)
    
    # Calculates the log PPD
    log_pd[i] <- log(dt(y[i], mu[i], tau, k))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }
  
  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
  mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
  mu_b1 ~ dnorm(-0.25, 1)  
  mu_b2 ~ dnorm(-0.6, 1)   
  b3 ~ dnorm(0, 1)      
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  sigma_b3 ~ dunif(0, 10)
  # shape = 2; rate = 0.1; hist(rgamma(n = 1000, shape = shape, rate = rate, scale = 1/rate))
  k ~ dgamma(2, 0.1)   # this is the degrees of freedom parameter, prior from here: https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/
  
  # Derived quantiles
  tau <- 1/(sigma*sigma)
  tau_b0 <- 1/(sigma_b0*sigma_b0)
  tau_b1 <- 1/(sigma_b1*sigma_b1)
  tau_b2 <- 1/(sigma_b2*sigma_b2)
  
}"

met_model <- textConnection(met_model)


# Manually set initial values, because otherwise all the chains get the same
inits_met = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    b3 = 1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm_met = jags.model(met_model,
                    data = met_data, 
                    n.adapt = 5000, 
                    n.chains = 3,
                    inits = inits_met)

update(jm_met, n.iter = n.iter) 



# D. MODEL VALIDATION ==============================================================
# Metabolic rate ===================================================================
# CODA - Nice for getting the raw posteriors
cs_met <- coda.samples(jm_met,
                       variable.names = c("b0", "b1", "b2", "b3",
                                          "mu_b0", "mu_b1", "mu_b2",
                                          "sigma_b0", "sigma_b1", "sigma_b2",
                                          "sigma"), 
                       n.iter = n.iter, 
                       thin = thin)

# summary(cs_met)
# 2. Quantiles for each variable:
#   
#           2.5%       25%        50%      75%     97.5%
# b3        0.003867  0.010427  0.0139547  0.01757  0.024181
# mu_b0    -2.364561 -2.233227 -2.1672409 -2.10229 -1.969999
# mu_b1    -0.254998 -0.221676 -0.2039408 -0.18716 -0.152763
# mu_b2    -0.666476 -0.629984 -0.6118622 -0.59380 -0.560223


# Convert to ggplottable data frame
cs_met_df <- ggs(cs_met)

#**** Species intercepts (1/2 because to many species) =============================
# Plot posterior densities of species intercepts
unique(cs_met_df$Parameter)

p1a <- cs_met_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",
                          "b0[14]", "b0[15]", "b0[16]", "b0[17]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord1a <- p1a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p2a <- cs_met_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",
                          "b0[14]", "b0[15]", "b0[16]", "b0[17]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord2a <- p2a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord1a + pWord2a
ggsave("figures/supp/log_linear/met_con/validation_met_intercepts_1.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species intercepts (2/2 because to many species) =============================
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1b <- cs_met_df %>% 
  filter(Parameter %in% c("b0[18]", "b0[19]",
                          "b0[20]", "b0[21]", "b0[22]", "b0[23]", "b0[24]", "b0[25]",
                          "b0[25]", "b0[26]", "b0[27]", "b0[28]", "b0[29]", "b0[30]",
                          "b0[31]", "b0[32]", "b0[33]", "b0[34]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord1b <- p1b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p2b <- cs_met_df %>% 
  filter(Parameter %in% c("b0[18]", "b0[19]",
                          "b0[20]", "b0[21]", "b0[22]", "b0[23]", "b0[24]", "b0[25]",
                          "b0[25]", "b0[26]", "b0[27]", "b0[28]", "b0[29]", "b0[30]",
                          "b0[31]", "b0[32]", "b0[33]", "b0[34]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord2b <- p2b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord1b + pWord2b
ggsave("figures/supp/log_linear/met_con/validation_met_intercepts_2.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species mass-effects (1/2 because to many species) ===========================
# Plot posterior densities of species mass-effects
p3a <- cs_met_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b10[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",
                          "b1[14]", "b1[15]", "b1[16]", "b1[17]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord3a <- p3a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))


# Traceplot for evaluating chain convergence
p4a <- cs_met_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b10[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",
                          "b1[14]", "b1[15]", "b1[16]", "b1[17]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord4a <- p4a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord3a + pWord4a
ggsave("figures/supp/log_linear/met_con/validation_met_mass_1.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species mass-effects (2/2 because to many species) ===========================
# Plot posterior densities of species mass-effects
p3b <- cs_met_df %>% 
  filter(Parameter %in% c("b1[18]", "b1[19]",
                          "b1[20]", "b1[21]", "b1[22]", "b1[23]", "b1[24]", "b1[25]",
                          "b1[25]", "b1[26]", "b1[27]", "b1[28]", "b1[29]", "b1[30]",
                          "b1[31]", "b1[32]", "b1[33]", "b1[34]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord3b <- p3b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p4b <- cs_met_df %>% 
  filter(Parameter %in% c("b1[18]", "b1[19]",
                          "b1[20]", "b1[21]", "b1[22]", "b1[23]", "b1[24]", "b1[25]",
                          "b1[25]", "b1[26]", "b1[27]", "b1[28]", "b1[29]", "b1[30]",
                          "b1[31]", "b1[32]", "b1[33]", "b1[34]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord4b <- p4b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord3b + pWord4b
ggsave("figures/supp/log_linear/met_con/validation_met_mass_2.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species temperature-effects (1/2 because to many species) ====================
# Plot posterior densities of temperature-effects
p5a <- cs_met_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]",
                          "b2[14]", "b2[15]", "b2[16]", "b2[17]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord5a <- p5a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p6a <- cs_met_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]",
                          "b2[14]", "b2[15]", "b2[16]", "b2[17]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord6a <- p6a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord5a + pWord6a
ggsave("figures/supp/log_linear/met_con/validation_met_temp_1.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species temperature-effects (2/2 because to many species) ====================
# Plot posterior densities of temperature-effects
p5b <- cs_met_df %>% 
  filter(Parameter %in% c("b2[18]", "b2[19]",
                          "b2[20]", "b2[21]", "b2[22]", "b2[23]", "b2[24]", "b2[25]",
                          "b2[25]", "b2[26]", "b2[27]", "b2[28]", "b2[29]", "b2[30]",
                          "b2[31]", "b2[32]", "b2[33]", "b2[34]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord5b <- p5b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p6b <- cs_met_df %>% 
  filter(Parameter %in% c("b2[18]", "b2[19]",
                          "b2[20]", "b2[21]", "b2[22]", "b2[23]", "b2[24]", "b2[25]",
                          "b2[25]", "b2[26]", "b2[27]", "b2[28]", "b2[29]", "b2[30]",
                          "b2[31]", "b2[32]", "b2[33]", "b2[34]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord6b <- p6b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord5b + pWord6b
ggsave("figures/supp/log_linear/met_con/validation_met_temp_2.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities of group-level means and standard deviations
p7 <- cs_met_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "b3",
                          "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p8 <- cs_met_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "b3",
                          "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord8 <- p8 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord7 + pWord8
ggsave("figures/supp/log_linear/met_con/validation_met.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p9 <- cs_met_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.004) +
  geom_point(size = 1) +
  NULL
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
ggsave("figures/supp/log_linear/met_con/validation_rhat_met.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(0, 0.04)   # mean of all species intercept
# mu_b1 ~ dnorm(-0.25, 1)  # mean of all species mass-exponent
# mu_b2 ~ dnorm(-0.6, 1)   # mean of all species temperature coefficient
# b3 ~ dnorm(0, 1)         # global interaction

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# tau <- 1/sigma^2
# from tau to sigma 
# sigma <- sqrt(1/tau)

# Define priors for plot
tau <- 1
tau_int <- 0.04

mu_b0 <- rnorm(25000, 0, sqrt(1/tau_int))  
mu_b1 <- rnorm(25000, -0.25, sqrt(1/tau))  
mu_b2 <- rnorm(25000, -0.6, sqrt(1/tau))   
b3 <- rnorm(25000, 0, sqrt(1/tau))         

PR <- as.matrix(cbind(mu_b0, mu_b1, mu_b2, b3))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/log_linear/met_con/validation_prior_post_met.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs_met,
          params = c("mu_b0", "mu_b1", "mu_b2", "b3"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE,
          type = "density")   # removes the trace plot

dev.off()



# Evaluate model fit & residuals ===================================================
# https://rpubs.com/Niko/332320

#**** Metabolism ===================================================================
# Extract generated data and data
cs_fit_met = coda.samples(jm_met, n.iter = n.iter, thin = thin,
                          variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df_met <- data.frame(as.matrix(cs_fit_met))

# Model fit
p_fit_m <- ggplot(cs_fit_df_met, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df_met)))) +
  geom_histogram() +
  geom_vline(xintercept = cs_fit_df_met$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL

ggsave("figures/supp/log_linear/met_con/fit_met_mean.png", width = 6.5, height = 6.5, dpi = 600)

# Residuals
# Extract posteriors for each data point for calculation of residuals
resid <- coda.samples(jm_met, variable.names = c("y_sim"), n.iter = n.iter, thin = thin, )

# Tidy-up
resid_df <- ggs(resid)

resid_df <- resid_df %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = met_data$y,
         resid = y - yhat)

# Check linearity
p_lin <- ggplot(resid_df, aes(yhat, resid)) +
  geom_point() +
  ggtitle("Linearity")

# Check normality
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle("QQ")

p_combo <- (p_lin + p_qq) + plot_annotation(title = 'Log-linear metabolism')
p_combo & theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)

ggsave("figures/supp/log_linear/met_con/resid_met.png", width = 6.5, height = 6.5, dpi = 600)


#**** Consumption ===================================================================
# Extract generated data and data
cs_fit_con = coda.samples(jm_con, n.iter = n.iter, thin = thin,
                          variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))

# Model fit
p_fit_c <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df_con)))) +
  geom_histogram() +
  geom_vline(xintercept = cs_fit_df_con$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL

ggsave("figures/supp/log_linear/met_con/fit_con_mean.png", width = 6.5, height = 6.5, dpi = 600)

# Residuals
# Extract posteriors for each data point for calculation of residuals
resid <- coda.samples(jm_con, variable.names = c("y_sim"), n.iter = n.iter, thin = thin, )

# Tidy-up
resid_df <- ggs(resid)

resid_df <- resid_df %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = con_data$y,
         resid = y - yhat)

# Check linearity
p_lin <- ggplot(resid_df, aes(yhat, resid)) +
  geom_point() +
  ggtitle("Linearity")

# Check normality
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle("QQ")

p_combo <- (p_lin + p_qq) + plot_annotation(title = 'Log-linear consumption')
p_combo & theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)

ggsave("figures/supp/log_linear/met_con/resid_con.png", width = 6.5, height = 6.5, dpi = 600)

