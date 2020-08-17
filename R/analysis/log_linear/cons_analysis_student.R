#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of maximum consumption rate as a function of temperature and mass
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
con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Filter data points at below optimum temperatures
con <- con %>% filter(above_peak_temp == "N")

# Rename species factor for JAGS (must be numbered 1:n)
con$species_n <- as.numeric(as.factor(con$species_ab))

# Mean-center predictor variables
con$log_mass_ct <- con$log_mass - mean(con$log_mass)
con$temp_arr_ct <- con$temp_arr - mean(con$temp_arr)

# Use mass-specific values
con$y_spec <- con$y / con$mass_g

# Prepare data for JAGS
con_data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
con_data = list(
  y = log(con$y_spec), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_ct,
  temp = con$temp_arr_ct)


# C. FIT MODEL =====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

# Maximum consumption rate =========================================================
# Select model with lowest WAIC (see con_model_selection.R)
con_model = "model{
  
  for(i in 1:n_obs){
    
    # Likelihood
    y[i] ~ dt(mu[i], tau, k)            
    mu[i] <- b0[species_n[i]] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i]
    
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
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  # shape = 2; rate = 0.1; hist(rgamma(n = 1000, shape = shape, rate = rate, scale = 1/rate))
  k ~ dgamma(2, 0.1)   # this is the degrees of freedom parameter, prior from here: https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/
  
  # Derived quantiles
  tau <- 1/(sigma*sigma)
  tau_b0 <- 1/(sigma_b0*sigma_b0)
  tau_b1 <- 1/(sigma_b1*sigma_b1)
  tau_b2 <- 1/(sigma_b2*sigma_b2)
  
}"

con_model <- textConnection(con_model)


# Manually set initial values, because otherwise all the chains get the same
inits_con = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
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
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm_con = jags.model(con_model,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3,
                    inits = inits_con)

update(jm_con, n.iter = n.iter) 



# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
# cs_con <- coda.samples(jm_con,
#                        variable.names = c("b0", "b1", "b2", 
#                                           "mu_b0", "mu_b1", "mu_b2",
#                                           "sigma_b0", "sigma_b1", "sigma_b2",
#                                           "sigma"), 
#                        n.iter = n.iter, 
#                        thin = thin)
# 
# # summary(cs_con)
# # 2. Quantiles for each variable:
# #   
# #          2.5%     25%     50%     75%    97.5%
# # mu_b0    -3.46443 -3.1237 -2.9576 -2.7944 -2.44343
# # mu_b1    -0.45399 -0.3997 -0.3750 -0.3489 -0.29574
# # mu_b2    -0.84680 -0.7429 -0.6945 -0.6430 -0.54058
# 
# # Convert to ggplottable data frame
# cs_con_df <- ggs(cs_con)
# 
# #**** Species intercepts ===========================================================
# # Plot posterior densities of species intercepts
# unique(cs_con_df$Parameter)
# 
# p1 <- cs_con_df %>% 
#   filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
#                           "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",
#                           "b0[14]", "b0[15]", "b0[16]", "b0[17]", "b0[18]", "b0[19]", "b0[20]")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 8),
#                                        axis.text = element_text(size = 5))
# 
# # Traceplot for evaluating chain convergence
# p2 <- cs_con_df %>% 
#   filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
#                           "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",
#                           "b0[14]", "b0[15]", "b0[16]", "b0[17]", "b0[18]", "b0[19]", "b0[20]")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   NULL
# pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 8),
#                                        axis.text = element_text(size = 5))
# pWord1 + pWord2
# ggsave("figures/supp/log_linear/met_con/validation_con_intercepts_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Species mass-effects =========================================================
# # Plot posterior densities of species mass-effects
# p3 <- cs_con_df %>% 
#   filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
#                           "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",
#                           "b1[14]", "b1[15]", "b1[16]", "b1[17]", "b1[18]", "b1[19]", "b0[20]")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 8),
#                                        axis.text = element_text(size = 5))
# 
# 
# # Traceplot for evaluating chain convergence
# p4 <- cs_con_df %>% 
#   filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
#                           "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",
#                           "b1[14]", "b1[15]", "b1[16]", "b1[17]", "b1[18]", "b1[19]", "b0[20]")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   NULL
# pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 8),
#                                        axis.text = element_text(size = 5))
# pWord3 + pWord4
# ggsave("figures/supp/log_linear/met_con/validation_con_mass_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Species temperature-effects ==================================================
# # Plot posterior densities of temperature-effects
# p5 <- cs_con_df %>% 
#   filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
#                           "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]",
#                           "b2[14]", "b2[15]", "b2[16]", "b2[17]", "b2[18]", "b2[19]", "b0[20]")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 8),
#                                        axis.text = element_text(size = 5))
# 
# # Traceplot for evaluating chain convergence
# p6 <- cs_con_df %>% 
#   filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
#                           "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]",
#                           "b2[14]", "b2[15]", "b2[16]", "b2[17]", "b2[18]", "b2[19]", "b0[20]")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   NULL
# pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 8),
#                                        axis.text = element_text(size = 5))
# pWord5 + pWord6
# ggsave("figures/supp/log_linear/met_con/validation_con_temp_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Group-level means ============================================================
# # Plot posterior densities of group-level means and standard deviations
# p7 <- cs_con_df %>% 
#   filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", 
#                           "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# 
# # Traceplot for evaluating chain convergence
# p8 <- cs_con_df %>% 
#   filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", 
#                           "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   NULL
# pWord8 <- p8 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# pWord7 + pWord8
# ggsave("figures/supp/log_linear/met_con/validation_con_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Chain convergencve (Rhat) ====================================================
# p9 <- cs_con_df %>% 
#   ggs_Rhat(.) + 
#   xlab("R_hat") +
#   xlim(0.999, 1.003) +
#   geom_point(size = 1) +
#   NULL
# pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# ggsave("figures/supp/log_linear/met_con/validation_rhat_con_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Prior vs posterior ===========================================================
# # https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html
# 
# # Priors from JAGS
# # mu_b0 ~ dnorm(0, 0.04)   # mean of all species intercept
# # mu_b1 ~ dnorm(-0.25, 1)  # mean of all species mass-exponent
# # mu_b2 ~ dnorm(-0.6, 1)   # mean of all species temperature coefficient
# 
# # Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# # tau = 1/variance
# # from sigma to tau 
# # sigma = 1
# # tau <- 1/sigma^2
# 
# # Define priors for plot
# tau <- 1
# tau_int <- 0.04
# 
# mu_b0 <- rnorm(25000, 0, sqrt(1/tau_int))
# mu_b1 <- rnorm(25000, -0.25, sqrt(1/tau))
# mu_b2 <- rnorm(25000, -0.6, sqrt(1/tau)) 
# 
# PR <- as.matrix(cbind(mu_b0, mu_b1, mu_b2))
# 
# # This is not a ggplot...
# png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/log_linear/met_con/validation_prior_post_con_student.png", 
#     units = "px", width = 1800, height = 1800, res = 300)
# 
# MCMCtrace(cs_con,
#           params = c("mu_b0", "mu_b1", "mu_b2"),
#           ISB = FALSE,
#           priors = PR,
#           pdf = FALSE,
#           Rhat = TRUE,
#           n.eff = TRUE,
#           type = "density")   # removes the trace plot
# 
# dev.off()



# Evaluate model fit & residuals ===================================================
# https://rpubs.com/Niko/332320

#**** Consumption ==================================================================
# Extract generated data and data
# cs_fit_con = coda.samples(jm_con, n.iter = n.iter, thin = thin,
#                           variable.names = c("mean_y", "mean_y_sim", "p_mean"))
# 
# # Convert to data frames
# cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))
# 
# # Model fit
# p_fit_c <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
#   coord_cartesian(expand = 0) +
#   geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df_con)))) +
#   geom_histogram() +
#   geom_vline(xintercept = cs_fit_df_con$mean_y, color = "white", 
#              linetype = 2, size = 0.4) +
#   labs(x = "mean simulated data", y = "count") +
#   theme_classic() +
#   theme(text = element_text(size = 12), aspect.ratio = 1) +
#   NULL
# 
# ggsave("figures/supp/log_linear/met_con/fit_con_mean.png", width = 6.5, height = 6.5, dpi = 600)

# Residuals
# Extract posteriors for each data point for calculation of residuals
resid <- coda.samples(jm_con, variable.names = c("y_sim", "k"), n.iter = n.iter, thin = thin, )

# Tidy-up
resid_df <- ggs(resid)

# Extract degrees of freedom parameter:
k <- resid_df %>% 
  filter(Parameter == "k") %>% 
  summarize(median = median(value)) %>% 
  as.vector()

# Get residuals
resid_df <- resid_df %>%
  ungroup() %>%
  group_by(Parameter) %>%
  filter(!Parameter == "k") %>% 
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
  stat_qq(distribution = stats::qt, dparams = list(df = k$median)) +
  stat_qq_line(distribution = stats::qt, dparams = list(df = k$median)) +
  ggtitle("QQ")

p_combo <- (p_lin + p_qq) + plot_annotation(title = 'Log-linear metabolism')
p_combo & theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)

ggsave("figures/supp/log_linear/met_con/resid_con_student.png", width = 6.5, height = 6.5, dpi = 600)
