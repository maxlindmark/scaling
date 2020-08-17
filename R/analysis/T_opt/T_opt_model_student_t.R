#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of optimum growth temperatures a function of 
# normalized body mass with different group-effects and compare DIC and WAIC
# 
#   Here we use a student-t distribution for the likelihood and compare residuals
#   with the normal model
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit model
# 
# D. Model validation (convergence, fit, residuals)
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

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file
dat <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/topt_analysis.csv"))

# Count data
length(unique(dat$common_name))
nrow(dat)

# Rename species factor for JAGS (must be numbered 1:n)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Predictor variables
# (temperature already centered by subtracting the mean optimum within species in exploratory script)
# (masses are normalized to maturation mass in exploratory script)
dat$log_mass_norm_mat <- log(dat$mass_norm_mat)

# Mean-center predictor variable
dat$log_mass_norm_mat_ct <- dat$log_mass_norm_mat - mean(dat$log_mass_norm_mat)

# Masses for prediction
mass_pred <-  seq(from = min(dat$log_mass_norm_mat_ct), 
                  to = max(dat$log_mass_norm_mat_ct),
                  length.out = 100)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
data = list(
  y = dat$opt_temp_c_ct, 
  n_obs = length(dat$opt_temp_c_ct), 
  species_n = dat$species_n,
  mass = dat$log_mass_norm_mat_ct#,
 #mass_pred = mass_pred
)


# C. FIT MODEL =====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample


#**** M2: Random intercept with a student t likelihood =============================
model2 <- "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dt(mu[i], tau, k)
    mu[i] <- b0[species_n[i]] + b1*mass[i]    # varying intercept and slope

  # Simulate for comparison with data
    y_sim[i] ~ dt(mu[i], tau, k)
    
  # Add log likelihood computation for each observation
    pd[i] <- dt(y[i], mu[i], tau, k)
  
  # Calculates the log PPD
    log_pd[i] <- log(dt(y[i], mu[i], tau, k))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
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
    # shape = 2; rate = 0.1; hist(rgamma(n = 1000, shape = shape, rate = rate, scale = 1/rate))
    k ~ dgamma(2, 0.1)   # this is the degrees of freedom parameter, prior from here: https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/
  
  # Derived quantiles
    tau <- 1/(sigma*sigma)
    tau_b0 <- 1/(sigma_b0*sigma_b0)
  
  }
"
model2 <- textConnection(model2)

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    sigma = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    b1 = 0.001,
    sigma = 0.001,
    sigma_b0 = 0.001,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm2 = jags.model(model2,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm2, n.iter = n.iter) 


# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm2, n.iter = n.iter, thin = thin,
                   variable.names = c("b0", "b1", "mu_b0", "sigma_b0",
                                      "sigma", "k"))

summary(cs)
# 2. Quantiles for each variable:
#             2.5%     25%       50%       75%    97.5%
# ...
# b1       -0.72557 -0.53345 -0.43241 -0.3367 -0.1572
# k         1.83299  4.92901  9.30399 17.3193 44.2271
# mu_b0    -0.38973 -0.06287  0.12387  0.3140  0.7311
# sigma     0.71445  1.17887  1.37781  1.5485  1.8900
# sigma_b0  0.02018  0.16699  0.34289  0.6034  1.3013


#** Evaluate convergence ===========================================================
# # Convert to ggplottable data frame
# cs_df <- ggs(cs)
# 
# #**** Species intercepts ===========================================================
# # Plot posterior densities of species intercepts
# unique(cs_df$Parameter)
# 
# p1 <- cs_df %>% 
#   filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
#                           "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# 
# # Traceplot for evaluating chain convergence
# p2 <- cs_df %>% 
#   filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
#                           "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   NULL
# pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# pWord1 + pWord2
# ggsave("figures/supp/T_opt/validation_topt_intercepts_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Group-level means ============================================================
# # Plot posterior densities of group-level means and standard deviations
# p3 <- cs_df %>% 
#   filter(Parameter %in% c("mu_b0", "mu_b1",
#                           "sigma_b0", "sigma_b1", "sigma", "k")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# 
# # Traceplot for evaluating chain convergence
# p4 <- cs_df %>% 
#   filter(Parameter %in% c("mu_b0", "mu_b1",
#                           "sigma_b0", "sigma_b1", "sigma")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   NULL
# pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# pWord3 + pWord4
# ggsave("figures/supp/T_opt/validation_topt_student.png", width = 6.5, height = 6.5, dpi = 600)
# 
# 
# #**** Chain convergence (Rhat) ====================================================
# p7 <- cs_df %>% 
#   ggs_Rhat(.) + 
#   xlab("R_hat") +
#   xlim(0.999, 1.002) +
#   geom_point(size = 2) +
#   NULL
# pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
#                                        axis.text = element_text(size = 5))
# pWord7
# ggsave("figures/supp/T_opt/validation_rhat_topt_student.png", width = 6.5, height = 6.5, dpi = 600)


#** Evaluate model fit & residuals =================================================
# https://rpubs.com/Niko/332320

# Extract generated data and data
cs_fit = coda.samples(jm2, n.iter = n.iter, thin = thin,
                      variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

# Model fit
p_fit <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  xlim(-1.5, 1.5) + 
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df)))) +
  geom_histogram() +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL

ggsave("figures/supp/T_opt/fit_topt_student_m2.png", width = 6.5, height = 6.5, dpi = 600)

# Residuals
# Extract posteriors for each data point for calculation of residuals
resid <- coda.samples(jm2, variable.names = c("y_sim"), n.iter = n.iter, thin = thin, )

# Extract degrees of freedom parameter from the output
k <- data.frame(summary(cs)[2])
k$param <- rownames(k)
k <- k %>% filter(param == "k")
k <- k$quantiles.50.

# Tidy-up
resid_df <- ggs(resid)

resid_df <- resid_df %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value),
            sd = sd(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = data$y,
         resid = y - yhat,
         index = 1:length(data$y))

# Check linearity
p_lin <- ggplot(resid_df, aes(yhat, resid)) +
  geom_point(shape = 21, fill = "black", color = "white", size = 2) +
  ggtitle("Linearity")

# Check normality
# ** MUST COMPARE TO STUDENT, NOT NORMAL QQ
p_qq <- ggplot(resid_df, aes(sample = resid)) +  # Create QQplot with ggplot2 package
  stat_qq(shape = 21, fill = "black", color = "white", size = 2, distribution = stats::qt, dparams = list(df = k)) +
  stat_qq_line(col = "red", distribution = stats::qt, dparams = list(df = k)) +
  ggtitle("QQ")

p_combo <- (p_lin + p_qq) + plot_annotation(title = 'Growth T_opt model STUDENT')
p_combo & theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)

ggsave("figures/supp/T_opt/resid_topt_student_m2.png", width = 6.5, height = 6.5, dpi = 600)
