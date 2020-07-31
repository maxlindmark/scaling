#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of growth rate as a function of 
# temperature and mass
# 
#   Here we use a student-t distribution for the likelihood and compare residuals
#   with the normal model
#
# A. Load libraries
#
# B. Read data
#
# C. Fit models
#
# D. Model validation (convergence, fit, residuals)
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

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file
dat <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

str(dat)

# Count data
length(unique(dat$common_name))
nrow(dat)

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_peak_temp == "N")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Mean-center predictor variables
dat$log_mass_ct <- dat$log_mass - mean(dat$log_mass)
dat$temp_arr_ct <- dat$temp_arr - mean(dat$temp_arr)

# Use only positive growth rates
dat <- dat %>% filter(y > 0)

# Masses for prediction
mass_pred <-  seq(from = min(dat$log_mass_ct), 
                  to = max(dat$log_mass_ct),
                  length.out = 100)

# Temperature for prediction
temp_pred <- 0 # This means we use mean temperature as it is centered 

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
data = list(
  y = log(dat$y), 
  n_obs = length(dat$y), 
  species_n = dat$species_n,
  mass = dat$log_mass_ct,
  temp = dat$temp_arr_ct#,
  # mass_pred = mass_pred,
  # temp_pred = temp_pred
)

mean(dat$temp_c)

# C. FIT MODEL =====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

# Select model with lowest WAIC (see grow_model_selection.R)
model1 = "model{
  
  for(i in 1:n_obs){
  
  # Likelihood
    y[i] ~ dt(mu[i], tau, k)            
    mu[i] <- b0[species_n[i]] + b1[species_n[i]]*mass[i] + b2[species_n[i]]*temp[i] + b3[species_n[i]]*mass[i]*temp[i]  
    
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
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Model fit
    mean_y <- mean(y[])
    mean_y_sim <- mean(y_sim[])
    p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  # Priors	
    mu_b0 ~ dnorm(0, 0.04) # remember the second argument is precision (1/variance)   
    mu_b1 ~ dnorm(-0.25, 1)  
    mu_b2 ~ dnorm(-0.6, 1)   
    mu_b3 ~ dnorm(0, 1)      
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
    tau_b3 <- 1/(sigma_b3*sigma_b3)
  
}"

model1 <- textConnection(model1)
  
# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    mu_b3 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    sigma_b3 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    mu_b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    sigma_b3 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    mu_b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    sigma_b3 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm = jags.model(model1,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3,
                inits = inits)

update(jm, n.iter = n.iter) 


# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm,
                   variable.names = c("b0", "b1", "b2", "b3", 
                                      "mu_b0", "mu_b1", "mu_b2", "mu_b3",
                                      "sigma_b0", "sigma_b1", "sigma_b2", "sigma_b3",
                                      "sigma", "k"), 
                   n.iter = n.iter, 
                   thin = thin)

summary(cs)
# 2. Quantiles for each variable:
#   
#           2.5%      25%       50%      75%      97.5%
# ...
# k         1.813826  2.582085  3.1553275  3.956740  6.77335
# mu_b0     0.181386  0.483368  0.6119490  0.739065  1.01114
# mu_b1    -0.463698 -0.380662 -0.3437053 -0.307054 -0.22326
# mu_b2    -0.778933 -0.672176 -0.6257087 -0.582567 -0.49957
# mu_b3    -0.043146 -0.004392  0.0137486  0.031135  0.07109
# sigma     0.123137  0.147122  0.1611916  0.176357  0.20787
# sigma_b0  0.400038  0.546277  0.6521094  0.783043  1.14734
# sigma_b1  0.074714  0.130433  0.1684221  0.214368  0.34491
# sigma_b2  0.015392  0.085574  0.1317332  0.183879  0.32011
# sigma_b3  0.003243  0.024151  0.0441224  0.068236  0.12722

#** Evaluate convergence ===========================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

#**** Species intercepts ===========================================================
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p2 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord1 + pWord2
ggsave("figures/supp/log_linear/growth/validation_gro_intercepts_student.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species mass-effects =========================================================
# Plot posterior densities of species mass-effects
p3 <- cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))


# Traceplot for evaluating chain convergence
p4 <- cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord3 + pWord4
ggsave("figures/supp/log_linear/growth/validation_gro_mass_student.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species temperature-effects ==================================================
# Plot posterior densities of temperature-effects
p5 <- cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p6 <- cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord5 + pWord6
ggsave("figures/supp/log_linear/growth/validation_gro_temp_student.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities of group-level means and standard deviations
p7 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "mu_b3",
                          "sigma_b0", "sigma_b1", "sigma_b2", "sigma_b3")) %>% 
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
p8 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "mu_b3",
                          "sigma_b0", "sigma_b1", "sigma_b2", "sigma_b3")) %>% 
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
ggsave("figures/supp/log_linear/growth/validation_gro_student.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p9 <- cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.003) +
  geom_point(size = 2) +
  NULL
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
ggsave("figures/supp/log_linear/growth/validation_rhat_gro_student.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(0, 0.04)   # mean of all species intercept
# mu_b1 ~ dnorm(-0.25, 1)  # mean of all species mass-exponent
# mu_b2 ~ dnorm(-0.6, 1)   # mean of all species temperature coefficient
# mu_b3 ~ dnorm(0, 1)      # mean of all species interaction

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# sigma = 1
# tau <- 1/sigma^2

# Define priors for plot
tau <- 1
tau_int <- 0.04

mu_b0 <- rnorm(25000, 0, sqrt(1/tau_int))
mu_b1 <- rnorm(25000, -0.25, sqrt(1/tau))
mu_b2 <- rnorm(25000, -0.6, sqrt(1/tau))
mu_b3 <- rnorm(25000, 0, sqrt(1/tau))

PR <- as.matrix(cbind(mu_b0, mu_b1, mu_b2, mu_b3))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/log_linear/growth/validation_prior_post_growth_student.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs,
          params = c("mu_b0", "mu_b1", "mu_b2", "mu_b3"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE,
          type = "density")   # removes the trace plot

dev.off()


#** Evaluate model fit & residuals =================================================
# https://rpubs.com/Niko/332320

# Extract generated data and data
cs_fit = coda.samples(jm, n.iter = n.iter, thin = thin,
                      variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

# Model fit
p_fit <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df)))) +
  geom_histogram() +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL

ggsave("figures/supp/log_linear/growth/fit_gro_mean_student.png", width = 6.5, height = 6.5, dpi = 600)


# Residuals
# Extract posteriors for each data point for calculation of residuals
resid <- coda.samples(jm, variable.names = c("y_sim"), n.iter = n.iter, thin = thin, )

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
         resid = y - yhat)

# Check linearity
p_lin <- ggplot(resid_df, aes(yhat, resid)) +
  geom_point() +
  ggtitle("Linearity")

# Check normality
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq(distribution = stats::qt, dparams = list(df = k)) +
  stat_qq_line(distribution = stats::qt, dparams = list(df = k)) +
  ggtitle("QQ")

p_combo <- (p_lin + p_qq) + plot_annotation(title = 'Log-linear growth STUDENT')
p_combo & theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)

ggsave("figures/supp/log_linear/growth/resid_growth_student.png", width = 6.5, height = 6.5, dpi = 600)
