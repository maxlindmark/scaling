#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.06: Max Lindmark
#
# Code to fit and evaluate non-linear (polynomial) models to consumption rates
# 
# A. Load libraries
#
# B. Read data
#
# C. Model selection
#
# D. Model validation (convergence, fit, residuals)
# 
# E. Plot predictions
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
con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Which species have data above optimum?
spec <- unique(filter(con, above_peak_temp == "Y"))$species

# length(unique(spec))

con <- con %>% filter(species %in% spec) %>% droplevels()

# Rename species factor for JAGS (must be numbered 1:n)
str(con)
con$species_n <- as.numeric(as.factor(con$species))

# Print numeric factor level vs species name
sort(unique(con$species_ab))
unique(con$species_n)

# Center temperature relative to mean in environment by species
con$temp_env_ct <- con$temp_c - con$median_temp

# Mean-center mass
con$mass_g_ct <- con$mass_g - mean(con$mass_g)

# Add mass-specific consumption
con$y_spec <- con$y / con$mass_g

# Express consumption rate as fraction of mean within species
con$y_mean_species <- ave(con$y_spec, con$species_ab)

con$y_ct <- con$y_spec / con$y_mean_species

# Temp-range used for prediction
temp_pred = seq(from = min(con$temp_env_ct), 
                to = max(con$temp_env_ct),
                length.out = 100)

summary(con$mass_g_ct)

# Mass-range for prediction (note this is on centered scale)
summary(con$mass_g)
mass_pred_s = -118
mass_pred_l = 77 

# Plot to compare with real scale
ggplot(con, aes(mass_g_ct, mass_g)) +
  geom_point() +
  geom_hline(yintercept = 200) +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = -118) +
  geom_vline(xintercept = 77) +
  NULL

#con %>% filter(mass_g_ct > -117 & mass_g_ct < -113)  

# Inspect data final time before fitting (this was spotted when looking the residuals)
ggplot(con, aes(mass_g_ct, y_ct)) + geom_point()

# Remove outlier
con <- con %>% filter(y_ct < 4)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
data = list(
  y = con$y_ct, 
  n_obs = length(con$y_ct), 
  species_n = con$species_n,
  temp = con$temp_env_ct,
  temp_pred = temp_pred,
  mass = con$mass_g_ct,
  mass_pred_s = mass_pred_s,
  mass_pred_l = mass_pred_l)

#length(unique(con$species))

# C. MODEL SELECTION ===============================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample


#**** M1: 2nd degree polynomial ====================================================
model1 = "R/analysis/JAGS_models/non_linear_model/non_linear_model1.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    b1 = 0.001,
    b2 = 0.001,
    b3 = 0.001,
    sigma = 0.001,
    sigma_b0 = 0.001,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm1 = jags.model(model1,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm1, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj1 = jags.samples(jm1, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd1 <- -2*sum(log(summary(zj1$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC1 <- sum((summary(zj1$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m1 <- lppd1 + 2*pd.WAIC1


#**** M1b: 2nd degree polynomial with sigma increasing with mass ===================
model1b = "R/analysis/JAGS_models/non_linear_model/non_linear_model1b.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    b1 = 0.001,
    b2 = 0.001,
    b3 = 0.001,
    sigma_b0 = 0.001,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm1b = jags.model(model1b,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm1b, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj1b = jags.samples(jm1b, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd1b <- -2*sum(log(summary(zj1b$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC1b <- sum((summary(zj1b$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m1b <- lppd1b + 2*pd.WAIC1b


#**** M1c: 2nd degree polynomial with sigma increasing with temp ===================
model1c = "R/analysis/JAGS_models/non_linear_model/non_linear_model1c.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    b1 = 0.001,
    b2 = 0.001,
    b3 = 0.001,
    sigma_b0 = 0.001,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm1c = jags.model(model1c,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm1c, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj1c = jags.samples(jm1c, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd1c <- -2*sum(log(summary(zj1c$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC1c <- sum((summary(zj1c$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m1c <- lppd1c + 2*pd.WAIC1c


#**** M1d: 2nd degree polynomial with sigma increasing with mass & temp ============
model1d = "R/analysis/JAGS_models/non_linear_model/non_linear_model1d.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    b1 = 0.001,
    b2 = 0.001,
    b3 = 0.001,
    sigma_b0 = 0.001,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm1d = jags.model(model1d,
                  data = data, 
                  n.adapt = 5000, 
                  n.chains = 3,
                  inits = inits)

update(jm1d, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj1d = jags.samples(jm1d, 
                    variable.names = c("pd", "log_pd"), 
                    n.iter = n.iter, 
                    thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd1d <- -2*sum(log(summary(zj1d$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC1d <- sum((summary(zj1d$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m1d <- lppd1d + 2*pd.WAIC1d


#**** M2: 3rd degree polynomial ====================================================
model2 = "R/analysis/JAGS_models/non_linear_model/non_linear_model2.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    b1 = 0.1,
    b2 = 0.1,
    b3 = 0.1,
    b4 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    b1 = 1,
    b2 = 1,
    b3 = 1,
    b4 = 1,
    sigma = 1,
    sigma_b0 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    b1 = 0.001,
    b2 = 0.001,
    b3 = 0.001,
    b4 = 0.001,
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

# Monitor the likelihood to calculate WAIC
zj2 = jags.samples(jm2, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = n.iter, 
                   thin = thin)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd2 <- -2*sum(log(summary(zj2$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC2 <- sum((summary(zj2$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
waic_m2 <- lppd2 + 2*pd.WAIC2


#**** Compare WAIC =================================================================
waic_m1
waic_m1b
waic_m1c
waic_m1d
waic_m2

# > waic_m1
# [1] 619.5119
# > waic_m1b
# [1] 593.4074
# > waic_m1c
# [1] 588.1218
# > waic_m1d
# [1] 1829.666
# > waic_m2
# [1] 615.7727

waic_m1-waic_m1c
waic_m1b-waic_m1c
waic_m1c-waic_m1c
waic_m1d-waic_m1c
waic_m2-waic_m1c


# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors

cs <- coda.samples(jm1c, n.iter = n.iter, thin = thin,
                   variable.names = c("mu_b0", "b0", "b1", "b2", "b3",
                                      "sigma_b0", "alpha.sigma", "b1.sigma"))

summary(cs)
# 2. Quantiles for each variable:
#   
#              2.5%       25%        50%        75%        97.5%
# alpha.sigma  0.369307  0.3858199  0.3951141  4.049e-01  0.4257428
# ...
# b1          -0.002087 -0.0017880 -0.0016305 -1.479e-03 -0.0011969
# b1.sigma     0.008750  0.0111996  0.0124513  1.360e-02  0.0156307
# b2           0.036207  0.0399303  0.0419299  4.393e-02  0.0479142
# b3          -0.000756 -0.0004276 -0.0002519 -8.267e-05  0.0002351
# mu_b0        0.694740  0.7750862  0.8097589  8.454e-01  0.9266065
# sigma_b0     0.060548  0.1142524  0.1483631  1.887e-01  0.3053840


#** Evaluate convergence ===========================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

#**** Species intercepts (b0) ======================================================
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]")) %>% 
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
                          "b0[8]", "b0[9]", "b0[10]")) %>% 
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
ggsave("figures/supp/non_linear/validation_non_linear_intercepts_sigma.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities of group-level means and standard deviations
p3 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "b1", "b2", "b3",
                          "sigma_b0", "alpha.sigma", "b1.sigma")) %>% 
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
  filter(Parameter %in% c("mu_b0", "b1", "b2", "b3", 
                          "sigma_b0", "alpha.sigma", "b1.sigma")) %>% 
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
ggsave("figures/supp/non_linear/validation_non_linear_sigma.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p5 <- cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.001) +
  geom_point(size = 2) +
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

ggsave("figures/supp/non_linear/validation_rhat_non_linear_sigma.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(0, 0.04)
# b1 ~ dnorm(0, 0.04)
# b2 ~ dnorm(0, 0.04)
# b3 ~ dnorm(0, 0.04)

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# sigma = 1
# tau <- 1/sigma^2

# Define priors for plot
tau <- 0.04

mu_b0 <- rnorm(25000, 0, sqrt(1/tau))
b1 <- rnorm(25000, 0, sqrt(1/tau))
b2 <- rnorm(25000, 0, sqrt(1/tau))
b3 <- rnorm(25000, 0, sqrt(1/tau))

PR <- as.matrix(cbind(mu_b0, b1, b2, b3))

PR <- as.matrix(cbind(mu_b0, b2, b3))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/non_linear/validation_prior_post_non_linear_sigma.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs,
          params = c("mu_b0", "b2", "b3"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE,
          type = "density")   # removes the trace plot

dev.off()

# Some bug here, can't use b1????


#** Evaluate model fit & residuals =================================================
# https://rpubs.com/Niko/332320

# Extract generated data and data
cs_fit = coda.samples(jm1c, n.iter = n.iter, thin = thin,
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

ggsave("figures/supp/non_linear/fit_nonlinear_sigma.png", width = 6.5, height = 6.5, dpi = 600)


# Residuals
# Extract posteriors for each data point for calculation of residuals
resid <- coda.samples(jm1c, variable.names = c("y_sim"), n.iter = n.iter, thin = thin, )

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

# Now we in addition need to calculate the standardized residuals, else we cant detect any improvements
# I simply calculate them as the residual divided by the standard deviation. https://www.isixsigma.com/dictionary/standardized-residual/
# Because sigma depends on mass, we now need add mass in the data
# Here are the median estimates for the intercept and slope for the model sigma ~ mass
sigma_pars <- coda.samples(jm1c, n.iter = n.iter, thin = thin, variable.names = c("alpha.sigma", "b1.sigma"))

# Get coefficients from quantiles 
alpha.sigma <- data.frame(summary(sigma_pars)[2])$quantiles.50.[1]
b1.sigma <- data.frame(summary(sigma_pars)[2])$quantiles.50.[2]

# Now add in the standardized residuals
resid_df <- resid_df %>% 
  mutate(temp = data$temp, # First we need to add temp, so that we can calculate the true sigma  
         sigma = alpha.sigma + temp*b1.sigma, # calculate sigma given mass
         resid_st = resid/sigma) # calculate standardized resid

# Check linearity
p_lin <- ggplot(resid_df, aes(yhat, resid_st)) + # NOTE STANDARDIZED RESID! CHANGE IF ANOTHER MODEL
  geom_point() +
  ggtitle("Linearity")

# Check normality
p_qq <- ggplot(resid_df, aes(sample = resid_st)) + # NOTE STANDARDIZED RESID! CHANGE IF ANOTHER MODEL
  stat_qq() +
  stat_qq_line(col = "red") +
  ggtitle("QQ")

p_combo <- (p_lin + p_qq) + plot_annotation(title = 'Non-linear consumption')
p_combo & theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)

ggsave("figures/supp/non_linear/resid_non_linear_sigma.png", width = 6.5, height = 6.5, dpi = 600)


# E. PLOT PREDICTIONS ==============================================================
#** Fits and data ==================================================================
# jags.samples - Nice for summaries and predictions
# Extract the prediction at each x including credible interval
js = jags.samples(jm1c, 
                  variable.names = c("pred_small",
                                     "pred_large"), 
                  n.iter = n.iter, 
                  thin = thin)

# Save quantiles
pred_s <- summary(js$pred_small, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat
pred_l <- summary(js$pred_large, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
pred_s_df <- data.frame(lwr_95 = pred_s[1, ],
                        lwr_80 = pred_s[2, ],
                        median = pred_s[3, ],
                        upr_80 = pred_s[4, ],
                        upr_95 = pred_s[5, ],
                        temp_pred = temp_pred,
                        size = 2) 

pred_l_df <- data.frame(lwr_95 = pred_l[1, ],
                        lwr_80 = pred_l[2, ],
                        median = pred_l[3, ],
                        upr_80 = pred_l[4, ],
                        upr_95 = pred_l[5, ],
                        temp_pred = temp_pred,
                        size = 200) 

pred_df <- rbind(pred_l_df, pred_s_df)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
# Extend color palette
colourCount = length(unique(con$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p7 <- ggplot(pred_df, aes(temp_pred, median)) +
  geom_ribbon(data = filter(pred_df, size == 2), aes(x = temp_pred, ymin = lwr_95, ymax = upr_95),
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = filter(pred_df, size == 2), aes(x = temp_pred, ymin = lwr_80, ymax = upr_80),
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_ribbon(data = filter(pred_df, size == 200), aes(x = temp_pred, ymin = lwr_95, ymax = upr_95),
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = filter(pred_df, size == 200), aes(x = temp_pred, ymin = lwr_80, ymax = upr_80),
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(data = pred_df, aes(temp_pred, median, linetype = factor(size)), size = 1, alpha = 1, col = "black") +
  geom_point(data = con, aes(temp_env_ct, y_ct, fill = species),
             shape = 21, 
             size = 2, 
             alpha = 0.8, 
             color = "white") +
  #ylim(0, 3.3) +
  coord_cartesian(ylim = c(0, 3.3)) +
  scale_fill_manual(values = pal, guide = guide_legend(label.theme = element_text(angle = 0,
                                                                                  face = "italic",
                                                                                  size = 8))) +
  scale_size(range = c(2, 8), breaks = c(0, 1, 10, 100, 1000)) +
  guides(#fill = FALSE,
         size = guide_legend(override.aes = list(fill = "black",
                                                 color = "black"))) +
  labs(x = "Rescaled temperature",
       y = "ln(rescaled consumptium rate)",
       fill = "Species",
       linetype = "Size (g)") +
  NULL

pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 12),
                                       aspect.ratio = 4/5,
                                       #legend.position = "bottom", 
                                       legend.title = element_text(size = 9),
                                       legend.text = element_text(size = 8))

pWord7

ggsave("figures/supp/non_linear/non_linear_con.png", width = 6.5, height = 6.5, dpi = 600)

