#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of optimum growth temperatureas a function of 
# normalized body mass with different group-effects and compare DIC and WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Model selection (WAIC)
# 
# D. Model validation (convergence)
#
# E. Evaluate model fit
#
# F. Plot predictions
#
# G. Additional calculations on the posterior
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
  mass = dat$log_mass_norm_mat_ct,
  mass_pred = mass_pred
)


# C. MODEL SELECTION ===============================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample


#**** M1: Random intercept and slope ===============================================
model1 = "R/analysis/JAGS_models/T_opt/m1_T_opt_pred_fit.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 0.001,
    mu_b1 = 0.001,
    sigma = 0.001,
    sigma_b0 = 0.001,
    sigma_b1 = 0.001,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm1 = jags.model(model1,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3,
                 inits = inits)

update(jm1, n.iter = n.iter) 

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


#**** M2: Random intercept =========================================================
model2 = "R/analysis/JAGS_models/T_opt/m2_T_opt_pred_fit.txt"

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


#**** Compare both WAIC ============================================================
waic_m1
waic_m2

# > waic_m1
# [1] 180.6992
# > waic_m2
# [1] 178.417


# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm2,
                   variable.names = c("b0", "b1",
                                      "mu_b0", 
                                      "sigma_b0", "sigma"), 
                   n.iter = n.iter, 
                   thin = thin)

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
  theme(axis.text.x = element_text(size = 6)) +
  NULL
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord1 + pWord2
ggsave("figures/supp/T_opt/validation_topt_intercepts.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities of group-level means and standard deviations
p3 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "b1",
                          "sigma_b0", "sigma")) %>% 
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
  filter(Parameter %in% c("mu_b0", "b1",
                          "sigma_b0", "sigma")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord3 + pWord4
ggsave("figures/supp/T_opt/validation_topt.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p5 <- cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.002) +
  geom_point(size = 2) +
  theme(aspect.ratio = 1)+
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
ggsave("figures/supp/T_opt/validation_rhat_topt.png", width = 6.5, height = 6.5, dpi = 600)


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
tau_int <- 0.04
tau_slope <- 0.04

mu_b0 <- rnorm(25000, 0, sqrt(1/tau_int))
b1 <- rnorm(25000, 0, sqrt(1/tau_slope))

PR <- as.matrix(cbind(mu_b0, b1))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/T_opt/validation_prior_post_topt.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs,
          params = c("mu_b0", "b1"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE,
          type = "density")   # removes the trace plot

dev.off()


# E. EVALUATE MODEL FIT ============================================================
# * Note I'm not plotting how the cv's compare in data and simulations, because
# I use centered data, the cv's are extremely large...
# CODA - Nice for getting the raw posteriors
cs_fit = coda.samples(jm2,
                      variable.names = c("mean_y",
                                         "mean_y_sim", 
                                         "p_mean",
                                         "cv_y",
                                         "cv_y_sim",
                                         "p_cv"), 
                      n.iter = n.iter, 
                      thin = thin)

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df)))

# Growth
p6 <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df$p_mean), digits = 3))) +
  labs(x = "Mean simulated T_opt", y = "count") +
  coord_cartesian(expand = 0) + 
  NULL
pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)
pWord6
ggsave("figures/supp/T_opt/fit_topt_mean.png", width = 6.5, height = 6.5, dpi = 600)


# F. PLOT PREDICTIONS ==============================================================
#** Fits and data ==================================================================
# jags.samples - Nice for summaries and predictions
# Extract the prediction at each x including credible interaval
js = jags.samples(jm2, 
                  variable.names = c("pred"), 
                  n.iter = n.iter, 
                  thin = thin)

# Save quantiles
pred <- summary(js$pred, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
pred_df <- data.frame(lwr_95 = pred[1, ],
                      lwr_80 = pred[2, ],
                      median = pred[3, ],
                      upr_80 = pred[4, ],
                      upr_95 = pred[5, ],
                      mass = mass_pred + mean(dat$log_mass_norm_mat))

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
# Extend color palette
colourCount = length(unique(dat$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p7 <- ggplot(pred_df, aes(mass, median)) +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 1, alpha = 1, col = "black") +
  geom_point(data = dat, aes(log_mass_norm_mat, opt_temp_c_ct, fill = species, size = mass_g),
             shape = 21, 
             alpha = 0.8, 
             color = "white") +
  scale_fill_manual(values = pal) +
  scale_size(range = c(2, 8), breaks = c(0, 1, 10, 100, 1000)) +
  guides(fill = FALSE,
         size = guide_legend(override.aes = list(fill = "black",
                                                 color = "black"))) +
  labs(x = "ln(rescaled mass)",
       y = "Rescaled optimum growth temperature",
       size = "Mass [g]") +
  NULL

pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 12),
                                       aspect.ratio = 4/5,
                                       legend.position = "bottom", 
                                       legend.title = element_text(size = 10))
pWord7
ggsave("figures/T_opt_scatter.png", width = 6.5, height = 6.5, dpi = 600)


#** Parameter estimates ============================================================
# CODA - Nice for getting the raw posteriors
# Maximum consumption rate
cs2 <- coda.samples(jm2,
                    variable.names = c("b1"), 
                    n.iter = n.iter, 
                    thin = thin)

color_scheme_set("gray")

# Slope
p8 <- cs2 %>% 
  mcmc_dens(pars = "b1") +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Slope") +
  ggtitle("Optimum temperature ~ body mass") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  NULL
pWord8 <- p8 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))
ggsave("figures/supp/T_opt/posterior_main_param.png", width = 6.5, height = 6.5, dpi = 600)


# G. ADDITINAL CALCULATIONS ON THE POSTERIOR =======================================
# Calculate the proportion of the posterior of activation energy that is less than zero
js2 = jags.samples(jm2,
                   variable.names = c("b1"),
                   n.iter = n.iter,
                   thin = thin)

ecdf(js2$b1)(0) 
# [1] 0.9943333
