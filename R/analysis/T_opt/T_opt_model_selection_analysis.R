#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of optimum growth temperatures a function of 
# normalized body mass with different group-effects and compare DIC and WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Model selection (WAIC)
# 
# D. Model validation (convergence, fit, residuals)
#
# E. Plot predictions
#
# F. Additional calculations on the posterior
# 
# G. Plot temperature data
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
  mass = dat$log_mass_norm_mat_ct,
  mass_pred = mass_pred
)


# C. MODEL SELECTION ===============================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample


#**** M1: Random intercept and slope ===============================================
model1 = "JAGS_models/T_opt/m1_T_opt_pred_fit.txt"

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
model2 = "JAGS_models/T_opt/m2_T_opt_pred_fit.txt"

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

waic_m2-waic_m1

# > waic_m1
# [1] 177.3295
# > waic_m2
# [1] 178.2945
# > waic_m2-waic_m1
# [1] 0.9649845

# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm1,
                   variable.names = c("b0", "b1",
                                      "mu_b0", "mu_b1",
                                      "sigma_b0", "sigma"), 
                   n.iter = n.iter, 
                   thin = thin)

summary(cs)
# 2. Quantiles for each variable:
#               2.5%     25%       50%     75%   97.5%
# mu_b0    -0.62893 -0.2653 -0.07414  0.11706  0.48268
# mu_b1    -0.75438 -0.4495 -0.31358 -0.18053  0.15475
# sigma     1.17733  1.3901  1.51670  1.66081  1.97838
# sigma_b0  0.01205  0.1274  0.25706  0.44737  0.95766

# Evaluate convergence =============================================================
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
ggsave("figures/supp/T_opt/validation_topt_intercepts.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species slopes ================================================================
#Plot posterior densities of species intercepts
unique(cs_df$Parameter)

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
ggsave("figures/supp/T_opt/validation_topt_slopes.png", width = 6.5, height = 6.5, dpi = 600)


#**** Global means =================================================================
# Plot posterior densities of group-level means and standard deviations
p5 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1",
                          "sigma_b0", "sigma_b1", "sigma")) %>% 
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
  filter(Parameter %in% c("mu_b0", "mu_b1",
                          "sigma_b0", "sigma_b1", "sigma")) %>% 
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
ggsave("figures/supp/T_opt/validation_topt.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p7 <- cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.002) +
  geom_point(size = 2) +
  NULL
pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
ggsave("figures/supp/T_opt/validation_rhat_topt.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(0, 0.04)    
# mu_b1 ~ dnorm(0, 0.04)       
# sigma ~ dunif(0, 10) 
# sigma_b0 ~ dunif(0, 10)

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# sigma = 1
# tau <- 1/sigma^2

# Define priors for plot
tau_int <- 0.04
tau_slope <- 0.04

set.seed(42)

mu_b0 <- rnorm(25000, 0, sqrt(1/tau_int))
mu_b1 <- rnorm(25000, 0, sqrt(1/tau_slope))

PR <- as.matrix(cbind(mu_b0, mu_b1))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/T_opt/validation_prior_post_topt.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs,
          params = c("mu_b0", "mu_b1"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = TRUE,
          type = "density",
          sz_txt = 1)   # removes the trace plot

dev.off()


# Evaluate model fit & residuals =================================================
# https://rpubs.com/Niko/332320

# Extract generated data and data
cs_fit = coda.samples(jm1, n.iter = n.iter, thin = thin,
                      variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

#-- Model fit
p_fit <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df)))) +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  annotate("text", -Inf, Inf, label = round(mean(cs_fit_df$p_mean), digits = 3), 
           size = 3, hjust = -0.5, vjust = 1.3) +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL


#-- Posterior predictive distributions
# https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/#posterior-predictive-analysis

# Extract posteriors for each data point for calculation of residuals
y_sim <- coda.samples(jm1, variable.names = c("y_sim"), n.iter = n.iter, thin = thin)

# Tidy-up
df_y_sim <- ggs(y_sim)

pal <- brewer.pal(n = 3, name = "Dark2")

pp <- ggplot() +
  geom_density(data = df_y_sim, aes(value, fill = 'Posterior\nPredictive'), alpha = 0.6) +
  geom_density(data = dat, aes(opt_temp_c_ct, fill = 'Observed'), alpha = 0.6) +
  scale_fill_manual(values = pal[c(3,2)]) +
  coord_cartesian(expand = 0) +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1,
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.2, 0.95),
        legend.key.size = unit(0.3, "cm"))

#-- Residual vs fitted
df_y_sim <- df_y_sim %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = dat$opt_temp_c_ct,
         resid = y - yhat)

p_resid <- ggplot(df_y_sim, aes(yhat, resid)) +
  geom_point(fill = "black", color = "white", shape = 21) + 
  theme_classic() +
  theme(text = element_text(size = 12)) 

(p_fit | pp) / p_resid + plot_annotation(tag_levels = 'A')

ggsave("figures/supp/T_opt/fit_pp_resid_T_opt.png", width = 7, height = 7, dpi = 600)


# E. PLOT PREDICTIONS ==============================================================
#** Fits and data ==================================================================
# jags.samples - Nice for summaries and predictions
# Extract the prediction at each x including credible interval
js = jags.samples(jm1, 
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

p9 <- ggplot(pred_df, aes(mass, median)) +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.5, alpha = 1, col = "black") +
  geom_point(data = dat, aes(log_mass_norm_mat, opt_temp_c_ct, fill = species_ab, size = mass_g),
             shape = 21, alpha = 0.8, color = "grey10", stroke = 0.2) +
  scale_fill_manual(values = pal, name = "Species") +
  scale_size(range = c(1, 4), breaks = c(0, 1, 10, 100, 1000)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 1)),
         size = guide_legend(override.aes = list(fill = "black",
                                                 color = "black"))) +
  annotate("text", -8, -4, label = paste("n=", nrow(dat), sep = ""), size = 1.5,
           hjust = -0.5, vjust = 1.3) +
  labs(x = "ln(rescaled mass)",
       y = expression(paste("Rescaled ", italic(T[opt]))),
       size = "Mass [g]") +
  NULL

pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 8),
                                       aspect.ratio = 1,
                                       legend.spacing.x = unit(0, 'cm'),
                                       legend.margin = margin(-0.3, 0, 0, 0, unit = "cm"),
                                       legend.box.spacing = unit(0, 'cm'),
                                       legend.key.size = unit(0.1, 'cm'), 
                                       legend.title = element_text(size = 4),
                                       legend.text = element_text(size = 4, face = "italic"))

pWord9

ggsave("figures/T_opt_scatter.png", width = 6, height = 6, dpi = 600, unit = "cm")

# Large version of the figure
p10 <- ggplot(pred_df, aes(mass, median)) +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.5, alpha = 1, col = "black") +
  geom_point(data = dat, aes(log_mass_norm_mat, opt_temp_c_ct, fill = species_ab, size = mass_g),
             shape = 21, alpha = 0.8, color = "grey10", stroke = 0.2) +
  scale_fill_manual(values = pal, name = "Species") +
  scale_size(range = c(1.5, 6), breaks = c(0, 1, 10, 100, 1000)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 1)),
         size = guide_legend(override.aes = list(fill = "black",
                                                 color = "black"))) +
  annotate("text", -8, -4, label = paste("n=", nrow(dat), sep = ""), size = 3,
           hjust = -0.5, vjust = 1.3) +
  labs(x = "ln(rescaled mass)",
       y = expression(paste("Rescaled ", italic(T[opt]))),
       size = "Mass [g]") +
  NULL

pWord10 <- p10 + theme_classic() + theme(text = element_text(size = 12),
                                         aspect.ratio = 1,
                                         legend.spacing.x = unit(0, 'cm'),
                                         legend.margin = margin(-0.3, 0, 0, 0, unit = "cm"),
                                         legend.box.spacing = unit(0, 'cm'),
                                         legend.key.size = unit(0.1, 'cm'), 
                                         legend.title = element_text(size = 7),
                                         legend.text = element_text(size = 7, face = "italic"))

pWord10

ggsave("figures/T_opt_scatter_v2.png", width = 10, height = 10, dpi = 600, unit = "cm")



# F. ADDITINAL CALCULATIONS ON THE POSTERIOR =======================================
# Calculate the proportion of the posterior of activation energy that is less than zero
js2 = jags.samples(jm1,
                   variable.names = c("mu_b0", "mu_b1"),
                   n.iter = n.iter,
                   thin = thin)

ecdf(js2$mu_b1)(0) 
# [1] 0.926


# G. PLOT TEMPERATURE DATA =========================================================
# Plot optimum growth compared to experimental and environmental temperature
# Read in growth data as that has experimental temperature
dat2 <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

dat2$env_temp_max <- as.numeric(as.character(dat2$env_temp_max))
dat2$env_temp_min <- as.numeric(as.character(dat2$env_temp_min))

specs <- sort(unique(dat2$species_ab))
dat2$env_temp_min

# Separate environmental temperature sources (temperature in habitat or preferred)
dat2$temp_source <- 1

dat2$env_temp_min[is.na(dat2$env_temp_min)] <- -9

dat2$temp_source <- ifelse(dat2$env_temp_min == -9,
                           2,
                           dat2$temp_source)

unique(filter(dat2, env_temp_min == -9)$species_ab)

# In the data set I defined min and max of the environent temperatures as the ones
# describing the environment they live in (from Fishbase). For many species, this 
# information did not exist. For those species, I instead here use the minimum and 
# maximum temperature of the preffered temperatures, also from FishBase. These values
# where double checked 2020.04.28, prior to pre-print submission.
# Environmental temperature from Joh et al. (2013). 
dat2$env_temp_min <- ifelse(dat2$species_ab == "P.yokohamae",
                            3,
                            dat2$env_temp_min)
# Environmental temperature from Joh et al. (2013). 
dat2$env_temp_max <- ifelse(dat2$species_ab == "P.yokohamae",
                            24,
                            dat2$env_temp_max)

dat2$env_temp_min <- ifelse(dat2$species_ab == "C.lumpus",
                            0.6,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "C.lumpus",
                            11.4,
                            dat2$env_temp_max)

dat2$env_temp_min <- ifelse(dat2$species_ab == "P.olivaceus",
                            8.6,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "P.olivaceus",
                            25,
                            dat2$env_temp_max)

dat2$env_temp_min <- ifelse(dat2$species_ab == "H.hippoglossus",
                            0.4,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "H.hippoglossus",
                            7.9,
                            dat2$env_temp_max)

dat2$env_temp_min <- ifelse(dat2$species_ab == "S.maximus",
                            5.9,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "S.maximus",
                            11.9,
                            dat2$env_temp_max)

dat2$env_temp_min <- ifelse(dat2$species_ab == "A.minor",
                            0.6,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "A.minor",
                            7.6,
                            dat2$env_temp_max)


# Prepare experimental data and convert to long data frame with "source"
# indicating if it's experimental or environmental temperature
sub <- dat2 %>% 
  select(temp_c, median_temp, species_ab, log_mass, temp_source) %>% 
  gather(source, temp, 1:2) %>%
  ungroup()

filter(sub, species_ab == "A.minor")

# Raw T_opt data (using the data set we fitted models to in this script)
sub2 <- dat %>% 
  mutate(source = "Optimum (data)",
         temp = opt_temp_c,
         log_mass = log(mass_g)) %>% 
  select(temp, species_ab, source, log_mass)

head(sub2)
head(sub)

# Now combine the environmental and optimum growth data
sub <- bind_rows(sub2, sub)

# Change labels for plotting
sub$source <- ifelse(sub$source == "temp_c",
                     "Experimental\ntemperature",
                     sub$source)

sub$source <- ifelse(sub$source == "median_temp",
                     "Mid-point env. temperature",
                     sub$source)

# Because this plot will overplot, I order the sources as I want them using
# a new column with that info
sub$source2 <- 1
sub$source2 <- ifelse(sub$source == "Optimum (data)", 3, sub$source2)
sub$source2 <- ifelse(sub$source == "Mid-point env. temperature", 2, sub$source2)

str(dat2)

# Now I need to summarize the data so that one row = one species
dat3 <- dat2 %>% 
  group_by(species_ab) %>% 
  tidyr::drop_na(env_temp_max) %>% 
  tidyr::drop_na(env_temp_min) %>% 
  summarise(upper = mean(env_temp_max), 
            lower = mean(env_temp_min),
            temp = mean(temp_c),
            temp_source = mean(temp_source)) %>% 
  ungroup() %>% 
  arrange(species_ab)

dat3

# ... And do the same for T_opt data
sub3 <- dat %>%
  group_by(species_ab) %>%
  summarize(temp_opt = mean(opt_temp_c),
            source   = "Optimum (data)") %>%
  ungroup()

sub3

dat4 <- full_join(sub3, dat3) %>% drop_na() %>% arrange()

dat4

# Hard to choose which variable to sort the plot on.. doing a manual one here
unique(dat4$species_ab)
dat4$sort <- 1
dat4$sort <- ifelse(dat4$species_ab == "A.minor", 2, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "H.hippoglossus", 3, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "C.lumpus", 4, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "S.salar", 5, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "G.morhua", 6, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "S.alpinus", 7, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "S.maximus", 8, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "P.yokohamae", 9, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "P.olivaceus", 10, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "P.fulvidraco", 11, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "L.calcarifer", 12, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "R.canadum", 13, dat4$sort)

sub$sort <- 1
sub$sort <- ifelse(sub$species_ab == "A.minor", 2, sub$sort)
sub$sort <- ifelse(sub$species_ab == "H.hippoglossus", 3, sub$sort)
sub$sort <- ifelse(sub$species_ab == "C.lumpus", 4, sub$sort)
sub$sort <- ifelse(sub$species_ab == "S.salar", 5, sub$sort)
sub$sort <- ifelse(sub$species_ab == "G.morhua", 6, sub$sort)
sub$sort <- ifelse(sub$species_ab == "S.alpinus", 7, sub$sort)
sub$sort <- ifelse(sub$species_ab == "S.maximus", 8, sub$sort)
sub$sort <- ifelse(sub$species_ab == "P.yokohamae", 9, sub$sort)
sub$sort <- ifelse(sub$species_ab == "P.olivaceus", 10, sub$sort)
sub$sort <- ifelse(sub$species_ab == "P.fulvidraco", 11, sub$sort)
sub$sort <- ifelse(sub$species_ab == "L.calcarifer", 12, sub$sort)
sub$sort <- ifelse(sub$species_ab == "R.canadum", 13, sub$sort)

# Plot
pal2 <- RColorBrewer::brewer.pal("Dark2", n = 5)

p10 <- ggplot() +
  geom_point(data = sub, aes(x = reorder(species_ab, sort), y = temp,
                              fill = source, alpha = source), size = 1,
             shape = 21, color = "white", position = position_dodge(width = 1)) +
  geom_point(data = sub, aes(x = reorder(species_ab, sort), y = temp,
                             fill = source, alpha = source), size = 3,
             shape = 21, color = "white", position = position_dodge(width = 1)) +
  coord_flip() +
  geom_point(data = filter(sub, temp_source == 1 & source == "Mid-point env. temperature"),
             aes(x = reorder(species_ab, sort), y = temp),
             size = 3.5, shape = 25, fill = pal2[1], color = "white", position = position_dodge(width = 1)) +
  geom_errorbar(data = dat4, aes(x = reorder(species_ab, sort),
                                 ymin = lower, ymax = upper, linetype = factor(temp_source)),
                width = 0.5, color = pal2[1], shape = 23) +
  geom_count(data = filter(sub, source == "Optimum (data)"),
             aes(x = reorder(species_ab, sort), y = temp),
             alpha = 0.8, shape = 21, color = pal2[2], fill = pal2[2]) +
  scale_fill_manual(values = c("grey75", pal2[1], pal2[2], pal2[2])) +
  scale_alpha_manual(values = c(0.8, 0.8, 0)) +
  scale_linetype_manual(values = c(1,1)) +
  #scale_size_area() +
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,1),
                                                 color = "white")),
         linetype = FALSE,
         shape = FALSE) +
  xlab("") + 
  ylab(expression(paste("Temperature [", degree*C, "]"))) + 
  NULL 

pWord10 <- p10 + theme_classic() + theme(legend.position = "bottom",
                                         legend.direction = "vertical",
                                         axis.text.y = element_text(face = "italic"),
                                         legend.text = element_text(size = 10),
                                         aspect.ratio = 6/7,
                                         text = element_text(size = 12))

ggsave("figures/supp/T_opt/env_exp_temp.png", width = 6.5, height = 6.5, dpi = 600)

