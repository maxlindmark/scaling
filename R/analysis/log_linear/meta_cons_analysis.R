#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of metabolic and maximum consumption
# rate as a function of temperature and mass
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

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Count data
length(unique(met$common_name))
length(unique(con$common_name))
nrow(met)
nrow(con)
mean(met$temp_c)
mean(con$temp_c)
n_stand <- nrow(filter(met, type == "Standard"))
n_rout_rest <- nrow(filter(met, type %in% c("Resting", "Routine")))
n_stand / (n_stand + n_rout_rest)
1 - (n_stand / (n_stand + n_rout_rest))

summary(met$mass_g)
summary(con$mass_g)


# Filter data points at below optimum temperatures
met <- met %>% filter(above_peak_temp == "N")
con <- con %>% filter(above_peak_temp == "N")

# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species_ab))
con$species_n <- as.numeric(as.factor(con$species_ab))

# Mean-center predictor variables
met$log_mass_ct <- met$log_mass - mean(met$log_mass)
met$temp_arr_ct <- met$temp_arr - mean(met$temp_arr)

con$log_mass_ct <- con$log_mass - mean(con$log_mass)
con$temp_arr_ct <- con$temp_arr - mean(con$temp_arr)

# Use mass-specific values
met$y_spec <- met$y / met$mass_g
con$y_spec <- con$y / con$mass_g

# Masses for prediction
mass_pred_met <-  seq(from = min(met$log_mass_ct), 
                      to = max(met$log_mass_ct),
                      length.out = 100)

mass_pred_con <-  seq(from = min(con$log_mass_ct), 
                      to = max(con$log_mass_ct),
                      length.out = 100)

# Temperature for prediction
temp_pred_met <- 0 # This means we use mean temperature as it is centered 
temp_pred_con <- 0 # This means we use mean temperature as it is centered 

# Prepare data for JAGS
met_data = NULL # Clear any old data lists that might confuse things
con_data = NULL

# Data in list-format for JAGS
met_data = list(
  y = log(met$y_spec), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_ct,
  temp = met$temp_arr_ct,
  mass_pred = mass_pred_met,
  temp_pred = temp_pred_met)

# Data in list-format for JAGS
con_data = list(
  y = log(con$y_spec), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_ct,
  temp = con$temp_arr_ct,
  mass_pred = mass_pred_con,
  temp_pred = temp_pred_con)


# C. FIT MODELS ====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

# Metabolic rate ===================================================================
# Select model with lowest WAIC (see met_model_selection.R)
met_model = "R/analysis/JAGS_models/log_linear/selected_models/m2_pred_fit.txt"

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

# Maximum consumption rate =========================================================
# Select model with lowest WAIC (see con_model_selection.R)
con_model = "R/analysis/JAGS_models/log_linear/selected_models/m5_pred_fit.txt"

# Manually set initial values, because otherwise all the chains get the same
# NOTE I don't do it for all parameters...
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


# Maximum consumption rate =========================================================
# CODA - Nice for getting the raw posteriors
cs_con <- coda.samples(jm_con,
                   variable.names = c("b0", "b1", "b2", 
                                      "mu_b0", "mu_b1", "mu_b2",
                                      "sigma_b0", "sigma_b1", "sigma_b2",
                                      "sigma"), 
                   n.iter = n.iter, 
                   thin = thin)

# summary(cs_con)
# 2. Quantiles for each variable:
#   
#          2.5%     25%     50%     75%    97.5%
# mu_b0    -3.46443 -3.1237 -2.9576 -2.7944 -2.44343
# mu_b1    -0.45399 -0.3997 -0.3750 -0.3489 -0.29574
# mu_b2    -0.84680 -0.7429 -0.6945 -0.6430 -0.54058

# Convert to ggplottable data frame
cs_con_df <- ggs(cs_con)

#**** Species intercepts ===========================================================
# Plot posterior densities of species intercepts
unique(cs_con_df$Parameter)

p1 <- cs_con_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",
                          "b0[14]", "b0[15]", "b0[16]", "b0[17]", "b0[18]", "b0[19]", "b0[20]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 8),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p2 <- cs_con_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",
                          "b0[14]", "b0[15]", "b0[16]", "b0[17]", "b0[18]", "b0[19]", "b0[20]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 8),
                                       axis.text = element_text(size = 5))
pWord1 + pWord2
ggsave("figures/supp/log_linear/met_con/validation_con_intercepts.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species mass-effects =========================================================
# Plot posterior densities of species mass-effects
p3 <- cs_con_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",
                          "b1[14]", "b1[15]", "b1[16]", "b1[17]", "b1[18]", "b1[19]", "b0[20]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 8),
                                       axis.text = element_text(size = 5))


# Traceplot for evaluating chain convergence
p4 <- cs_con_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",
                          "b1[14]", "b1[15]", "b1[16]", "b1[17]", "b1[18]", "b1[19]", "b0[20]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 8),
                                       axis.text = element_text(size = 5))
pWord3 + pWord4
ggsave("figures/supp/log_linear/met_con/validation_con_mass.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species temperature-effects ==================================================
# Plot posterior densities of temperature-effects
p5 <- cs_con_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]",
                          "b2[14]", "b2[15]", "b2[16]", "b2[17]", "b2[18]", "b2[19]", "b0[20]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 8),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p6 <- cs_con_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]",
                          "b2[14]", "b2[15]", "b2[16]", "b2[17]", "b2[18]", "b2[19]", "b0[20]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 8),
                                       axis.text = element_text(size = 5))
pWord5 + pWord6
ggsave("figures/supp/log_linear/met_con/validation_con_temp.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities of group-level means and standard deviations
p7 <- cs_con_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", 
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
p8 <- cs_con_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", 
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
ggsave("figures/supp/log_linear/met_con/validation_con.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p9 <- cs_con_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.003) +
  geom_point(size = 1) +
  NULL
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
ggsave("figures/supp/log_linear/met_con/validation_rhat_con.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(0, 0.04)   # mean of all species intercept
# mu_b1 ~ dnorm(-0.25, 1)  # mean of all species mass-exponent
# mu_b2 ~ dnorm(-0.6, 1)   # mean of all species temperature coefficient

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

PR <- as.matrix(cbind(mu_b0, mu_b1, mu_b2))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/log_linear/met_con/validation_prior_post_con.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs_con,
          params = c("mu_b0", "mu_b1", "mu_b2"),
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





# E. PLOT PREDICTIONS ==============================================================
#** Fits and data ==================================================================
# jags.samples - Nice for summaries and predictions
# Extract the prediction at each x including credible interaval
# For nice labels and ln-axis:
# https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual

# Metabolism
js_met = jags.samples(jm_met, 
                      variable.names = "pred", 
                      n.iter = n.iter, 
                      thin = thin)

# Save quantiles
m_pred <- summary(js_met$pred, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
m_pred_df <- data.frame(lwr_95 = m_pred[1, ],
                        lwr_80 = m_pred[2, ],
                        median = m_pred[3, ],
                        upr_80 = m_pred[4, ],
                        upr_95 = m_pred[5, ],
                        mass = mass_pred_met)

# Consumption
js_con = jags.samples(jm_con, 
                      variable.names = "pred", 
                      n.iter = n.iter, 
                      thin = thin)

# Save quantiles
c_pred <- summary(js_con$pred, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_df <- data.frame(lwr_95 = c_pred[1, ],
                        lwr_80 = c_pred[2, ],
                        median = c_pred[3, ],
                        upr_80 = c_pred[4, ],
                        upr_95 = c_pred[5, ],
                        mass = mass_pred_con)

# This is the mean temperature for the predictions:
# round(mean(con$temp_arr))  
# round(mean(con$temp_c)) 
# round(mean(met$temp_arr))  
# round(mean(met$temp_c)) 

# Plot with mass on x and logarithmic axis instead...
# First add back the mean
m_pred_df$mass_met_non_ct <- m_pred_df$mass + mean(met$log_mass)
c_pred_df$mass_con_non_ct <- c_pred_df$mass + mean(con$log_mass)

# Then exponentiate
m_pred_df$mass_g <- exp(m_pred_df$mass_met_non_ct)
c_pred_df$mass_g <- exp(c_pred_df$mass_con_non_ct)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
# Expand color palette
colourCount = length(unique(con$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p14 <- ggplot(c_pred_df, aes(mass_g, median)) +
  geom_point(data = con, aes(mass_g, log(con$y_spec), fill = species_ab),
             size = 2, shape = 21, alpha = 0.8, color = "white") +
  geom_ribbon(data = c_pred_df, aes(x = mass_g, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = c_pred_df, aes(x = mass_g, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(trans = scales::log_trans(),
                     labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = FALSE) +
  labs(x = "mass [g]",
       y = "ln(maximum consumption rate [g/g/day])") +
  annotate("text", 0.05, 1.2, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  # annotate("text", -Inf, Inf, label = "B", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) + # This solution doesn't play with the ln axis...
  NULL

pWord14 <- p14 + theme_classic() + theme(text = element_text(size = 12))


colourCount = length(unique(met$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p15 <- ggplot(m_pred_df, aes(mass_g, median)) +
  geom_point(data = met, aes(mass_g, log(met$y_spec), fill = species_ab),
             size = 2, shape = 21, alpha = 0.8, color = "white") +
  geom_ribbon(data = m_pred_df, aes(x = mass_g, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = m_pred_df, aes(x = mass_g, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(trans = scales::log_trans(),
                     labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = FALSE) +
  labs(x = "mass [g]",
       y = "ln(metabolic rate [mg 02/g/day])") +
  annotate("text", 0.022, 0.3, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  # annotate("text", -Inf, Inf, label = "A", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) + # This solution doesn't play with the ln axis...
  NULL
pWord15 <- p15 + theme_classic() + theme(text = element_text(size = 12))

pWord14 / pWord15
ggsave("figures/pred_con_met.png", width = 6.5, height = 6.5, dpi = 600)


#** Parameter estimates ============================================================
# CODA - Nice for getting the raw posteriors
# Maximum consumption rate
cs_con <- coda.samples(jm_con,
                       variable.names = c("b0", "b1", "b2", 
                                          "mu_b0", "mu_b1", "mu_b2", 
                                          "sigma_b0", "sigma_b1", "sigma_b2",
                                          "sigma"), 
                       n.iter = n.iter, 
                       thin = thin)

# Metabolic rate - NOTE I'm also sampling b3 here! That parameter is not in consumption
cs_met <- coda.samples(jm_met,
                       variable.names = c("b0", "b1", "b2", "b3",
                                          "mu_b0", "mu_b1", "mu_b2", 
                                          "sigma_b0", "sigma_b1", "sigma_b2",
                                          "sigma"),
                       n.iter = n.iter, 
                       thin = thin)


#**** Plot species-predictions =====================================================
# Maxiumum consumption
# First get the species names in order as they appear in the output. This is based on 
# the order of the factor level (1:n), based on the species names in alphabetical order,
# which is not the same as the order they appear in the data.
con_spec <- unique(arrange(con, species_n)$species_ab)

con_df <- data.frame(summary(cs_con)[2]) # Extract quantiles
con_df$Parameter <- rownames(con_df)
con_df$Parameter_sub <- factor(substring(con_df$Parameter, 1, 2))

std_con <- data.frame(summary(cs_con)[1])
std_con$Parameter <- rownames(std_con)

#** Mass exponent
con_b <- con_df %>% filter(Parameter_sub == "b1")
con_b$Species <- con_spec
con_b$Rate <- "Maximum consumption rate"
con_b$Parameter_mte <- "Mass exponent"
con_b$pred <- filter(con_df, Parameter == "mu_b1")$quantiles.50.
con_b$pred_sd <- filter(std_con, Parameter == "mu_b1")$statistics.SD

#** Activation energy
con_e <- con_df %>% filter(Parameter_sub == "b2")
con_e$Species <- con_spec
con_e$Rate <- "Maximum consumption rate"
con_e$Parameter_mte <- "Activation energy"
con_e$pred <- filter(con_df, Parameter == "mu_b2")$quantiles.50.
con_e$pred_sd <- filter(std_con, Parameter == "mu_b2")$statistics.SD

# Metabolism
# First get the species names in order as they appear in the output. This is based on 
# the order of the factor level (1:n), based on the species names in alphabetical order,
# which is not the same as the order they appear in the data.
met_spec <- unique(arrange(met, species_n)$species_ab)

met_df <- data.frame(summary(cs_met)[2])
met_df$Parameter <- rownames(met_df)
met_df$Parameter_sub <- factor(substring(met_df$Parameter, 1, 2))

std_met <- data.frame(summary(cs_met)[1])
std_met$Parameter <- rownames(std_met)

#** Mass exponent
met_b <- met_df %>% filter(Parameter_sub == "b1")
met_b$Species <- met_spec
met_b$Rate <- "Metabolic rate"
met_b$Parameter_mte <- "Mass exponent"
met_b$pred <- filter(met_df, Parameter == "mu_b1")$quantiles.50.
met_b$pred_sd <- filter(std_met, Parameter == "mu_b1")$statistics.SD

#** Activation energy
met_e <- met_df %>% filter(Parameter_sub == "b2")
met_e$Species <- met_spec
met_e$Rate <- "Metabolic rate"
met_e$Parameter_mte <- "Activation energy"
met_e$pred <- filter(met_df, Parameter == "mu_b2")$quantiles.50.
met_e$pred_sd <- filter(std_met, Parameter == "mu_b2")$statistics.SD

#** M*T interaction
met_c <- met_df %>% filter(Parameter == "b3")
met_c$Species <- NA
met_c$Rate <- "Metabolic rate"
met_c$Parameter_mte <- "M*T interaction"
met_c$pred <- filter(met_df, Parameter == "b3")$quantiles.50.
met_c$pred_sd <- NA

# Merge data frames
df <- rbind(con_b, con_e, met_b, met_e, met_c)

# Define color palettes
pal <- brewer.pal("Dark2", n = 5)[c(1,3)]

# Convert temperature ceofficient to activation energy by multiplying with -1
df$quantiles.2.5. <- ifelse(df$Parameter_mte == "Activation energy",
                            df$quantiles.2.5. * -1,
                            df$quantiles.2.5.)

df$quantiles.25. <- ifelse(df$Parameter_mte == "Activation energy",
                           df$quantiles.25. * -1,
                           df$quantiles.25.)

df$quantiles.50. <- ifelse(df$Parameter_mte == "Activation energy",
                           df$quantiles.50. * -1,
                           df$quantiles.50.)

df$quantiles.75. <- ifelse(df$Parameter_mte == "Activation energy",
                           df$quantiles.75. * -1,
                           df$quantiles.75.)

df$quantiles.97.5. <- ifelse(df$Parameter_mte == "Activation energy",
                             df$quantiles.97.5. * -1,
                             df$quantiles.97.5.)

df$pred <- ifelse(df$Parameter_mte == "Activation energy",
                  df$pred * -1,
                  df$pred)

df$pred_sd <- ifelse(df$Parameter_mte == "Activation energy",
                     df$pred_sd * -1,
                     df$pred_sd)

# Create data frame for rectangles (prediction +/- 2*standard deviation)
df_std <- df[!duplicated(df$pred_sd), ]
df_std$ymax <- df_std$pred + 2*df_std$pred_sd
df_std$ymin <- df_std$pred - 2*df_std$pred_sd

# Plot all species varying estimates and global mean
p16 <- df %>% 
  filter(Parameter_mte %in% c("Activation energy", "Mass exponent")) %>% 
  ggplot(., aes(Species, quantiles.50., color = Rate, shape = Rate)) +
  #facet_grid(~ Parameter_mte, scales = "free") +
  facet_grid(Rate ~ Parameter_mte, scales = "free") +
  #guides(color = FALSE, fill = FALSE, shape = FALSE) +
  scale_color_manual(values = pal[1:2]) +
  scale_fill_manual(values = pal[1:2]) +
  scale_shape_manual(values = c(21, 24)) +
  geom_hline(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass exponent")), 
             aes(yintercept = pred, color = Rate),
             size = 0.6, alpha = 1, linetype = "dashed") +
  geom_rect(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass exponent")), 
            inherit.aes = FALSE, aes(ymin = ymin, ymax = ymax, fill = Rate), xmin = 0, xmax = 50, 
            alpha = 0.2) +
  coord_flip() +
  geom_errorbar(aes(Species, quantiles.50., color = Rate, 
                    ymin = quantiles.2.5., ymax = quantiles.97.5.),
                size = 1, width = 0, alpha = 0.4) +
  geom_errorbar(aes(Species, quantiles.50., color = Rate, 
                    ymin = quantiles.25., ymax = quantiles.75.), 
                size = 1.5, width = 0, alpha = 0.7) +
  geom_point(size = 1.5, fill = "white") +
  labs(x = "Species", y = "Prediction") + 
  NULL

pWord16 <- p16 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 7, face = "italic"),
                                         aspect.ratio = 2/1,
                                         legend.position = "bottom",
                                         legend.title = element_blank())
ggsave("figures/species_b_ea.png", width = 6.5, height = 6.5, dpi = 600)


#**** Plot global-predictions ======================================================
color_scheme_set("gray")

# Mass exponent
p17 <- cs_met %>% 
  mcmc_dens(pars = "mu_b1") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Mass exponent" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass exponent") +
  ggtitle("") +
  xlim(-0.6, -0.1) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Mass exponent" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL
pWord17 <- p17 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12))

# Temperature coefficient
p18 <- cs_met %>% 
  mcmc_dens(pars = "mu_b2") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Activation energy" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Temperature coefficient") +
  ggtitle("Metabolic rate") +
  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Activation energy" & Rate == "Metabolic rate")$pred*-1, # To get back to coefficient rather than activation energy because the distribution is for the coefficient 
             linetype = "dashed", color = "white") +
  NULL
pWord18 <- p18 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12))

# Mass-temperature interaction
p19 <- cs_met %>% 
  mcmc_dens(pars = "b3") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "M*T interaction" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "M*T interaction") +
  ggtitle("") +
  geom_vline(xintercept = filter(df, Parameter_mte == "M*T interaction" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "red") +
  NULL
pWord19 <- p19 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12))

# Mass exponent
p20 <- cs_con %>% 
  mcmc_dens(pars = "mu_b1") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Mass exponent" & Rate == "Maximum consumption rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass exponent") +
  ggtitle("") +
  xlim(-0.6, -0.1) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Mass exponent" & Rate == "Maximum consumption rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL
pWord20 <- p20 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12))

# Temperature coefficient
p21 <- cs_con %>% 
  mcmc_dens(pars = "mu_b2") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Activation energy" & Rate == "Maximum consumption rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Temperature coefficient") +
  ggtitle("Maximum consumption rate") +
  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Activation energy" & Rate == "Maximum consumption rate")$pred*-1, # To get back to coefficient rather than activation energy because the distribution is for the coefficient 
             linetype = "dashed", color = "white") +
  NULL
pWord21 <- p21 + theme_classic() + theme(text = element_text(size = 12),
                                         axis.text = element_text(size = 12))

pWord17 + pWord18 + pWord19 + pWord20 + pWord21 + plot_layout(ncol = 3)

ggsave("figures/supp/log_linear/met_con/posterior_main_param.png", width = 6.5, height = 6.5, dpi = 600)


# F. ADDITINAL CALCULATIONS ON THE POSTERIOR =======================================
# Calculate the proportion of the posterior of activation energy that is less than zero
js = jags.samples(jm_met,
                  variable.names = c("mu_b1", "mu_b2", "b3"),
                  n.iter = n.iter,
                  thin = thin)

ecdf(js$mu_b1)(-0.25) 
# [1] 0.8321111

1-ecdf(js$b3)(0) # how much is below?

# How much does the mass exponent decline per change in unit T?
#summary(cs)

# Coefficient is 0.016
# head(dat)
# dat$b_a <- 0.016 * dat$temp_norm_arr_ct
# 
# summary(lm(b_a ~ temp_norm_arr_ct, data = dat))

# Now fit the same exponents to C
#summary(lm(b_a ~ temp_norm, data = dat))



# Calculate the proportion of the posterior of activation energy that is less than zero
js = jags.samples(jm_con,
                  variable.names = c("mu_b1", "mu_b2"),
                  n.iter = n.iter,
                  thin = thin)

ecdf(js$mu_b1)(-0.25) 
# [1] 0.8321111

#ecdf(js$b3)(0) # how much is below?

# How much does the mass exponent decline per change in unit C?
#summary(cs)

dat <- data.frame(temp = seq(0, 20, 1))
dat$temp_arr <- 1/((dat$temp + 273.15) * 8.617332e-05)

# Coefficient is 0.014. It means per unit Arrhenius temp, the exponent declines with 0.014
dat$b_a <- 0.75 + (0.014 * dat$temp_arr)
summary(lm(b_a ~ temp_arr, data = dat))

# Now fit the same exponents to C
summary(lm(b_a ~ temp, data = dat))



