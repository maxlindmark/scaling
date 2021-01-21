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
library(bayesplot)
library(tidylog)
library(tibble)

# > sessionInfo()
# other attached packages:
# [1] tibble_3.0.4       tidylog_1.0.2      scales_1.1.1       MCMCvis_0.14.0     bayesplot_1.7.2    patchwork_1.0.1    viridis_0.5.1     
# [8] viridisLite_0.3.0  magrittr_2.0.1     readxl_1.3.1       RCurl_1.98-1.2     ggmcmc_1.4.1       ggplot2_3.3.2      tidyr_1.1.2       
# [15] dplyr_1.0.2        RColorBrewer_1.1-2 rjags_4-10         coda_0.19-4


# B. READ IN DATA ==================================================================
# Read in your data file
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

length(unique(con$common_name))

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

# Add dummy variable for metabolic type (0 for routine/resting, 1 for standard)
met$met_type <- ifelse(met$type == "Standard", 1, 0)
met$met_s <- ifelse(met$type == "Standard", 1, 0)
met$met_r <- ifelse(!met$type == "Standard", 1, 0)

plot(met$met_s ~ met$met_r)

# Prepare data for JAGS
met_data = NULL # Clear any old data lists that might confuse things
con_data = NULL

# Arrange by species name so there's no confusion between the index and the name
met <- met %>% arrange(species_n)
con <- con %>% arrange(species_n)

# Data in list-format for JAGS
met_data = list(
  y = log(met$y_spec), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_ct,
  temp = met$temp_arr_ct,
  met_r = met$met_r,
  met_s = met$met_s,
  mass_pred = mass_pred_met,
  temp_pred = temp_pred_met,
  spec_s = distinct(met, species, .keep_all = TRUE)$met_s, # New variable to go in the 1:34 for loop 
  spec_r = distinct(met, species, .keep_all = TRUE)$met_r) # New variable to go in the 1:34 for loop

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
met_model = "JAGS_models/log_linear/selected_models/m1_pred_fit_met.txt"

# Manually set initial values, because otherwise all the chains get the same
inits_met = list(
  list(
    mu_b0_s = 0.1,
    mu_b0_r = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    mu_b3 = 1,
    sigma = 0.1,
    sigma_b0_r = 0.1,
    sigma_b0_s = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 1,
    mu_b0_r = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    mu_b3 = 1,
    sigma = 1,
    sigma_b0_r = 1,
    sigma_b0_s = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0_s = 2,
    mu_b0_r = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    mu_b3 = 2,
    sigma = 2,
    sigma_b0_r = 2,
    sigma_b0_s = 2,
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
con_model = "JAGS_models/log_linear/selected_models/m5_pred_fit_con.txt"

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

update(jm_con, n.iter = n.iter, inits = inits_con) 


# D. MODEL VALIDATION ==============================================================
# Metabolic rate ===================================================================
# CODA - Nice for getting the raw posteriors
cs_met <- coda.samples(jm_met,
                       variable.names = c("b0_r", "b0_s", "b1", "b2", "b3",
                                          "mu_b0_r", "mu_b0_s", "mu_b1", "mu_b2", "mu_b3",
                                          "sigma_b0_r", "sigma_b0_s", "sigma_b1", "sigma_b2", "sigma_b3",
                                          "sigma"), 
                       n.iter = n.iter, 
                       thin = thin)

# summary(cs_met)
# 2. Quantiles for each variable:
#              2.5%       25%        50%      75%     97.5%
# ...
# mu_b0_r    -2.1087401 -1.990029 -1.9310915 -1.869688 -1.741797
# mu_b0_s    -2.8085439 -2.591969 -2.4874513 -2.391457 -2.169679
# mu_b1      -0.2584772 -0.225095 -0.2088925 -0.193078 -0.160210
# mu_b2      -0.6731687 -0.638018 -0.6209557 -0.604735 -0.572977
# mu_b3       0.0008228  0.011839  0.0174556  0.023837  0.036746
# ...

# Convert to ggplottable data frame
cs_met_df <- ggs(cs_met)


#**** Species metabolic type effects standard (1/2 because to many species) ========
# Plot posterior densities of species intercepts
unique(cs_met_df$Parameter)

p1a <- cs_met_df %>% 
  filter(Parameter %in% c("b0_s[1]", "b0_s[2]", "b0_s[3]", "b0_s[4]", "b0_s[5]", "b0_s[6]", "b0_s[7]", 
                          "b0_s[8]", "b0_s[9]", "b0_s[10]", "b0_s[11]", "b0_s[12]", "b0_s[13]",
                          "b0_s[14]", "b0_s[15]", "b0_s[16]", "b0_s[17]")) %>% 
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
  filter(Parameter %in% c("b0_s[1]", "b0_s[2]", "b0_s[3]", "b0_s[4]", "b0_s[5]", "b0_s[6]", "b0_s[7]", 
                          "b0_s[8]", "b0_s[9]", "b0_s[10]", "b0_s[11]", "b0_s[12]", "b0_s[13]",
                          "b0_s[14]", "b0_s[15]", "b0_s[16]", "b0_s[17]")) %>% 
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
ggsave("figures/supp/log_linear/met_con/validation_met_type_effects_s_1.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species metabolic type effects standard (2/2 because to many species) ========
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1b <- cs_met_df %>% 
  filter(Parameter %in% c("b0_s[18]", "b0_s[19]",
                          "b0_s[20]", "b0_s[21]", "b0_s[22]", "b0_s[23]", "b0_s[24]", "b0_s[25]",
                          "b0_s[25]", "b0_s[26]", "b0_s[27]", "b0_s[28]", "b0_s[29]", "b0_s[30]",
                          "b0_s[31]", "b0_s[32]", "b0_s[33]", "b0_s[34]")) %>% 
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
  filter(Parameter %in% c("b0_s[18]", "b0_s[19]",
                          "b0_s[20]", "b0_s[21]", "b0_s[22]", "b0_s[23]", "b0_s[24]", "b0_s[25]",
                          "b0_s[25]", "b0_s[26]", "b0_s[27]", "b0_s[28]", "b0_s[29]", "b0_s[30]",
                          "b0_s[31]", "b0_s[32]", "b0_s[33]", "b0_s[34]")) %>% 
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
ggsave("figures/supp/log_linear/met_con/validation_met_type_effects_s_2.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species metabolic type effects routine (1/2 because to many species) =================
# Plot posterior densities of species intercepts
unique(cs_met_df$Parameter)

p1a <- cs_met_df %>% 
  filter(Parameter %in% c("b0_r[1]", "b0_r[2]", "b0_r[3]", "b0_r[4]", "b0_r[5]", "b0_r[6]", "b0_r[7]", 
                          "b0_r[8]", "b0_r[9]", "b0_r[10]", "b0_r[11]", "b0_r[12]", "b0_r[13]",
                          "b0_r[14]", "b0_r[15]", "b0_r[16]", "b0_r[17]")) %>% 
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
  filter(Parameter %in% c("b0_r[1]", "b0_r[2]", "b0_r[3]", "b0_r[4]", "b0_r[5]", "b0_r[6]", "b0_r[7]", 
                          "b0_r[8]", "b0_r[9]", "b0_r[10]", "b0_r[11]", "b0_r[12]", "b0_r[13]",
                          "b0_r[14]", "b0_r[15]", "b0_r[16]", "b0_r[17]")) %>% 
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
ggsave("figures/supp/log_linear/met_con/validation_met_type_effects_r_1.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species metabolic type effects routine (2/2 because to many species) =========
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1b <- cs_met_df %>% 
  filter(Parameter %in% c("b0_r[18]", "b0_r[19]",
                          "b0_r[20]", "b0_r[21]", "b0_r[22]", "b0_r[23]", "b0_r[24]", "b0_r[25]",
                          "b0_r[25]", "b0_r[26]", "b0_r[27]", "b0_r[28]", "b0_r[29]", "b0_r[30]",
                          "b0_r[31]", "b0_r[32]", "b0_r[33]", "b0_r[34]")) %>% 
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
  filter(Parameter %in% c("b0_r[18]", "b0_r[19]",
                          "b0_r[20]", "b0_r[21]", "b0_r[22]", "b0_r[23]", "b0_r[24]", "b0_r[25]",
                          "b0_r[25]", "b0_r[26]", "b0_r[27]", "b0_r[28]", "b0_r[29]", "b0_r[30]",
                          "b0_r[31]", "b0_r[32]", "b0_r[33]", "b0_r[34]")) %>% 
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
ggsave("figures/supp/log_linear/met_con/validation_met_type_effects_r_2.2.png", width = 6.5, height = 6.5, dpi = 600)


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


#**** Species interaction-effects (1/2 because to many species) ====================
# Plot posterior densities of temperature-effects
p7a <- cs_met_df %>% 
  filter(Parameter %in% c("b3[1]", "b3[2]", "b3[3]", "b3[4]", "b3[5]", "b3[6]", "b3[7]", 
                          "b3[8]", "b3[9]", "b3[10]", "b3[11]", "b3[12]", "b3[13]",
                          "b3[14]", "b3[15]", "b3[16]", "b3[17]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord7a <- p7a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p8a <- cs_met_df %>% 
  filter(Parameter %in% c("b3[1]", "b3[2]", "b3[3]", "b3[4]", "b3[5]", "b3[6]", "b3[7]", 
                          "b3[8]", "b3[9]", "b3[10]", "b3[11]", "b3[12]", "b3[13]",
                          "b3[14]", "b3[15]", "b3[16]", "b3[17]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord8a <- p8a + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord7a + pWord8a
ggsave("figures/supp/log_linear/met_con/validation_met_inter_1.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species interaction-effects (2/2 because to many species) ====================
# Plot posterior densities of temperature-effects
p7b <- cs_met_df %>% 
  filter(Parameter %in% c("b3[18]", "b3[19]",
                          "b3[20]", "b3[21]", "b3[22]", "b3[23]", "b3[24]", "b3[25]",
                          "b3[25]", "b3[26]", "b3[27]", "b3[28]", "b3[29]", "b3[30]",
                          "b3[31]", "b3[32]", "b3[33]", "b3[34]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord7b <- p7b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p8b <- cs_met_df %>% 
  filter(Parameter %in% c("b3[18]", "b3[19]",
                          "b3[20]", "b3[21]", "b3[22]", "b3[23]", "b3[24]", "b3[25]",
                          "b3[25]", "b3[26]", "b3[27]", "b3[28]", "b3[29]", "b3[30]",
                          "b3[31]", "b3[32]", "b3[33]", "b3[34]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord8b <- p8b + theme_classic() + theme(text = element_text(size = 8),
                                         axis.text = element_text(size = 5))
pWord7b + pWord8b
ggsave("figures/supp/log_linear/met_con/validation_met_inter_2.2.png", width = 6.5, height = 6.5, dpi = 600)


#**** Global means =================================================================
# Plot posterior densities of global means and standard deviations
p9 <- cs_met_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b0_r", "mu_b0_s", "mu_b1", "mu_b2", "mu_b3",
                          "sigma_b0_r", "sigma_b0_s", "sigma_b1", "sigma_b2", "sigma_b3")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p10 <- cs_met_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b0_r", "mu_b0_s", "mu_b1", "mu_b2", "mu_b3",
                          "sigma_b0_r", "sigma_b0_s", "sigma_b1", "sigma_b2", "sigma_b3")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord10 <- p10 + theme_classic() + theme(text = element_text(size = 10),
                                         axis.text = element_text(size = 5))
pWord9 + pWord10
ggsave("figures/supp/log_linear/met_con/validation_met.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergence (Rhat) ====================================================
p11 <- cs_met_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.02) +
  geom_point(size = 1) +
  NULL
pWord11 <- p11 + theme_classic() + theme(text = element_text(size = 10),
                                         axis.text = element_text(size = 5))
ggsave("figures/supp/log_linear/met_con/validation_rhat_met.png", width = 6.5, height = 12, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS (global means)
# mu_b0_s ~ dnorm(-2, 0.04)
# mu_b0_r ~ dnorm(-1, 0.04)
# mu_b1 ~ dnorm(-0.25, 1)  
# mu_b2 ~ dnorm(-0.6, 1)   
# mu_b3 ~ dnorm(0, 1)     

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# tau <- 1/sigma^2
# from tau to sigma 
# sigma <- sqrt(1/tau)

# Define priors for plot
tau <- 1
tau_int <- 0.04

mu_b0_s <- rnorm(25000, -2, sqrt(1/tau_int))
mu_b0_r <- rnorm(25000, -1, sqrt(1/tau_int)) 
mu_b1 <- rnorm(25000, -0.25, sqrt(1/tau)) # hist(rnorm(25000, -2, sqrt(1/tau)))
mu_b2 <- rnorm(25000, -0.6, sqrt(1/tau))   
mu_b3 <- rnorm(25000, 0, sqrt(1/tau))         

PR <- as.matrix(cbind(mu_b0_s, mu_b0_r, mu_b1, mu_b2, mu_b3))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/log_linear/met_con/validation_prior_post_met.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs_met,
          params = c("mu_b0_s", "mu_b0_r", "mu_b1", "mu_b2", "mu_b3"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = TRUE,
          type = "density")   # removes the trace plot

dev.off()

# hist(mu_b3, xlim = c(-0.3, -0.1)) # It's flat...


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
#           2.5%     25%     50%     75%    97.5%
# ...
# mu_b0    -3.46443 -3.1237 -2.9576 -2.7944 -2.44343
# mu_b1    -0.45399 -0.3997 -0.3750 -0.3489 -0.29574
# mu_b2    -0.84680 -0.7429 -0.6945 -0.6430 -0.54058
# ...

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


#**** Global means =================================================================
# Plot posterior densities of global means and standard deviations
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
  xlim(0.999, 1.005) +
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
          Rhat = FALSE,
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

#-- Model fit
p_fit_m <- ggplot(cs_fit_df_met, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df_met)))) +
  geom_vline(xintercept = cs_fit_df_met$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  annotate("text", -Inf, Inf, label = round(mean(cs_fit_df_met$p_mean), digits = 3), 
           size = 3, hjust = -0.5, vjust = 1.3) +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL


#-- Posterior predictive distributions
# https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/#posterior-predictive-analysis

# Extract posteriors for each data point for calculation of residuals
y_sim_met <- coda.samples(jm_met, variable.names = c("y_sim"), n.iter = n.iter, thin = thin)

# Tidy-up
df_y_sim_met <- ggs(y_sim_met)

pal <- brewer.pal(n = 3, name = "Dark2")

pp_m <- ggplot() +
  geom_density(data = df_y_sim_met, aes(value, fill = 'Posterior\nPredictive'), alpha = 0.6) +
  geom_density(data = met, aes(log(y_spec), fill = 'Observed'), alpha = 0.6) +
  scale_fill_manual(values = pal[c(3,2)]) +
  coord_cartesian(expand = 0) +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1,
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.2, 0.95),
        legend.key.size = unit(0.3, "cm"))

#-- Residual vs fitted
df_y_sim_met <- df_y_sim_met %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = met_data$y,
         resid = y - yhat)

p_resid_met <- ggplot(df_y_sim_met, aes(yhat, resid)) +
  geom_point(fill = "black", color = "white", shape = 21) + 
  theme_classic() +
  theme(text = element_text(size = 12)) 

(p_fit_m | pp_m) / p_resid_met + plot_annotation(tag_levels = 'A')

ggsave("figures/supp/log_linear/met_con/fit_pp_resid_met.png", width = 7, height = 7, dpi = 600)


#**** Consumption ===================================================================
# Extract generated data and data
cs_fit_con = coda.samples(jm_con, n.iter = n.iter, thin = thin,
                          variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))

#-- Model fit
p_fit_c <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df_con)))) +
  geom_vline(xintercept = cs_fit_df_con$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  annotate("text", -Inf, Inf, label = round(mean(cs_fit_df_con$p_mean), digits = 3), 
           size = 3, hjust = -0.5, vjust = 1.3) +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL


#-- Posterior predictive distributions
# https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/#posterior-predictive-analysis

# Extract posteriors for each data point for calculation of residuals
y_sim_con <- coda.samples(jm_con, variable.names = c("y_sim"), n.iter = n.iter, thin = thin)

# Tidy-up
df_y_sim_con <- ggs(y_sim_con)

pal <- brewer.pal(n = 3, name = "Dark2")

pp_c <- ggplot() +
  geom_density(data = df_y_sim_con, aes(value, fill = 'Posterior\nPredictive'), alpha = 0.6) +
  geom_density(data = con, aes(log(y_spec), fill = 'Observed'), alpha = 0.6) +
  scale_fill_manual(values = pal[c(3,2)]) +
  coord_cartesian(expand = 0) +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1,
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.2, 0.95),
        legend.key.size = unit(0.3, "cm"))


#-- Residual vs fitted
df_y_sim_con <- df_y_sim_con %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = con_data$y,
         resid = y - yhat)

p_resid_con <- ggplot(df_y_sim_con, aes(yhat, resid)) +
  geom_point(fill = "black", color = "white", shape = 21) + 
  theme_classic() +
  theme(text = element_text(size = 12)) 

(p_fit_c | pp_c) / p_resid_con + plot_annotation(tag_levels = 'A')

ggsave("figures/supp/log_linear/met_con/fit_pp_resid_con.png", width = 7, height = 7, dpi = 600)


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

# Next sample parameter estimates, as the total plot is combination 

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
                       variable.names = c("b0_s", "b0_r", "b1", "b2", "b3",
                                          "mu_b0_s", "mu_b0_r", "mu_b1", "mu_b2", "mu_b3",
                                          "sigma_b0_s", "sigma_b0_r", "sigma_b1", "sigma_b2", "sigma_b3",
                                          "sigma"),
                       n.iter = n.iter, 
                       thin = thin)


#**** Plot species-intercepts ======================================================
cs_met2 <- coda.samples(jm_met,
                        variable.names = c("b0_s", "mu_b0_s", "b0_r", "mu_b0_r"), 
                        n.iter = n.iter, 
                        thin = thin)

summary(cs_met2)[2]

species_intercepts <- data.frame(summary(cs_met2)[2]) %>% 
  rownames_to_column("param") %>% 
  mutate(param = factor(param)) %>% 
  mutate(standard = stringr::str_detect(param, 'b0_s')) %>% 
  mutate(met_type_est = ifelse(standard == TRUE, "standard", "routine/resting")) %>% 
  rownames_to_column("species_n") %>%  
  mutate(species_n = as.numeric(species_n)) %>% 
  mutate(species_n = factor(ifelse(species_n > 34, species_n-34, species_n)))

# Now find which species(index) have which data
spec_type <- met %>%
  group_by(species_n) %>%
  distinct(met_type) %>%
  ungroup() %>% 
  mutate(species_n = factor(species_n)) %>%
  mutate(met_type_data = factor(ifelse(met_type == 1, "standard", "routine/resting"))) %>%
  select(species_n, met_type_data) %>% 
  data.frame()

species_intercepts <- left_join(species_intercepts, spec_type)

pal = RColorBrewer::brewer.pal(n = 3, name = "Dark2")

# Add in the correct species names
met_names <- met %>%
  distinct(common_name, .keep_all = TRUE) %>% 
  dplyr::select(species_ab, species_n) %>% 
  mutate(species_n = factor(species_n))

species_intercepts <- left_join(species_intercepts, met_names, by = "species_n")

species_intercepts %>% 
  drop_na(met_type_data) %>% 
  filter(!quantiles.50. == 0) %>% 
  ggplot(., aes(factor(species_ab), quantiles.50., color = factor(met_type_data))) + 
  geom_errorbar(aes(ymin = quantiles.2.5., ymax = quantiles.97.5., x = factor(species_ab))) +
  geom_hline(yintercept = filter(species_intercepts, param %in% c("mu_b0_r"))$quantiles.50.,
             alpha = 0.5, linetype = 1, color = pal[1], size = 1) +
  geom_hline(yintercept = filter(species_intercepts, param %in% c("mu_b0_s"))$quantiles.50.,
             alpha = 0.5, linetype = 1, color = pal[2], size = 1) +
  geom_hline(yintercept = filter(species_intercepts, param %in% c("mu_b0_r"))$quantiles.97.5.,
             alpha = 0.5, linetype = 2, color = pal[1], size = 1) +
  geom_hline(yintercept = filter(species_intercepts, param %in% c("mu_b0_s"))$quantiles.97.5.,
             alpha = 0.5, linetype = 2, color = pal[2], size = 1) +
  geom_hline(yintercept = filter(species_intercepts, param %in% c("mu_b0_r"))$quantiles.2.5.,
             alpha = 0.5, linetype = 2, color = pal[1], size = 1) +
  geom_hline(yintercept = filter(species_intercepts, param %in% c("mu_b0_s"))$quantiles.2.5.,
             alpha = 0.5, linetype = 2, color = pal[2], size = 1) +
  scale_color_brewer(palette = "Dark2", name = "Metabolic\ntype") +
  geom_point(size = 3) +
  coord_flip() +
  #facet_wrap(~ met_type_est) +
  labs(y = "Intercept", x = "") + 
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("figures/supp/log_linear/met_con/metabolism_intercepts_type.png", width = 7, height = 7, dpi = 600)

# Test it's correct using a few species...
# head(filter(met, species_n == 4), 1) %>% dplyr::select(species_n, type, met_type, species_ab)
# head(filter(met, species_n == 17), 1) %>% dplyr::select(species_n, type, met_type, species_ab)
# head(filter(met, species_n == 23), 1) %>% dplyr::select(species_n, type, met_type, species_ab)
# head(filter(met, species_n == 31), 1) %>% dplyr::select(species_n, type, met_type, species_ab)
# head(filter(met, species_n == 9), 1) %>% dplyr::select(species_n, type, met_type, species_ab)
# head(filter(met, species_n == 22), 1) %>% dplyr::select(species_n, type, met_type, species_ab)


#**** Plot species-predictions =====================================================
# Maximum consumption
# First get the species names and match them to the parameter (1:n species)
con_spec <- unique(con$species_ab)

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
# First get the species names and match them to the parameter (1:n species)
met_spec <- unique(met$species_ab)

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
met_c <- met_df %>% filter(Parameter_sub == "b3")
met_c$Species <- met_spec
met_c$Rate <- "Metabolic rate"
met_c$Parameter_mte <- "M*T interaction"
met_c$pred <- filter(met_df, Parameter == "mu_b3")$quantiles.50.
met_c$pred_sd <- filter(std_met, Parameter == "mu_b3")$statistics.SD

# Merge data frames
df <- rbind(con_b, con_e, met_b, met_e, met_c)

# Define color palettes
pal <- brewer.pal("Dark2", n = 5)[c(1,3)]

# Convert temperature coefficient to activation energy by multiplying with -1
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


#** Make a combined plot ===========================================================
# Rates vs data and species-level intercepts (3 in one)
colourCount = length(unique(con$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p14 <- ggplot(c_pred_df, aes(mass_g, median)) +
  geom_point(data = con, aes(mass_g, log(y_spec), fill = species_ab),
             size = 2, shape = 21, alpha = 0.8, color = "white") +
  geom_ribbon(data = c_pred_df, aes(x = mass_g, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = c_pred_df, aes(x = mass_g, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.8, alpha = 0.8) +
  coord_cartesian(ylim = c(min(min(con_data$y), min(met_data$y)), max(max(con_data$y), max(met_data$y)))) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(trans = scales::log_trans(),
                     #labels = scales::number_format(accuracy = .1), # Use this to get evenly space ticks, and the below to round them up!
                     breaks = c(0.5, 7, 150)) +
  guides(fill = FALSE) +
  labs(x = "mass [g]",
       y = "ln(maximum consumption rate [g/g/day])") +
  annotate("text", 0.1, -8, label = paste("n=", nrow(con), sep = ""), size = 3,
           hjust = -0.5, vjust = 1.3) +
  ggtitle("A") +
  NULL

pWord14 <- p14 + theme_classic() + theme(text = element_text(size = 14),
                                         aspect.ratio = 1,
                                         plot.title = element_text(face = "bold"))

colourCount = length(unique(met$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p15 <- ggplot(m_pred_df, aes(mass_g, median)) +
  geom_point(data = met, aes(mass_g, log(y_spec), fill = species_ab),
             size = 2, shape = 21, alpha = 0.8, color = "white") +
  geom_ribbon(data = m_pred_df, aes(x = mass_g, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = m_pred_df, aes(x = mass_g, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.8, alpha = 0.8) +
  coord_cartesian(ylim = c(min(min(con_data$y), min(met_data$y)), max(max(con_data$y), max(met_data$y)))) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(trans = scales::log_trans(),
                     #labels = scales::number_format(accuracy = .1), # Use this to get evenly space ticks, and the below to round them up!
                     breaks = c(0.5, 20, 1100)) +
  guides(fill = FALSE) +
  labs(x = "mass [g]",
       y = expression(paste("ln(metabolic rate [mg ", O[2], "/g/day]"))) +
  # annotate("text", 0.05, 1.2, label = "B", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", 0.1, -8, label = paste("n=", nrow(met), sep = ""), size = 3,
           hjust = -0.5, vjust = 1.3) +
  ggtitle("B") +
  NULL

pWord15 <- p15 + theme_classic() + theme(text = element_text(size = 14),
                                         aspect.ratio = 1,
                                         plot.title = element_text(face = "bold"))

pWord14 | pWord15

# Add a fake datapoint to get identical axis ranges for both rates
ylim_activ <- c(max(filter(df, Parameter_mte == "Activation energy")$quantiles.2.5.),
                min(filter(df, Parameter_mte == "Activation energy")$quantiles.97.5.))

ylim_massexp <- c(min(filter(df, Parameter_mte == "Mass exponent")$quantiles.2.5.),
                  max(filter(df, Parameter_mte == "Mass exponent")$quantiles.97.5.))

plot_dat <- data.frame(quantiles.50. = c(ylim_activ, ylim_massexp),
                       Parameter_mte = rep(c("Activation energy", "Mass exponent"), each = 2),
                       Species = rep("O.mykiss", 4))

# Define color palettes
pal <- brewer.pal("Dark2", n = 5)[c(1,3)]

# Plot all species varying estimates and global mean
p16 <- df %>% 
  filter(Parameter_mte %in% c("Activation energy", "Mass exponent")) %>% 
  ggplot(., aes(Species, quantiles.50., color = Rate, shape = Rate)) +
  facet_grid(Parameter_mte ~ Rate, scales = "free") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_shape_manual(values = c(21, 24)) +
  geom_point(data = plot_dat, aes(Species, quantiles.50.), color = "white", inherit.aes = FALSE) +
  geom_hline(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass exponent")),
             aes(yintercept = pred, color = Rate),
             size = 0.6, alpha = 1, linetype = "dashed") +
  geom_rect(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass exponent")),
            inherit.aes = FALSE, aes(ymin = ymin, ymax = ymax, fill = Rate), xmin = 0, xmax = 50,
            alpha = 0.2) +
  geom_errorbar(aes(Species, quantiles.50., color = Rate, 
                    ymin = quantiles.2.5., ymax = quantiles.97.5.),
                size = 1, width = 0, alpha = 0.4) +
  geom_errorbar(aes(Species, quantiles.50., color = Rate, 
                    ymin = quantiles.25., ymax = quantiles.75.), 
                size = 1.5, width = 0, alpha = 0.7) +
  geom_point(size = 1.5, fill = "white") +
  labs(x = "", y = "") + 
  guides(color = FALSE, shape = FALSE, fill = FALSE) +
  ggtitle("C") +
  NULL 

pWord16 <- p16 + theme_classic() + theme(text = element_text(size = 14),
                                         axis.text.x = element_text(size = 7, face = "italic", angle = 90),
                                         title = element_text(face = "bold"),
                                         aspect.ratio = 1/2,
                                         legend.position = "bottom",
                                         legend.title = element_blank())

(pWord14 + pWord15) / pWord16 + 
  plot_layout(heights = unit(c(7.15, 1), c('cm', 'null'))) #+
  #plot_annotation(tag_levels = "A")

ggsave("figures/meta_cons_combined.png", width = 22, height = 22, dpi = 600, units = "cm")
ggsave("figures/meta_cons_combined.pdf", width = 22, height = 22, dpi = 600, units = "cm")


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
  mcmc_dens(pars = "mu_b3") +
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
                  variable.names = c("mu_b1", "mu_b2", "mu_b3"),
                  n.iter = n.iter,
                  thin = thin)

1-ecdf(js$mu_b1)(-0.25) 
# [1] 0.9505556

1-ecdf(js$mu_b3)(0) # how much is above 0??

# How much does the mass exponent decline per change in unit C?
#summary(cs)
dat <- data.frame(temp = seq(0, 20, 1))
dat$temp_arr <- 1/((dat$temp + 273.15) * 8.617332e-05)

# Coefficient is 0.018. It means per unit Arrhenius temp, the exponent declines with 0.018
dat$b_a <- 0.75 + (0.018 * dat$temp_arr)
summary(lm(b_a ~ temp_arr, data = dat))

# Now fit the same exponents to C
summary(lm(b_a ~ temp, data = dat))


# Calculate the proportion of the posterior of activation energy that is less than zero
js = jags.samples(jm_con,
                  variable.names = c("mu_b1", "mu_b2"),
                  n.iter = n.iter,
                  thin = thin)

ecdf(js$mu_b1)(-0.25) 
# [1] 0.9964444


