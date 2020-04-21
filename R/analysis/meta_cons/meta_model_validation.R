#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of maximum consumption rate as a function of 
# temperature with group-effects and evalute convergence
# 
# A. Load libraries
#
# B. Read data
#
# C. Model validation
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

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file(s)
dat <- read.csv("data/met_analysis.csv")

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_optimum == "N")

str(dat)
summary(dat)
glimpse(dat)

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Mass-range used for prediction
# mass_pred = seq(from = min(dat$log_mass_norm_ct), 
#                 to = max(dat$log_mass_norm_ct),
#                 length.out = 100)

mass_pred = seq(from = min(dat$log_mass_ct), 
                to = max(dat$log_mass_ct),
                length.out = 100)

# Data in list-format for JAGS
data = list(
  y = log(dat$y), 
  n_obs = length(dat$y), 
  species_n = dat$species_n,
  mass = dat$log_mass_ct,
  temp = dat$temp_arr - mean(dat$temp_arr)
)


# C. MODEL VALIDATION ==============================================================
model = "R/analysis/log-linear model/models/m2.txt"

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm, n.iter = burn.in) 


#** Sample from the posterior ======================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm,
                   variable.names = c("b0", "b1", "b2", "b3", 
                                      "mu_b0", "mu_b1", "mu_b2", 
                                      "sigma_b0", "sigma_b1", "sigma_b2",
                                      "sigma"), 
                   n.iter = samples*2, 
                   thin = 5)

summary(cs) # Get the mean estimate and SE and 95% CIs

cs_df <- data.frame(summary(cs)[1])
cs_df$Parameter <- row.names(cs_df)


#** Evaluate convergence ===========================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

test <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
# b3        0.01720 0.005241 6.766e-05      1.648e-04
# mu_b0     1.59616 0.153649 1.984e-03      1.951e-03
# mu_b1     0.76722 0.034897 4.505e-04      4.587e-04
# mu_b2    -0.61155 0.028394 3.666e-04      4.168e-04
# sigma     0.24951 0.003430 4.429e-05      4.579e-05
# sigma_b0  0.90802 0.121776 1.572e-03      2.287e-03
# sigma_b1  0.20002 0.027405 3.538e-04      4.394e-04
# sigma_b2  0.15189 0.024788 3.200e-04      4.530e-04

test2 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
# b3        0.01714 0.005228 2.134e-05      1.101e-04
# mu_b0     1.59944 0.155483 6.348e-04      7.468e-04
# mu_b1     0.76671 0.035087 1.432e-04      2.034e-04
# mu_b2    -0.61199 0.028450 1.161e-04      2.124e-04
# sigma     0.24952 0.003425 1.398e-05      1.949e-05
# sigma_b0  0.91009 0.120602 4.924e-04      1.058e-03
# sigma_b1  0.19882 0.026929 1.099e-04      2.410e-04
# sigma_b2  0.15193 0.024800 1.012e-04      2.788e-04



# Plot posterior densities of species intercepts
p1 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",  "b0[14]",
                          "b0[15]", "b0[16]", "b0[17]", "b0[18]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
p2 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]",  "b0[14]",
                          "b0[15]", "b0[16]", "b0[17]", "b0[18]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
p1+p2
#ggsave("figures/supp/model_validation_met_intercepts.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


# Plot posterior densities of species mass-effects
p3 <- cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",  "b1[14]",
                          "b1[15]", "b1[16]", "b1[17]", "b1[18]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
p4 <- cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]",  "b1[14]",
                          "b1[15]", "b1[16]", "b1[17]", "b1[18]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
p3+p4
#ggsave("figures/supp/model_validation_met_mass.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# Plot posterior densities of temperature-effects
p5 <- cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]", "b2[14]",
                          "b2[15]", "b2[16]", "b2[17]", "b2[18]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
p6 <- cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]", "b2[14]",
                          "b2[15]", "b2[16]", "b2[17]", "b2[18]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
p5+p6
#ggsave("figures/supp/model_validation_met_temp.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# Plot posterior densities of interaction-effects
# p7 <- cs_df %>% 
#   filter(Parameter %in% c("b3[1]", "b3[2]", "b3[3]", "b3[4]", "b3[5]", "b3[6]", "b3[7]", 
#                           "b3[8]", "b3[9]", "b3[10]", "b3[11]", "b3[12]", "b3[13]", "b3[14]",
#                           "b3[15]", "b3[16]", "b3[17]", "b3[18]")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   theme_classic(base_size = 11) + 
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL
# 
# # Traceplot for evaluating chain convergence
# p8 <- cs_df %>% 
#   filter(Parameter %in% c("b3[1]", "b3[2]", "b3[3]", "b3[4]", "b3[5]", "b3[6]", "b3[7]", 
#                           "b3[8]", "b3[9]", "b3[10]", "b3[11]", "b3[12]", "b3[13]", "b3[14]",
#                           "b3[15]", "b3[16]", "b3[17]", "b3[18]")) %>% 
#   ggs_traceplot(.) +
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   theme_classic(base_size = 11) + 
#   geom_line(alpha = 0.3) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   theme(axis.text.x = element_text(size = 6)) +
#   NULL
# p7+p8
#ggsave("figures/supp/model_validation_met_inter.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


# Plot posterior densities of group-level means and standard deviations
p9 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "b3",
                          "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
p10 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "b3",
                          "sigma_b0", "sigma_b1", "sigma_b2")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
p9+p10
#ggsave("figures/supp/model_validation_met.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# cs_df %>% 
#   filter(Parameter %in% c("b3")) %>% 
#   ggs_density(.) + 
#   facet_wrap(~ Parameter, ncol = 2, scales = "free") +
#   theme_classic(base_size = 11) + 
#   geom_density(alpha = 0.05) +
#   scale_color_brewer(palette = "Dark2") + 
#   scale_fill_brewer(palette = "Dark2") +
#   labs(x = "Value", y = "Density", fill = "Chain #") +
#   guides(color = FALSE, fill = FALSE) +
#   NULL


# Time series of running means
# cs_df %>% 
#   ggs_running(.) + 
#   facet_wrap(~ Parameter, ncol = 4, scales = "free") +
#   geom_line(size = 1.1, alpha = 0.8) +
#   theme_classic(base_size = 11) + 
#   scale_color_brewer(palette = rev("Dark2")) + 
#   labs(x = "Iteration", y = "Value", color = "Chain #") +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   theme(axis.text.x = element_text(size = 6)) +
#   NULL
#ggsave("figs/supplement/growth_running_mean.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)

# Rhat
cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.9998, 1.004) +
  theme_classic(base_size = 25) +
  geom_point(size = 2) +
  theme(aspect.ratio = 3/1)+
  NULL
#ggsave("figures/supp/rhat_met.pdf", plot = last_plot(), scale = 1, width = 55, height = 60, units = "cm", dpi = 300)
