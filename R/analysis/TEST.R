#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of maximum consumption rate as a function of 
# temperature with different group-effects and compare DIC and WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit models
# 
# E. Evaluate model fit
#
# D. Plot predicted mass-scaling slopes in different temperatures
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
met <- read.csv("data/met_analysis.csv")
con <- read.csv("data/con_analysis.csv")

# Filter data points at below optimum temperatures
met <- met %>% filter(above_optimum == "N")
con <- con %>% filter(above_optimum == "N")

# Data in list-format for JAGS
# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species))
con$species_n <- as.numeric(as.factor(con$species))

# Mass-range used for prediction
mass_pred_con = seq(from = min(con$log_mass_norm_ct), 
                    to = max(con$log_mass_norm_ct),
                    length.out = 100)

# Data list for consumption model
con_data = list(
  y = log(con$y), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_norm_ct,
  temp = con$temp_norm_arr_ct,
  mass_pred_con = mass_pred_con
)

# Mass-range used for prediction
mass_pred_met = seq(from = min(met$log_mass_norm_ct), 
                    to = max(met$log_mass_norm_ct),
                    length.out = 100)

# Data list for metabolism model
met_data = list(
  y = log(met$y), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_norm_ct,
  temp = met$temp_norm_arr_ct,
  mass_pred_met = mass_pred_met
)


# Test which sizes I use
p1 <- ggplot(con, aes(mass_norm, fill = species)) + 
  geom_histogram() + 
  scale_fill_viridis(discrete = TRUE) +
  coord_cartesian(expand = 0) + 
  ggtitle("Consumption") +
  labs(x = "Mass/Asymptotic mass") +
  theme_classic(base_size = 18) +
  theme(aspect.ratio = 1) +
  guides(fill = FALSE) +
  NULL

p2 <- met %>% 
  filter(mass_norm < 1) %>% 
  ggplot(., aes(mass_norm, fill = species)) + 
  geom_histogram() + 
  scale_fill_viridis(discrete = TRUE) +
  coord_cartesian(expand = 0) + 
  ggtitle("Metabolism") +
  labs(x = "Mass/Asymptotic mass") +
  theme_classic(base_size = 18) +
  theme(aspect.ratio = 1) +
  guides(fill = FALSE) +
  NULL

p1+p2


plot(met$mass_norm)
test <- filter(met, mass_norm > 1)
plot(test$mass_norm)



