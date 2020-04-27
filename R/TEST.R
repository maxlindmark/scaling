#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - TEST SCRIPT
# 
# - TEST order of factor levels...
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A. Test order of factor levels in JAGS output vs data ============================
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

# Read in your data file
con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Filter data points at below optimum temperatures
con <- con %>% filter(above_optimum == "N")

# Rename species factor for JAGS (must be numbered 1:n)
con$species_n <- as.numeric(as.factor(con$species_ab))

# Mean-center predictor variables
con$log_mass_ct <- con$log_mass - mean(con$log_mass)
con$temp_arr_ct <- con$temp_arr - mean(con$temp_arr)

# Use mass-specific values
con$y_spec <- con$y / con$mass_g

# Masses for prediction
mass_pred_con <-  seq(from = min(con$log_mass_ct), 
                      to = max(con$log_mass_ct),
                      length.out = 100)

# Temperature for prediction
temp_pred_con <- 0 # This means we use mean temperature as it is centered 

# Prepare data for JAGS
con_data = NULL

# TEST to track one species... given one species higher values and check output (random intercepts)
con$y_spec_test <- con$y_spec

con$y_spec_test <- ifelse(con$species_ab == "C.argus",
                          con$y_spec*100, 
                          con$y_spec)

con$y_spec_test <- ifelse(con$species_ab == "P.annularis",
                          con$y_spec*100, 
                          con$y_spec_test)

con$y_spec_test <- ifelse(con$species_ab == "G.affinis",
                          con$y_spec*100, 
                          con$y_spec_test)

con$y_spec_test <- ifelse(con$species_ab == "L.microlophus",
                          con$y_spec*100, 
                          con$y_spec_test)

ggplot(con, aes(y = log(con$y_spec_test), x = con$log_mass_ct, color = species_ab)) +
  geom_point()

# Data in list-format for JAGS
con_data = list(
  y = log(con$y_spec_test), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_ct,
  temp = con$temp_arr_ct,
  mass_pred = mass_pred_con,
  temp_pred = temp_pred_con
)


# C. FIT MODELS ====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

# Maximum consumption rate =========================================================
# Select model with lowest WAIC (see con_model_selection.R)
con_model = "R/analysis/JAGS_models/log_linear/selected_models/m5_pred_fit.txt"

jm_con = jags.model(con_model,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3)

update(jm_con, n.iter = n.iter) 

cs_con <- coda.samples(jm_con,
                       variable.names = c("b0"), 
                       n.iter = n.iter, 
                       thin = thin)

con_df <- data.frame(summary(cs_con)[2]) # Extract quantiles
con_df$Parameter <- rownames(con_df)

con_df %>% arrange(desc(quantiles.50.))

filter(con, species_ab %in% c("C.argus", "P.annularis", "G.affinis", "L.microlophus")) %>% 
  select(species_n, species_ab) %>% 
  distinct()

# OK, so here we see that the parameter name corresponds to "species_n", as suspected,
# because the species with highest intercepts are 11, 1, 5 and 8, which are the ones
# I increased the value for (increasing their intercept)
# Because the output is ordered by name, I need to give the species name based on their 
# order in the numerical column, not the order in which they appear in the data!



