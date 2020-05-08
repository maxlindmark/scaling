#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.04.29: Max Lindmark
#
# - Code to exemplify the different temperature dependencies of metabolism and 
#   conumsumption using already fitted models. Note that the models for metabolism
#   and consumption differ a lot in what variables go in (consumption uses
#   normalized consumption ~ mass_g and centered temperature in Celcius, whereas 
#   metabolism is modelled as log(y) ~ log(centered mass) + Arrhenius temperature)
#
#   I have therefore plotted normalized rates (0:1 within each rate), to illustrate 
#   the curves
# 
# A. Load libraries
#
# B. Generate data
#
# D. Plot
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(RColorBrewer)
library(RCurl)
library(readxl)
library(magrittr)
library(viridis)
library(patchwork)
library(scales)
library(tidylog)
library(ggplot2)
library(tidyr)
library(dplyr)


# B. GENERATE DATA =================================================================
# First read in data so that we can take the means etc
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))


#*** Metabolism ====================================================================
# This is what the model is fitted to:
# data = list(
#   y = log(dat$y_spec),
#   n_obs = length(dat$y),
#   species_n = dat$species_n,
#   mass = dat$log_mass_ct,
#   temp = dat$temp_arr_ct
# )

# Find out which log_mass_ct correpond to 10 and 100 g
met$log_mass_ct <- met$log_mass - mean(met$log_mass)

log_mass_ct <- c(filter(met, mass_g == 10)$log_mass_ct[1], 
                 filter(met, mass_g == 100)$log_mass_ct[1])

# Set temperature range. Use a little bit longer range than in consumption
# In order to more easily merge data later, use the same temperature vector.
# This vector should cover the ranges we use more than enough
temp_c_ct <- seq(-20, 30, 0.1)

# Now calcualte the corresponding non-centered temperature for the metabolism case.
# This will not be used for the last plot though, there I'll only use the centered data
# because consumption data is already centered in the model fitting within species
temp_c <- temp_c_ct + mean(met$temp_c)

# Create new prediction data frame
met_pred <- data.frame(expand.grid(log_mass_ct = log_mass_ct,
                                   temp_c = temp_c))

# Add which rate it is
met_pred$rate <- "Metabolism"

# Mean-center temperature (across species) (for plotting later!)
met_pred$temp_c_ct <- met_pred$temp_c - mean(met_pred$temp_c)

# Convert non-centered temperature to Arrhenius scale for calculating rates
met_pred$temp_arr <- 1/((met_pred$temp_c + 273.15) * 8.617332e-05)

# Mean center Arrhenius temperature (which is what went in the model)
met_pred$temp_arr_ct <- met_pred$temp_arr - mean(met_pred$temp_arr)

# Add log_mass (non centered) by adding back mean
met_pred$log_mass <- met_pred$log_mass_ct + mean(met$log_mass)

# Add mass in grams
met_pred$mass_g <- exp(met_pred$log_mass)

# Calculate log metabolic rate
# Regression coefficients: see meta_cons_analysis.R
b3 <- 0.0139866
mu_b0<- -2.1669819
mu_b1 <- -0.2042507
mu_b2 <- -0.6121238

met_pred$log_y <- mu_b0 + mu_b1*met_pred$log_mass_ct + mu_b2*met_pred$temp_arr_ct + b3*met_pred$temp_arr_ct*met_pred$log_mass_ct

# Exponentiate prediction - Note it's mass-specific
met_pred$y_spec <- exp(met_pred$log_y)

# To get it into a comparable scale as consumption, standardize to max
met_pred <- met_pred %>% mutate(y_stand = y_spec/max(y_spec))

# The common temperature column will be named Temperature
met_pred$temperature <- met_pred$temp_c_ct


#*** Consumption ===================================================================
# This is what the model is fitted to:
# data = list(
#   y = con$y_ct, 
#   n_obs = length(con$y_ct), 
#   species_n = con$species_n,
#   temp = con$temp_env_ct,
#   temp_pred = temp_pred,
#   mass = con$mass_g_ct)

# Find out which mass_g_c correpond to 10 and 100 g
con$mass_g_ct <- con$mass_g - mean(con$mass_g)

mass_g_ct <- c(10 - mean(con$mass_g),
               100 - mean(con$mass_g))

# Range in Celcius
# Note that 0 in the non-linear model corresponds to the median environmental temperature,
# so it's calculated by species. I cannot recreate that here in this data frame, because I

# Use the ones in metabolism... But note it has a different interpretation here!
# That's why i keep caliing it "temp_env_ct"  in the new data frame
temp_c_ct <- seq(-20, 30, 0.1)

# Create new prediction data frame
con_pred <- data.frame(expand.grid(mass_g_ct = mass_g_ct,
                                   temp_env_ct = temp_c_ct))

# Add which rate it is
con_pred$rate <- "Maximum consumption"

# Add mass in grams
con_pred$mass_g <- con_pred$mass_g_ct + mean(con$mass_g)

# Calculate response variable
# Regression coefficients
b1 <- -0.00190553
b2 <- 0.04812750
b3 <- -0.00022713
b4 <- -0.00003802
mu_b0 <- 0.83201273

con_pred$y_ct <- mu_b0 + 
                 b1 * con_pred$mass_g_ct + 
                 b2 * con_pred$temp_env_ct + 
                 b3 * con_pred$temp_env_ct*con_pred$temp_env_ct + 
                 b4 * con_pred$temp_env_ct*con_pred$temp_env_ct*con_pred$temp_env_ct

# Divide by maximum 
con_pred <- con_pred %>% mutate(y_stand = y_ct/max(y_ct))

# The common temperature column will be named Temperature
con_pred$temperature <- con_pred$temp_env_ct

#*** Plot ==========================================================================
# Both response variables have been predicted for body masses 10 and 100 g

# The temperature variable differ however. In the consumption data, the regression
# coefficients correspond to temperature centered WITHIN species, by subtracting
# the mean temperature in the environment. 0 is this the mean habitat temperature.
# This is because the peak occurs at an optimum, and each species has been measured 
# at a different distance from that optimum. 
# For metabolism, the temperature is centered ACROSS species, such that 0 represents 
# the mean temperature in the experiments.

# The response variable in metabolism is the metabolic rate, here rescaled by 
# dividing it with the maximum. Again, an ACROSS rescaling.
# The consumption model is fitted to data dividied by mean within species (and then
# divided by mas. So it's fitted to rescaled data WITHIN species.

summary(con_pred$temperature)
summary(met_pred$temperature)

# Metabolism
ggplot(met_pred, aes(temp_c_ct, y_stand, linetype = factor(mass_g))) + 
  geom_line()

ggplot(con_pred, aes(temp_env_ct, y_stand, linetype = factor(mass_g))) + 
  geom_line()

# Combine data
met_pred_sub <- met_pred %>% select(y_stand, mass_g, temperature, rate) %>% arrange(temperature)
con_pred_sub <- con_pred %>% select(y_stand, mass_g, temperature, rate) %>% arrange(temperature)

diff <- con_pred_sub$y_stand - met_pred_sub$y_stand

dat <- data.frame(mass = rep(met_pred_sub$mass_g, 3),
                  temperature = rep(met_pred_sub$temperature, 3),
                  rate = rep(c("metabolism", "consumption", "diff"), each = nrow(met_pred_sub)),
                  y_stand = c(met_pred_sub$y_stand, con_pred_sub$y_stand, diff))

# Find peak temperature
peak <- dat %>% 
  group_by(mass) %>% 
  filter(rate == "diff") %>% 
  filter(y_stand == max(y_stand))

# Set palette
pal <- RColorBrewer::brewer.pal("Dark2", n = 3)

# Now plot all 6 curves
p1 <- ggplot(dat, aes(temperature, y_stand, color = factor(rate), linetype = factor(mass))) +
  geom_line(size = 1.2) +
  coord_cartesian(ylim = c(0, 1.1),
                  xlim = c(-5, 20), 
                  expand = 0) + 
  scale_color_manual(values = c(pal[1], pal[3], pal[2]),
                     labels = c("Consumption", "Consumption-Metabolism", "Metabolism")) +
  scale_linetype_manual(values = c(1, 2), 
                        labels = c("10g", "100g")) +
  scale_alpha_manual(values = c(0.75, 0.75, 1)) +
  labs(x = expression(paste("Rescaled temperature [", degree*C, "]")),
       y = "Rescaled rates",
       color = "Rate",
       linetype = "Mass") +
  geom_segment(data = peak, aes(x = temperature, xend = temperature, y = c(peak$y_stand[1], peak$y_stand[2]), yend = 0),
               size = 1, arrow = arrow(length = unit(0.35, "cm")), show.legend = FALSE) +
  guides(linetype = guide_legend(override.aes = list(size = 0.6))) +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       aspect.ratio = 1,
                                       legend.title = element_text(size = 10))
ggsave("figures/concept.png", width = 6.5, height = 6.5, dpi = 600)
