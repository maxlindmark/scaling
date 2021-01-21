#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.12.21: Max Lindmark
#
# Short script to see what the mass-exponent of net gain is given metabolism and
# consumption, and how that relates to the exponent for growth rate
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load libraries, install if needed
library(dplyr)
library(RCurl)
library(ggplot2)

# First read in data so that we can take the means etc
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))


#**** Metabolism ===================================================================
# Set mass range
mass_range <- 1:500

# Create new prediction data frame
met_pred <- data.frame(mass_g = mass_range)

# Add which rate it is
met_pred$rate <- "Metabolism"

# Calculate log mass
met_pred$log_mass <- log(met_pred$mass_g)

# Calculate mean-centered log mass (as it's used for model fitting)
met_pred$log_mass_ct <- met_pred$log_mass - mean(met$log_mass)

# Calculate log metabolic rate
# Regression coefficients: see meta_cons_analysis.R
# Using the intercept for routine metabolic rate
# mu_b0 <- -1.9310915
# mu_b1 <- -0.2088925

# Temperature-independent allometric function (valid at 19C, which is the mean temperature)
met_pred$log_y_spec <- mu_b0 + mu_b1*met_pred$log_mass_ct

# Exponentiate prediction - Note it's mass-specific
met_pred$y_spec <- exp(met_pred$log_y_spec)

# Whole organism rate
met_pred$y <- met_pred$y_spec * met_pred$mass_g

# Now convert to g/d using the same values as Jan
# 1 kcal = 295 mg O2 
# 0.003389831 kcal = 1 mg 02
met_pred$y_kcal_h <- met_pred$y * 0.003389831

# Convert to cal
met_pred$y_cal_h <- met_pred$y_kcal_h*1000

# Convert to Joule
met_pred$y_J_h <- met_pred$y_cal_h * 4.1855

# Convert to g/h. Energy content of fish: 5600 J/g
met_pred$y_g_h <- met_pred$y_J_h / 5600

# Convert to g/d
met_pred$y_g_d <- met_pred$y_g_h * 24


#**** Consumption ==================================================================
# Create new prediction data frame
con_pred <- data.frame(mass_g = mass_range)

# Add which rate it is
con_pred$rate <- "Maximum consumption"

# Calculate log mass
con_pred$log_mass <- log(con_pred$mass_g)

# Calculate mean-centered mass (as its used for model fitting)
con_pred$log_mass_ct <- con_pred$log_mass - mean(con$log_mass)

# Calculate log maximum consumption rate
# Regression coefficients: see meta_cons_analysis.R
mu_b0 <- -2.9576
mu_b1 <- -0.3757

# Temperature-independent allometric function (valid at 19C, which is the mean temperature)
con_pred$log_y_spec <- mu_b0 + mu_b1*con_pred$log_mass_ct

# Exponentiate prediction - note it's mass-specific
con_pred$y_spec <- exp(con_pred$log_y_spec)

# Whole organism rate
con_pred$y <- con_pred$y * con_pred$mass_g

# Specify unit
con_pred$y_g_d <- con_pred$y


#**** Net gain =====================================================================
# Combine data
met_pred_sub <- met_pred %>% select(y_g_d, mass_g, rate)
con_pred_sub <- con_pred %>% select(y_g_d, mass_g, rate)

diff <- con_pred_sub

diff$rate <- "Net gain"

diff$y_g_d <- con_pred_sub$y_g_d - met_pred_sub$y_g_d

diff$y_g_g_d2 <- con_pred_sub$y_g_d/con_pred_sub$mass_g - met_pred_sub$y_g_d/met_pred_sub$mass_g

# Now calculate the slopes
# Net gain vs mass (log-log) in unit g per day
summary(lm(log(diff$y_g_d) ~ log(diff$mass_g)))
# mass-slope = 6.243e-01 = 0.6243
# The mass-specific slope is therefore 0.6243-1 = -0.3757
diff$y_g_g_d <- diff$y_g_d / diff$mass_g
summary(lm(log(diff$y_g_g_d) ~ log(diff$mass_g)))
# -3.757e-01 = -0.3757


