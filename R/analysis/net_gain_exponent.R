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
mu_b0 <- -2.1672409
mu_b1 <- -0.2039408
mu_b2 <- -0.6118622

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
mu_b1 <- -0.3750
mu_b2 <- -0.6945

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

summary(lm(log(diff$y_g_d) ~ log(diff$mass_g)))
summary(lm(log(con_pred_sub$y_g_d) ~ log(con_pred_sub$mass_g)))
summary(lm(log(met_pred_sub$y_g_d) ~ log(met_pred_sub$mass_g)))

dat <- rbind(met_pred_sub, con_pred_sub, diff)

ggplot(dat, aes(log(mass_g), log(y_g_d), color = rate)) +
  geom_line(size = 1.5) +
  scale_color_brewer(palette = "Dark2") + 
  theme_classic() + 
  theme(aspect.ratio = 1) +
  geom_abline(intercept = -2, slope = 0.64, linetype = 2, color = "grey")
  NULL

