#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.04.29: Max Lindmark
#
# - Code to exemplify the different temperature dependencies of metabolism and 
#   conssumption using already fitted models. Note that the models for metabolism
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

#** Use log-log models to calculate temperature-independent rates ==================
#**** Metabolism ===================================================================
# This is what the model is fitted to:
# data = list(
#   y = log(dat$y_spec),
#   n_obs = length(dat$y),
#   species_n = dat$species_n,
#   mass = dat$log_mass_ct,
#   temp = dat$temp_arr_ct)

# Set temperature & mass range. 
temp_c_ct <- seq(-20, 40, 0.1)
mass_range <- c(2, 200)

# This is the mass-range in the respective datasets
# > summary(met$mass_g)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.025     9.545    85.530   310.779   238.818 21340.920 
# > summary(con$mass_g)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0515    1.1794   25.9101  103.5015  161.6000 1167.8281 

# Create new prediction data frame
met_pred <- data.frame(expand.grid(mass_g = mass_range,
                                   temp_c_ct = temp_c_ct))

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

# Add normalized temperature scalars 
# Here we want scalars that collapse to 1 at the mean temperature
# For metabolism this is exponential, for consumption it is based on the polynomial model

# First add non-normalized temperatures
met_pred$temp_c <- met_pred$temp_c_ct + 19

# Calculate Kelvin temperatures
met_pred$temp_K <- met_pred$temp_c + 273.15

# Calculate temp scalar
K <- 8.617332e-05 # Boltzmann's constant
met_pred$temp_scalar <- exp(mu_b2/(K * met_pred$temp_K))
  
# Plot temp scalar
ggplot(met_pred, aes(temp_c,temp_scalar)) + geom_point()

# Rescale to equal 1 at 19C
scaling_factor_m <- filter(met_pred, temp_c == 19)$temp_scalar[1]

met_pred$temp_scalar_scaled <- met_pred$temp_scalar / scaling_factor_m

# Plot rescaled temp scalar
ggplot(met_pred, aes(temp_c, temp_scalar_scaled)) +
  geom_point() +
  geom_vline(xintercept = 19) +
  geom_hline(yintercept = 1)

# Add temperature-dependent rates
met_pred$y_g_d_temp <- met_pred$y_g_d * met_pred$temp_scalar_scaled

# Plot temperature-corrected metabolic rate
ggplot(met_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) + geom_line()
#unique(met_pred$y_g_d)


#**** Consumption ==================================================================
# Create new prediction data frame
con_pred <- data.frame(expand.grid(mass_g = mass_range,
                                   temp_c_ct = temp_c_ct))

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

# Add normalized temperature scalars 
# Here we want scalars that collapse to 1 at the mean temperature
# For metabolism this is exponential, for consumption it is based on the polynomial model

# First add non-normalized temperatures
con_pred$temp_c <- con_pred$temp_c_ct + 19

# Coefficients from polynomial model
mu_b0 <- 0.81
b1 <- -0.0016
b2 <- 0.042
b3 <- -0.00025

# Now calculate the scalar by predicting maximum consumption for the mean body size (i.e. 0)
# The model is fitted to environment-centered temperature and mean-centered mass
con_pred$temp_scalar <- mu_b0 + b1*0 + b2*con_pred$temp_c_ct + b3*con_pred$temp_c_ct^2

# Plot temp scalar
ggplot(con_pred, aes(temp_c, temp_scalar)) + geom_point()

# Rescale to equal 1 at 19C
scaling_factor_c <- filter(con_pred, temp_c == 19)$temp_scalar[1]

con_pred$temp_scalar_scaled <- con_pred$temp_scalar * (1/scaling_factor_c)

# Plot rescaled temp scalar
ggplot(con_pred, aes(temp_c, temp_scalar_scaled)) + 
  geom_point() +
  geom_vline(xintercept = 19) +
  geom_hline(yintercept = 1)

# Add temperature-dependent rates
con_pred$y_g_d_temp <- con_pred$y_g_d * con_pred$temp_scalar_scaled

# Plot temperature-corrected maximum consumption rate
ggplot(con_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) + geom_line()


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


# Combine data
met_pred_sub <- met_pred %>% select(y_g_d_temp, mass_g, temp_c_ct, temp_c, rate) %>% arrange(temp_c_ct)
con_pred_sub <- con_pred %>% select(y_g_d_temp, mass_g, temp_c_ct, temp_c, rate) %>% arrange(temp_c_ct)

# head(met_pred_sub)
# head(con_pred_sub)

diff <- con_pred_sub

diff$rate <- "Net gain"

diff$y_g_d_temp <- con_pred_sub$y_g_d_temp - met_pred_sub$y_g_d_temp

dat <- rbind(met_pred_sub, con_pred_sub, diff)

# Add variable for plotting
dat$mass_facet <- as.factor(dat$mass_g)
levels(dat$mass_facet) <- c("2 g", "200 g")

# Find peak temperature
peak <- dat %>% 
  group_by(mass_g) %>% 
  filter(rate == "Net gain") %>% 
  filter(y_g_d_temp == max(y_g_d_temp))

# Set palette
pal <- RColorBrewer::brewer.pal("Dark2", n = 3)

# Add a dummy for axis range (for facets individually...)
dummy <- data.frame(y_g_d_temp = c(0.5, 2), mass_g = c(2, 200), rate = "Maximum consumption",
                    mass_facet = c("2 g", "200 g"), temp_c_ct = 15)

met_ribbon <- filter(dat, rate == "Metabolism")
colnames(met_ribbon)[1] <- "Metabolic_rate"
con_ribbon <- filter(dat, rate == "Maximum consumption")
colnames(con_ribbon)[1] <- "Maximum_consumption_rate"

con_ribbon <- con_ribbon %>%
  select(Maximum_consumption_rate)

ribbon <- cbind(met_ribbon, con_ribbon)
ribbon$diff <- ribbon$Maximum_consumption_rate - ribbon$Metabolic_rate
ribbon <- ribbon %>% filter(diff > 0.02)
ribbon <- ribbon %>% filter(Metabolic_rate > 0 & Maximum_consumption_rate > 0 & temp_c_ct > -5 & temp_c_ct < 30)


# Plot all 6 curves
p1 <- dat %>% filter(y_g_d_temp > 0 & temp_c_ct > -5 & temp_c_ct < 30) %>% 
  ggplot(., aes(temp_c_ct, y_g_d_temp, color = factor(rate),
                linetype = factor(mass_g), alpha = factor(rate))) +
  geom_ribbon(data = ribbon, inherit.aes = FALSE,
              aes(x = temp_c_ct, ymin = Metabolic_rate, ymax = Maximum_consumption_rate),
              fill = "grey95", color = NA) +
  geom_line(size = 0.75) +
  coord_cartesian(expand = 0) + 
  scale_color_manual(values = pal) +
  scale_linetype_manual(values = c(1, 1)) +
  scale_alpha_manual(values = c(0.5, 0.5, 1)) +
  #scale_alpha_manual(values = c(1, 1, 1)) +
  labs(x = expression(paste("Rescaled temperature [", degree*C, "]")),
       y = "Rescaled rates [g/day]",
       color = "Rate") +
  geom_segment(data = peak, aes(x = temp_c_ct, xend = temp_c_ct, y = y_g_d_temp, yend = 0),
               size = 0.75, arrow = arrow(length = unit(0.25, "cm")), show.legend = FALSE) +
  guides(linetype = FALSE) +
  guides(alpha = FALSE) +
  guides(color = guide_legend(ncol = 1)) +
  facet_wrap(~ mass_facet, scales = "free_y", ncol = 1) + # new addition that is needed for the larger size range
  geom_point(data = dummy, aes(temp_c_ct, y_g_d_temp), inherit.aes = F, color = "white", fill = "white") + 
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 8), # 12
                                       aspect.ratio = 1,
                                       legend.title = element_text(size = 6), # 10
                                       legend.key.size = unit(0.2, 'cm'), 
                                       legend.text = element_text(size = 6), 
                                       #legend.position = c(0.32, 0.91),
                                       legend.position = "bottom",
                                       legend.background = element_rect(fill = NA))

pWord1

#ggsave("figures/concept.png", width = 3.5, height = 6.5, dpi = 600)
ggsave("figures/concept.png", width = 11, height = 11, dpi = 600, unit = "cm")
