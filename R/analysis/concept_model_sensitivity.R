#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2020.04.29: Max Lindmark
#
# - Code to exemplify the different temperature dependencies of metabolism and 
#   consumption using allometric relationships and the effect on net energy gain
#   Note that the models for metabolism and (non-linear) consumption use different
#   predictor and response variables (consumption uses normalized consumption ~ mass_g
#   and centered temperature in Celsius, whereas metabolism is modeled as
#   log(y) ~ log(centered mass) + Arrhenius temperature)
#
#   I have therefore plotted normalized rates (0:1 within each rate), to illustrate 
#   the curves
# 
#   In this script, I in addition use posterior samples to populate the vector of
#   mass-exponents, to visualize uncertainty in the prediction around the T_opt 
#   declining.
# 
# A. Load libraries
#
# B. Generate data
#
# C. Plot
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

set.seed(42)

# B. GENERATE DATA =================================================================
# First read in data so that we can take the means etc

met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Now read in the posterior samples (n=1000) that I save in the meta_cons_model_pred.R script
# https://github.com/maxlindmark/scaling/blob/master/R/analysis/log_linear/meta_cons_model_pred.R#L318-L336

met_exponent_samples <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/output/meta_exponent_samples.csv"))
met_exponent_samples_full <- met_exponent_samples
median_met_mub1 <- median(met_exponent_samples_full$value)

met_exponent_samples <- met_exponent_samples[sample(nrow(met_exponent_samples), 50), ]
median(met_exponent_samples$value)


con_exponent_samples <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/output/cons_exponent_samples.csv"))

con_exponent_samples_full <- con_exponent_samples
median_con_mub1 <- median(con_exponent_samples_full$value)

con_exponent_samples <- con_exponent_samples[sample(nrow(con_exponent_samples), 50), ]
median(con_exponent_samples$value)

#** Use log-log models to calculate temperature-independent rates ==================
#**** Metabolism ===================================================================
# Set temperature & mass range. 
temp_c_ct <- seq(-20, 20, 0.001)
mass_range <- c(5, 1000)

# This is the mass-range in the respective datasets
# > summary(met$mass_g)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.025     9.399    85.135   310.324   238.558 21340.920 
# > summary(con$mass_g)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0515    1.1794   25.9101  103.5015  161.6000 1167.8281 

# Create new prediction data frame
met_pred <- data.frame(expand.grid(mass_g = mass_range,
                                   mu_b1 = c(met_exponent_samples$value, median_met_mub1), # We only use the posterior for the mass-exponent - rest are kept at their medians
                                   temp_c_ct = temp_c_ct))
                 
# Add which rate it is
met_pred$rate <- "Metabolism"

# Calculate log mass
met_pred$log_mass <- log(met_pred$mass_g)

# Calculate mean-centered log mass (as it's used for model fitting)
met_pred$log_mass_ct <- met_pred$log_mass - mean(met$log_mass)

# Calculate log metabolic rate
# Regression coefficients: see meta_cons_analysis.R. Ignore the mass-temp interaction(!)
#              2.5%       25%        50%      75%     97.5%
# mu_b0_r     1.681766  1.7997409  1.8598459  1.920484  2.04625
# mu_b0_s     0.974255  1.1908261  1.2925437  1.395646  1.61421
# mu_b1       0.741257  0.7749627  0.7909194  0.807061  0.83924
# mu_b2      -0.671663 -0.6374733 -0.6209867 -0.604387 -0.57241

# We only use the posterior for the mass-exponent - rest are kept at their medians
mu_b0_r <- 1.8598459
# mu_b1 # this is a column in the data frame already with 1000 posterior samples
mu_b2 <- -0.6209867

# Temperature-independent allometric function (valid at 19C = mean temperature)
met_pred$log_y <- mu_b0_r + met_pred$mu_b1*met_pred$log_mass_ct

log_y_intercept <- mu_b0_r + median_met_mub1*0 # *(log(0) - mean(met$log_mass)) #### Calculate also for 1 g mass to get the value in unit g/d

# Exponentiate prediction
met_pred$y <- exp(met_pred$log_y)
y_intercept <- exp(log_y_intercept) #### these hashtags are for calculating the interept, no the conversion factor for a specific size

# Now convert to g/d using the same values as Jan in Ohlberger et al (2012) Oikos
# 1 kcal = 295 mg O2 
# 0.003389831 kcal = 1 mg 02
met_pred$y_kcal_h <- met_pred$y * 0.003389831
y_intercept_kcal_h <- y_intercept * 0.003389831 ####

# Convert to cal/h from kilo cal/h 
met_pred$y_cal_h <- met_pred$y_kcal_h * 1000
y_intercept_cal_h <- y_intercept_kcal_h * 1000 ####

# Convert to Joule/h
met_pred$y_J_h <- met_pred$y_cal_h * 4.1855
y_intercept_J_h <- y_intercept_cal_h * 4.1855 #### 

# Convert to g/h. Energy content of fish: 5600 J/g
met_pred$y_g_h <- met_pred$y_J_h / 5600
y_intercept_g_h <- y_intercept_J_h / 5600 #### 

# Convert to g/d
met_pred$y_g_d <- met_pred$y_g_h * 24
y_intercept_g_d <- y_intercept_g_h * 24 #### 

# Test
exp(1.8598459)*(0.003389831*1000*4.1855/5600)*24
y_intercept_g_d
met_pred %>%
  filter(mass_g == min(mass_g) & temp_c_ct == 0 & mu_b1 == median_met_mub1) %>%
  dplyr::select(y_g_d)

# Now that I have metabolism in unit g_d, check the conversion factor
met_pred %>%
  filter(temp_c_ct == 0 & mu_b1 == median_met_mub1) %>% 
  mutate(conv_factor = y_g_d / y)
  
# Add non-normalized temperatures
met_pred$temp_c <- met_pred$temp_c_ct + 19

# Calculate Kelvin temperatures
met_pred$temp_K <- met_pred$temp_c + 273.15

# Calculate temp scalar
K <- 8.617332e-05 # Boltzmann's constant
met_pred$temp_scalar <- exp(mu_b2/(K * met_pred$temp_K))

# Plot temp scalar
# ggplot(met_pred, aes(temp_c, temp_scalar)) + geom_point()

# Rescale to max
met_pred$temp_scalar_scaled <- met_pred$temp_scalar / max(met_pred$temp_scalar)

# Plot rescaled temp scalar
# ggplot(met_pred, aes(temp_c, temp_scalar_scaled)) +
#   geom_point() +
#   geom_vline(xintercept = 19) +
#   geom_hline(yintercept = 1)

# Make allometric rates temperature-dependent
met_pred$y_g_d_temp <- met_pred$y_g_d * met_pred$temp_scalar_scaled

# Plot "temperature-corrected" metabolic rate
# ggplot(met_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) + geom_line()

# Now we need to scale it back again so that the temperature-dependent rate equals
# the non-temperature-dependent rate at 19C. Find the scaling factor that does that.
met_factor <- filter(met_pred, temp_c_ct == 0)$y_g_d / 
  filter(met_pred, temp_c_ct == 0)$y_g_d_temp

met_pred %>% 
  filter(temp_c == 19) %>% 
  ggplot(., aes(y_g_d , y_g_d_temp*met_factor[1])) +
  geom_line(size = 3) + 
  geom_abline(color = "red", linetype = 2, size = 3) + 
  facet_wrap(~ mass_g, scales = "free") 

# Here, the horizontal lines correspond to the values at 19C without any temperature,
# and we see that they intersect the temperature-dependent curves at 19C
met_pred$y_g_d_temp <- met_pred$y_g_d_temp*met_factor[1]

met_pred %>% filter(temp_c_ct == 0)

# Now allometric metabolic rates are multiplied by temperature-factors that equal 1
# at the maximum temperature in the data frame, and then the whole rate is shifted 
# with another factor so that the whole function at 19C equals the non-temperature
# dependent metabolic rate
# Plot again without the correction factor:
# ggplot(met_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) +
#   geom_line() + 
#   geom_hline(yintercept = filter(met_pred, temp_c_ct == 0)$y_g_d) +
#   geom_vline(xintercept = 19) + 
#   coord_cartesian(ylim = c(0, 5))


#**** Consumption ==================================================================
# Create new prediction data frame
con_pred <- data.frame(expand.grid(mass_g = mass_range,
                                   mu_b1 = c(con_exponent_samples$value, median_con_mub1), # We only use the posterior for the mass-exponent - rest are kept at their medians
                                   temp_c_ct = temp_c_ct))

# Add which rate it is
con_pred$rate <- "Maximum consumption"

# Calculate log mass
con_pred$log_mass <- log(con_pred$mass_g)

# Calculate mean-centered mass (as its used for model fitting)
con_pred$log_mass_ct <- con_pred$log_mass - mean(con$log_mass)

# Calculate log maximum consumption rate
# Regression coefficients: see meta_cons_analysis.R
#           2.5%     25%       50%       75%       97.5%
# ...
# mu_b0    -0.85268 -0.511181 -0.34097 -0.17422  0.16574
# mu_b1     0.54904  0.600763  0.62663  0.65332  0.70743
# mu_b2    -0.84671 -0.743781 -0.69353 -0.64283 -0.53869

mu_b0 <- -0.34097
# mu_b1 # this is a column in the data frame already with 1000 posterior samples
mu_b2 <- -0.69353

# Temperature-independent allometric function (valid at 19C = mean temperature)
#con_pred$log_y_spec <- mu_b0 + mu_b1*con_pred$log_mass_ct
con_pred$log_y <- -0.34097 + con_pred$mu_b1*con_pred$log_mass_ct

# Exponentiate prediction
con_pred$y <- exp(con_pred$log_y)

# Specify unit
con_pred$y_g_d <- con_pred$y

# Add non-normalized temperatures
con_pred$temp_c <- con_pred$temp_c_ct + 19

# Calculate temp scalar
# Coefficients from Sharpe-Schoolfield model

# 2. Quantiles for each variable:
#   
#           2.5%   25%    50%    75%    97.5%
# Eh        1.6753 1.8099 1.8849 1.9578 2.0987
# Th       -0.8571 0.1837 0.7463 1.2769 2.3665
# mu_E      0.5353 0.6591 0.7275 0.7962 0.9408
# mu_b0     0.5769 0.7227 0.7880 0.8521 0.9908

# Median
Eh <- 1.8849
Th <- 0.7463
mu_E <- 0.7275
mu_b0 <- 0.7880

# Define the constants in the Sharpe-Schoolfield equations
tref <- -10
bk <- 8.62e-05

# Now calculate the scalar by predicting maximum consumption for the mean body size (i.e. 0)
# The model is fitted to environment-centered temperature and mean-centered mass
con_pred$temp_scalar <- (mu_b0 * exp(mu_E*((1/(bk*(tref + 273.15))) - 1/(bk*(con_pred$temp_c_ct + 273.15))))) / (1 + exp(Eh*((1/(bk*(Th + 273.15))) - (1/(bk*(con_pred$temp_c_ct + 273.15))))))

# Plot temp scalar
# ggplot(con_pred, aes(temp_c_ct, temp_scalar)) + geom_point()

# Rescale to max
con_pred$temp_scalar_scaled <- con_pred$temp_scalar / max(con_pred$temp_scalar)

# Plot rescaled temp scalar
# ggplot(con_pred, aes(temp_c, temp_scalar_scaled)) +
#   geom_point() +
#   geom_vline(xintercept = 19) +
#   geom_hline(yintercept = 1)

# Make allometric rates temperature-dependent
con_pred$y_g_d_temp <- con_pred$y_g_d * con_pred$temp_scalar_scaled

# Plot temperature-corrected maximum consumption rate
ggplot(con_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) + geom_line()

# Now we need to scale it back again so that the temperature-dependent rate equals
# the non-temperature-dependent rate at 19C. Find the scaling factor that does that.
con_factor <- filter(con_pred, temp_c_ct == 0)$y_g_d /
  filter(con_pred, temp_c_ct == 0)$y_g_d_temp

# Same for both sizes because there's no mass-temp interaction here
con_pred %>% 
  filter(temp_c == 19) %>% 
  ggplot(., aes(y_g_d , y_g_d_temp*con_factor[1])) +
  geom_line(size = 3) + 
  geom_abline(color = "red", linetype = 2, size = 3) + 
  facet_wrap(~ mass_g, scales = "free") 

# Here, the vertical lines correspond to the values at 19C without any temperature,
# and we see that they intersect the temperature-dependent curves at 19C
con_pred$y_g_d_temp <- con_pred$y_g_d_temp*con_factor[1]

con_pred %>% filter(temp_c_ct == 0)

# Now allometric consumption rates are multiplied by temperature-factors that equal 1
# at the maximum temperature in the data frame, and then the whole rate is shifted 
# with another factor so that the whole function at 19C equals the non-temperature
# dependent consumption rate
# Plot again without the correction factor:
# ggplot(con_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) +
#   geom_line() + 
#   geom_hline(yintercept = filter(con_pred, temp_c_ct == 0)$y_g_d) +
#   geom_vline(xintercept = 19) + 
#   coord_cartesian(ylim = c(0, 12))

# Now that we have both rates, we can calculate the difference and make the conceptual figure:


# C. PLOT ==========================================================================
#**** Combine data =================================================================
met_pred_sub <- met_pred %>%
  dplyr::select(y_g_d_temp, mass_g, temp_c_ct, temp_c, rate, mu_b1) %>%
  arrange(temp_c_ct)

con_pred_sub <- con_pred %>%
  dplyr::select(y_g_d_temp, mass_g, temp_c_ct, temp_c, rate, mu_b1) %>%
  arrange(temp_c_ct)

diff <- con_pred_sub

diff$rate <- "Net gain"

diff$y_g_d_temp <- con_pred_sub$y_g_d_temp - met_pred_sub$y_g_d_temp

dat <- rbind(met_pred_sub, con_pred_sub, diff)
dat %>% group_by(rate) %>% summarise(max(y_g_d_temp))

# Add variable for plotting
dat$mass_facet <- as.factor(dat$mass_g)
levels(dat$mass_facet) <- c("5 g", "1000 g")

# Find peak temperature
peak <- dat %>% 
  group_by(mass_g, mu_b1) %>% 
  filter(rate == "Net gain") %>% 
  filter(y_g_d_temp == max(y_g_d_temp))

peak %>%
  filter(y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max) %>% 
  ggplot(., aes(temp_c_ct, fill = mass_facet)) + geom_histogram()


# How many unique peak temps per mass and mu_b1?
#peak %>% group_by(mu_b1, mass_g) %>% mutate(n = n()) %>% ungroup() %>% distinct(n)

ggplot(filter(peak, mass_g == 5), aes(temp_c_ct)) + geom_density()
ggplot(filter(peak, mass_g == 5), aes(temp_c_ct)) + geom_histogram()

# Set palettes
dark2 <- RColorBrewer::brewer.pal("Dark2", n = 8)
pal <- dark2[1:3]
pal2 <- dark2[7:8]

# Set plotting range
t_max <- 3.1
t_min <- -5.5

# Filter data to only show positive rates
pdat <- dat %>% filter(y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max)

# Make the plots!
p1a <- pdat %>%
  filter(mass_g == 5) %>% 
  ggplot(., aes(temp_c_ct, y_g_d_temp, color = factor(rate), linetype = factor(mu_b1))) +
  guides(linetype = FALSE, color = guide_legend(ncol = 1, override.aes = list(alpha = 0.8))) +
  geom_line(alpha = 0.15) +
  geom_line(data = filter(pdat, mass_g == 5 & mu_b1 %in% c(median_met_mub1, median_con_mub1)), # bad code, luckily they are different enough
            aes(temp_c_ct, y_g_d_temp, color = factor(rate)),
            size = 1, alpha = 1, linetype = 1) + 
  coord_cartesian(expand = 0, ylim = c(0, 0.42)) + 
  geom_rug(data = filter(peak, y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max & mass_g == 5),
           aes(x = temp_c_ct, y = y_g_d_temp), alpha = 0.3, sides = "b",
           length = unit(0.03, "npc"), color = pal2[1]) +
  scale_color_manual(values = pal) +
  scale_linetype_manual(values = c(rep(1, length(unique(pdat$mu_b1))))) + 
  labs(x = "",
       y = "Rescaled rates",
       color = "") +
  ggtitle("5 g") +
  theme_classic() +
  theme(text = element_text(size = 8),
        plot.title = element_text(size = 8),
        aspect.ratio = 1,
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, 'cm'), 
        legend.text = element_text(size = 7), 
        legend.margin = margin(-0.3, 0, 0, 0, unit="cm"),
        legend.position = c(0.45, 0.35),
        legend.background = element_rect(fill = NA),
        plot.tag = element_text(face = 'bold'))

p1b <- pdat %>%
  filter(mass_g == 1000) %>% 
  ggplot(., aes(temp_c_ct, y_g_d_temp, color = factor(rate), linetype = factor(mu_b1))) +
  geom_line(alpha = 0.15) +
  geom_line(data = filter(pdat, mass_g == 1000 & mu_b1 %in% c(median_met_mub1, median_con_mub1)), # bad code, luckily they are different enough
            aes(temp_c_ct, y_g_d_temp, color = factor(rate)),
            size = 1, alpha = 1, linetype = 1) + 
  coord_cartesian(expand = 0, ylim = c(0, 17)) + 
  geom_rug(data = filter(peak, y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max & mass_g == 1000),
           aes(x = temp_c_ct, y = y_g_d_temp), alpha = 0.3, sides = "b", 
           length = unit(0.03, "npc"), color = pal2[2]) +
  scale_color_manual(values = pal) +
  scale_linetype_manual(values = c(rep(1, length(unique(pdat$mu_b1))))) + 
  labs(x = expression(paste("Rescaled temperature [", degree*C, "]")),
       y = "Rescaled rates",
       color = "") +
  guides(linetype = FALSE, color = FALSE) +
  ggtitle("1000 g") +
  theme_classic() + 
  theme(text = element_text(size = 8), 
        plot.title = element_text(size = 8),
        aspect.ratio = 1,
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, 'cm'), 
        legend.text = element_text(size = 7), 
        legend.margin = margin(-0.3, 0, 0, 0, unit="cm"),
        legend.position = "bottom",
        legend.background = element_rect(fill = NA),
        plot.tag = element_text(face = 'bold'))

p2 <- peak %>%
  filter(y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max) %>% 
  ggplot(., aes(temp_c_ct, fill = mass_facet, color = mass_facet)) +
  geom_histogram(aes(y=..density..), color = NA, alpha = 0.4, bins = 18) + 
  geom_density(alpha = 0.4) +
  theme_classic() +
  coord_cartesian(expand = 0) +
  scale_fill_manual(values = pal2, name = "") +
  scale_color_manual(values = pal2, name = "") +
  ggtitle("Optimum growth temperature") + 
  labs(x = "",
       y = "Density") +
  theme_classic() +
  theme(text = element_text(size = 8), 
        plot.title = element_text(size = 8),
        aspect.ratio = 1,
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, 'cm'), 
        legend.text = element_text(size = 7), 
        legend.margin = margin(-0.3, 0, 0, 0, unit = "cm"),
        legend.position = c(0.2, 0.9),
        legend.background = element_rect(fill = NA),
        plot.tag = element_text(face = 'bold'))

(p1a | p1b | p2) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("figures/concept_sensitivity.png", width = 15, height = 15, dpi = 600, unit = "cm")

