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
# Set temperature & mass range. 
temp_c_ct <- seq(-20, 20, 0.1)
mass_range <- c(5, 1000)

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

# Temperature-independent allometric function (valid at 19C = mean temperature)
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

# Add non-normalized temperatures
met_pred$temp_c <- met_pred$temp_c_ct + 19

# Calculate Kelvin temperatures
met_pred$temp_K <- met_pred$temp_c + 273.15

# Calculate temp scalar
K <- 8.617332e-05 # Boltzmann's constant
met_pred$temp_scalar <- exp(mu_b2/(K * met_pred$temp_K))

# Plot temp scalar
ggplot(met_pred, aes(temp_c, temp_scalar)) + geom_point()

# Rescale to max
met_pred$temp_scalar_scaled <- met_pred$temp_scalar / max(met_pred$temp_scalar)

# Plot rescaled temp scalar
ggplot(met_pred, aes(temp_c, temp_scalar_scaled)) +
  geom_point() +
  geom_vline(xintercept = 19) +
  geom_hline(yintercept = 1)

# Add temperature-dependent rates
met_pred$y_g_d_temp <- met_pred$y_g_d * met_pred$temp_scalar_scaled

# Plot "temperature-corrected" metabolic rate
ggplot(met_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) + geom_line()

# Now we need to scale it back again so that the temperature-dependent rate equals
# the non-temperature-dependent rate at 19C. Find the scaling factor that does that.
met_factor <- filter(met_pred, temp_c_ct == 0)$y_g_d / 
  filter(met_pred, temp_c_ct == 0)$y_g_d_temp

ggplot(met_pred, aes(temp_c, y_g_d_temp*met_factor[1], color = factor(mass_g))) +
  geom_line() + 
  geom_hline(yintercept = filter(met_pred, temp_c_ct == 0)$y_g_d) +
  geom_vline(xintercept = 19) + 
  coord_cartesian(ylim = c(0, 1.2))

# Here, the vertical lines correspond to the values at 19C without any temperature,
# and we see that they intersect the temperature-dependent curves at 19C

met_pred$y_g_d_temp <- met_pred$y_g_d_temp*met_factor[1]

met_pred %>% filter(temp_c_ct == 0)

# Now allometric metabolic rates are multiplied by temperature-factors that equal 1
# at the maximum temperature in the data frame, and then the whole rate is shifted 
# with another factor so that the whole function at 19C equals the non-temperature
# dependent metabolic rate


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

# Temperature-independent allometric function (valid at 19C = mean temperature)
con_pred$log_y_spec <- mu_b0 + mu_b1*con_pred$log_mass_ct

# Exponentiate prediction - note it's mass-specific
con_pred$y_spec <- exp(con_pred$log_y_spec)

# Whole organism rate
con_pred$y <- con_pred$y * con_pred$mass_g

# Specify unit
con_pred$y_g_d <- con_pred$y

# Add non-normalized temperatures
con_pred$temp_c <- con_pred$temp_c_ct + 19

# Calculate temp scalar
# Coefficients from Sharpe-Schoolfield model

# With peak_temp_ct
# 2. Quantiles for each variable:
#   
#          2.5%    25%    50%    75%  97.5%
# Eh       2.1712 2.4675 2.6374 2.8278 3.2211
# Th       2.8158 3.6412 4.0288 4.3708 4.9798
# mu_E     0.4454 0.5308 0.5770 0.6254 0.7359
# mu_b0    0.5204 0.6463 0.7029 0.7611 0.8853
# sigma_b0 0.1498 0.2046 0.2463 0.3012 0.4666

# Median
Eh <- 2.6374
Th <- 4.0288
mu_E <- 0.5770
mu_b0 <- 0.7029

# Add in constants as well
tref <- -10
bk <- 8.62e-05

# Now calculate the scalar by predicting maximum consumption for the mean body size (i.e. 0)
# The model is fitted to environment-centered temperature and mean-centered mass
con_pred$temp_scalar <- (mu_b0 * exp(mu_E*((1/(bk*(tref + 273.15))) - 1/(bk*(con_pred$temp_c_ct + 273.15))))) / (1 + exp(Eh*((1/(bk*(Th + 273.15))) - (1/(bk*(con_pred$temp_c_ct + 273.15))))))

# Plot temp scalar
ggplot(con_pred, aes(temp_c_ct, temp_scalar)) + geom_point()

# Rescale to max
con_pred$temp_scalar_scaled <- con_pred$temp_scalar / max(con_pred$temp_scalar)

# Plot rescaled temp scalar
ggplot(con_pred, aes(temp_c, temp_scalar_scaled)) +
  geom_point() +
  geom_vline(xintercept = 19) +
  geom_hline(yintercept = 1)

# Add temperature-dependent rates
con_pred$y_g_d_temp <- con_pred$y_g_d * con_pred$temp_scalar_scaled

# Plot temperature-corrected maximum consumption rate
ggplot(con_pred, aes(temp_c, y_g_d_temp, color = factor(mass_g))) + geom_line()

# Now we need to scale it back again so that the temperature-dependent rate equals
# the non-temperature-dependent rate at 19C. Find the scaling factor that does that.
con_factor <- filter(con_pred, temp_c_ct == 0)$y_g_d /
  filter(con_pred, temp_c_ct == 0)$y_g_d_temp

ggplot(con_pred, aes(temp_c, y_g_d_temp*con_factor[1], color = factor(mass_g))) +
  geom_line() + 
  geom_hline(yintercept = filter(con_pred, temp_c_ct == 0)$y_g_d) +
  geom_vline(xintercept = 19) + 
  coord_cartesian(ylim = c(0, 5))

# Here, the vertical lines correspond to the values at 19C without any temperature,
# and we see that they intersect the temperature-dependent curves at 19C

con_pred$y_g_d_temp <- con_pred$y_g_d_temp*con_factor[1]

con_pred %>% filter(temp_c_ct == 0)

# Now allometric consumption rates are multiplied by temperature-factors that equal 1
# at the maximum temperature in the data frame, and then the whole rate is shifted 
# with another factor so that the whole function at 19C equals the non-temperature
# dependent consumption rate

# Now that we have both rates, we can calculate the difference and make the conceptual figure:

  
#*** Plot ==========================================================================
#**** Combine data =================================================================
met_pred_sub <- met_pred %>%
  dplyr::select(y_g_d_temp, mass_g, temp_c_ct, temp_c, rate) %>%
  arrange(temp_c_ct)

con_pred_sub <- con_pred %>%
  dplyr::select(y_g_d_temp, mass_g, temp_c_ct, temp_c, rate) %>%
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
  group_by(mass_g) %>% 
  filter(rate == "Net gain") %>% 
  filter(y_g_d_temp == max(y_g_d_temp))

# Set palette
pal <- RColorBrewer::brewer.pal("Dark2", n = 3)

# Set plotting range
t_max <- 3.1
t_min <- -5.5

# Add a dummy for axis range (for facets individually so that the y-axis range can be 5% more than maximum for visibility)
m5_rate <- max(filter(dat, y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max & mass_g == 5)$y_g_d_temp)

m1000_rate <- max(filter(dat, y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max & mass_g == 1000)$y_g_d_temp)  

dummy <- data.frame(y_g_d_temp = c(m5_rate*1.05, m1000_rate*1.05), mass_g = c(5, 1000), rate = "Maximum consumption",
                    mass_facet = c("5 g", "1000 g"), temp_c_ct = 0)

met_ribbon <- filter(dat, rate == "Metabolism")
colnames(met_ribbon)[1] <- "Metabolic_rate"
net_ribbon <- filter(dat, rate == "Net gain")
colnames(net_ribbon)[1] <- "Net_gain"

net_ribbon <- net_ribbon %>% select(Net_gain)

ribbon <- cbind(met_ribbon, net_ribbon)
ribbon$diff <- ribbon$Net_gain - ribbon$Metabolic_rate
#ribbon <- ribbon %>% filter(diff > 0.02)
ribbon <- ribbon %>%
  filter(Metabolic_rate > 0 & Net_gain > 0 & temp_c_ct > t_min & temp_c_ct < t_max)

# Plot all 6 curves
p1 <- dat %>% filter(y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max) %>% 
  ggplot(., aes(temp_c_ct, y_g_d_temp, color = factor(rate),
                linetype = factor(mass_g), alpha = factor(rate))) +
  geom_ribbon(data = filter(ribbon, temp_c_ct > t_min & temp_c_ct < t_max & diff > 0),
              inherit.aes = FALSE,
              aes(x = temp_c_ct, ymin = Metabolic_rate, ymax = Net_gain),
              fill = "grey95", color = NA) +
  geom_line(size = 0.5) +
  coord_cartesian(expand = 0) + 
  scale_color_manual(values = pal) +
  scale_linetype_manual(values = c(1, 1)) +
  scale_alpha_manual(values = c(0.5, 0.5, 1)) +
  labs(x = expression(paste("Rescaled temperature [", degree*C, "]")),
       y = "Rescaled rates [g/day]",
       color = "Rate") +
  geom_segment(data = filter(peak, y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max),
               aes(x = temp_c_ct, xend = temp_c_ct, y = y_g_d_temp, yend = 0),
               size = 0.5, arrow = arrow(length = unit(0.25, "cm")), show.legend = FALSE) +
  guides(linetype = FALSE) +
  guides(alpha = FALSE) +
  guides(color = guide_legend(ncol = 1)) +
  facet_wrap(~ mass_facet, scales = "free_y", ncol = 1) + # new addition that is needed for the larger size range
  geom_point(data = filter(dummy, y_g_d_temp > 0 & temp_c_ct > t_min & temp_c_ct < t_max),
             aes(temp_c_ct, y_g_d_temp), inherit.aes = F, color = "white", fill = "white")

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 8), # 12
                                       aspect.ratio = 1,
                                       legend.title = element_text(size = 6), # 10
                                       legend.key.size = unit(0.2, 'cm'), 
                                       legend.text = element_text(size = 6), 
                                       legend.position = c(0.85, 0.55),
                                       #legend.position = "bottom",
                                       legend.background = element_rect(fill = NA))

pWord1

ggsave("figures/concept.png", width = 11, height = 11, dpi = 600, unit = "cm")

# Just a quick calculation here...
# The difference in peak temperature here is 1.1C between a 5 g and 1000 g fish
# Our empirical analysis of T_opt gives this relationship:
# T_opt=-0.074-0.31×m 
# where m is ln rescaled body mass (relative to mass at maturation)
# Let's do some calculations on that using the ranges in the figure
test <- data.frame(m = seq(-8, 1, 0.1)) # log scale
# Note the intercept is for the mean-centered model. Add the mean back:
# > mean(dat$log_mass_norm_mat)
# [1] -3.324591
test$m_ct <- test$m + 3.324591
test$T_opt <- -0.074 - 0.31*test$m_ct
plot(test$T_opt ~ test$m)
# exponentiate
test$m_linear <- exp(test$m)
plot(test$T_opt ~ test$m_linear)
# Here 1 is mass at maturation
# Assume a fish reaches maturity at 500g
# 5 g is then 5/500 0.01 on the scale above, and 1000 is 2
# Plot again
plot(test$T_opt ~ test$m_linear, xlim = c(0.01, 2))
abline(v = c(0.01, 2), col = "red")
abline(h = c(0.35, -1.32), col = "red")
small_opt <- filter(test, m_linear > 0.009 & m_linear < 0.01)$T_opt
large_opt <- filter(test, m_linear > 1.9 & m_linear < 2.1)$T_opt

small_opt - large_opt

max(test$m_linear)

