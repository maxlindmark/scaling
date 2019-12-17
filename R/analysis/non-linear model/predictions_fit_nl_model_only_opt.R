#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.06: Max Lindmark
#
# - Code to fit non-linear (polynomial) models to consumption using ONLY optimum data
# 
# A. Load libraries
#
# B. Read data
#
# C. Model fit
# 
# D. Plot predictions
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
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       
# readxl_1.3.1      [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1
# tidyr_1.0.0        dplyr_0.8.3       [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file(s)
con <- read.csv("data/con_analysis.csv")

# There is a lot of variation in rates between species (see exploratory script. 
# Find species with optimum data
opt <- con %>% filter(above_optimum == "Y") %>% droplevels()

opt_spec <- unique(opt$species)

# There is a lot of variation in rates between species (see exploratory script. 
# Here we will fit models where y is relative to max for that species
con <- con %>% 
  dplyr::group_by(species) %>% 
  #dplyr::mutate(y_norm = y/max(y)) %>% 
  dplyr::mutate(y_norm = y/mean(y)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(temp_norm_ct = temp_norm - mean(temp_norm)) %>% 
  dplyr::filter(species %in% opt_spec) %>% 
  dplyr::mutate(species_n = as.numeric(as.factor(species))) %>% 
  dplyr::filter(temp_norm_ct < 12) %>% 
  droplevels()

unique(con$species)

# Temp-range used for prediction
temp_pred_con = seq(from = min(con$temp_norm_ct), 
                    to = max(con$temp_norm_ct),
                    length.out = 50)

# Which masses to use for contrasting predictions?
ggplot(con, aes(log_mass_norm_ct)) + 
  geom_histogram() +
  geom_vline(xintercept = c(-2.5, 2.5), color = "red")

ggplot(con, aes(log_mass_norm_ct, mass_norm)) + 
  geom_point() +
  geom_vline(xintercept = c(-2.5, 2.5), color = "red")

small <- -2.5
medium <- 0
large <- 2.5

# Data list for consumption model
data_con = list(
  # y = con$y_norm,
  # y = con$y,
  y = log(con$y_norm), 
  #y = log(con$y), 
  n_obs = length(con$y_norm), 
  mass = con$log_mass_norm_ct,
  temp = con$temp_norm_ct,
  temp_pred = temp_pred_con,
  species_n = con$species_n,
  small = small,
  medium = medium,
  large = large
)


# C. MODEL FIT =====================================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# First fit models. See "cons_nl_model_selection_only_opt"
model = "R/analysis/non-linear model/models/m1.txt"

jm_con = jags.model(model,
                    data = data_con, 
                    n.adapt = 5000, 
                    n.chains = 3)

# Plot mean of simulated data vs mean of observed data
# First convert your matrix 
cs_fit_con = coda.samples(jm_con,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

# Make sure cs now samples simulated data and the mean of the data
cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df_con)))

# Consumption
p1 <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_con$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df_con$p_mean), digits = 3))) +
  labs(x = "Mean simulated consumption", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  #ggtitle("Consumption") +
  NULL

p2 <- ggplot(cs_fit_df_con, aes(cv_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_con$cv_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "D", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -2, vjust = 3.3, 
           label = paste("P =", round(mean(cs_fit_df_con$p_cv), digits = 3))) +
  labs(x = "cv simulated consumption", y = "") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) +
  NULL

p1
#ggsave("figures/supp/nl_cv_mean_fit_only_opt.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# D. PLOT PREDICTIONS ==============================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

js_con = jags.samples(jm_con,
                      variable.names = c("pred_large", "pred_medium", "pred_small"), 
                      n.iter = samples, 
                      thin = n.thin)

# Generate medians and quantiles that can be used for storing info, plotting etc.
# Large size:
c_pred_large <- summary(js_con$pred_large, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_large_df <- data.frame(lwr_95 = c_pred_large[1, ],
                              lwr_80 = c_pred_large[2, ],
                              median = c_pred_large[3, ],
                              upr_80 = c_pred_large[4, ],
                              upr_95 = c_pred_large[5, ],
                              mass = large,
                              temp = temp_pred_con)

# Medium size:
c_pred_medium <- summary(js_con$pred_medium, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_medium_df <- data.frame(lwr_95 = c_pred_medium[1, ],
                               lwr_80 = c_pred_medium[2, ],
                               median = c_pred_medium[3, ],
                               upr_80 = c_pred_medium[4, ],
                               upr_95 = c_pred_medium[5, ],
                               mass = medium,
                               temp = temp_pred_con)

# Small size:
c_pred_small <- summary(js_con$pred_small, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_small_df <- data.frame(lwr_95 = c_pred_small[1, ],
                              lwr_80 = c_pred_small[2, ],
                              median = c_pred_small[3, ],
                              upr_80 = c_pred_small[4, ],
                              upr_95 = c_pred_small[5, ],
                              mass = small,
                              temp = temp_pred_con)


#** Plot ===========================================================================
# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
pal <- brewer.pal("Dark2", n = 5)

x_min <- min(con$temp_norm_ct, con$temp_norm_ct)
x_max <- max(con$temp_norm_ct, con$temp_norm_ct)

c_pdat <- rbind(c_pred_large_df, c_pred_medium_df, c_pred_small_df)

# Calcualte stuf for the arrows in the below figure
max_rates_con <- c_pdat %>% 
  group_by(factor(mass)) %>% 
  filter(median == max(median)) %>% 
  arrange(mass)

y_end <- min(log(con$y_norm))
x_s_con <- max_rates_con$temp[1]
x_m_con <- max_rates_con$temp[2]
x_l_con <- max_rates_con$temp[3]

p3 <- ggplot(c_pdat, aes(temp, median, color = factor(mass))) +
  geom_point(data = con, aes(temp_norm_ct, log(y_norm)), size = 2.8, shape = 21,
             alpha = 0.2, color = "white", fill = "grey40") +
  # geom_point(data = con, aes(temp_norm_ct, log(y)), size = 2.8, shape = 21, 
  #            alpha = 0.2, color = "white", fill = "grey40") +
  geom_ribbon(data = c_pdat, aes(x = temp, ymin = lwr_95, ymax = upr_95, fill = factor(mass)), 
              size = 0.6, alpha = 0.25, inherit.aes = FALSE) +
  geom_ribbon(data = c_pdat, aes(x = temp, ymin = lwr_80, ymax = upr_80, fill = factor(mass)), 
              size = 0.6, alpha = 0.4, inherit.aes = FALSE) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  geom_line(size = 0.6, alpha = 1) +
  theme_classic(base_size = 14) + 
  labs(x = expression("Centered Temperature " ( degree*C)),
       y = "ln(standardized\nconsumpption rate)",
       color = "ln(mass)\n(centered)") +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  theme(aspect.ratio = 3/4,
        legend.position = "bottom") +
  guides(fill = FALSE) +
  xlim(x_min, x_max) +
  # geom_segment(aes(x = x_s_con, xend = x_s_con, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[1]) +
  # geom_segment(aes(x = x_m_con, xend = x_m_con, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[2]) +
  # geom_segment(aes(x = x_l_con, xend = x_l_con, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[3]) +
  NULL

#**** Together ==================================================================
p3
#ggsave("figures/nl_model_only_opt.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# NOTES:
# LOG Y NORM: OPTIMUM DECLINES WITH SIZE
# Y NORM: OPTIMUM INCREASES WITH SIZE
# Y: NO OPTIMUM! HERE WE DON*T LOOK AT RELATIVE CONSUMPTION AND THIS WILL AFFECT THINGS IF SPECIES HAVE DIFFERENT THERMAL RANGES AND SIZE COMBINATIONS
# LOG Y: ALL HAVE THE SAME OPTIMUM (OR ACTUALLY NO OPTIMUM EVIDENT)
