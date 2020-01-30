#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.06: Max Lindmark
#
# Code to fit quadratic model for a species with good data to exemplify difference
# between Cmax (this model) and the general Arrhenius equation
# 
# A. Load libraries
#
# B. Read Cmax data
#
# C. Define quadratic model for Cmax 
#
# D. Simulate data for plotting, 
#
# E. Plot
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
# Instead of fitting an hierarchial model here (for now at least), we will fit models
# of relative rates, i.e. relative to max for that species

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(con$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", con$species )

con$species_ab <- paste(sp1, sp2, sep = ".")

con <- con %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(y_norm = y/mean(y)) %>% 
  dplyr::mutate(log_mass_ct = log(mass_g) - mean(log(mass_g))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(temp_norm_ct = temp_norm - mean(temp_norm))

ggplot(con, aes(log_mass_ct)) + 
  geom_histogram() + 
  facet_wrap(~species)

# Which species have data above optimum?
spec <- unique(filter(con, above_optimum == "Y"))$species

con %>% 
  filter(species %in% spec) %>% 
  #ggplot(., aes(temp_norm_ct, y_norm, color = species)) +
  ggplot(., aes(temp_norm_ct, y, color = species)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  stat_smooth(se = FALSE) +
  labs(x = "Rescaled temperature", y = "Rescaled consumption") +
  scale_color_viridis(discrete = TRUE, option = "magma") +
  NULL

con %>% 
  filter(species %in% spec) %>% 
  #ggplot(., aes(temp_norm_ct, y_norm, color = log_mass_norm)) +
  ggplot(., aes(temp_norm_ct, y, color = log_mass_norm)) + 
  #ggplot(., aes(temp_norm, y_norm, color = log_mass_norm)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  stat_smooth(se = FALSE) +
  labs(x = "Rescaled temperature", y = "Rescaled consumption") +
  scale_color_viridis(option = "magma") +
  facet_wrap(~ species, scales = "free") +
  NULL

con <- con %>% filter(species %in% spec)


# C. DEFINE MODEL =================================================================
# See non_linear_model.R

nl_model = "R/analysis/non-linear model/non-linear model.txt"

df <- filter(con, species_ab == "C.argus")
  
  jdat = list(
    y = df$y_norm,
    n_obs = length(df$y), 
    mass = df$log_mass_ct,
    temp = df$temp_norm_ct,
    temp_pred = seq(min(df$temp_norm_ct), max(df$temp_norm_ct), 0.1))
  
  jm = jags.model(nl_model,
                  data = jdat, 
                  n.adapt = 5000,
                  n.chains = 3) 
  
  burn.in = 10000
  
  update(jm, n.iter = burn.in) 
  
  # Create data frame for the predictions
  pred_dat <- data.frame(lwr_95 = t[1, ],
                         lwr_80 = t[2, ],
                         median = t[3, ],
                         upr_80 = t[4, ],
                         upr_95 = t[5, ],
                         mass = 0,
                         temp = seq(min(df$temp_norm_ct), max(df$temp_norm_ct), 0.1))  

# Create list from data-frame
pred_dat_df <- dplyr::bind_rows(pred_dat)


# Plot "raw" predictions
pred_dat_df %>% 
  ggplot(., aes(temp, median)) +
  geom_point(data = con, aes(temp_norm_ct, y_norm),
             size = 3, alpha = 0.8, shape = 21, fill = "grey20", color = "white") +
  geom_line(size = 1, alpha = 0.6) +
  geom_ribbon(data = filter(pred_dat_df, median > 0), 
              aes(x = temp, ymax = upr_80, ymin = lwr_80), 
              size = 1, alpha = 0.25, fill = "red") +
  theme_classic(base_size = 14) + 
  guides(fill = FALSE) +
  labs(x = "Rescaled temperature",
       y = "Rescaled consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4,
        legend.position = c(0.13, 0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12)) +
  NULL

# Now check parameter estimates:

summary(cs)
              # Mean  SD        Naive SE       Time-series SE
# b0          0.82444 0.0773428 9.985e-04      1.050e-03
# b1          0.45867 0.0841213 1.086e-03      1.086e-03
# b2          0.11493 0.0134257 1.733e-04      2.768e-04
# b3         -0.00502 0.0009126 1.178e-05      1.902e-05


# C. SIMULATE DATA FOR PLITTING ====================================================
summary(df$log_mass_ct)
summary(df$temp_norm_ct)

mass <- c(0, 0.5)
temp <- seq(-10, 20, 0.1)

# What is 0.5 here?
#filter(df, log_mass_ct, log_mass_ct > 0.4 & log_mass_ct < 0.6)
ggplot(df, aes(log_mass_ct, mass_g)) + 
  geom_line() +
  coord_cartesian(xlim = c(0, 0.5)) +
  NULL

dat <- expand.grid(mass = mass,
                   temp = temp)

dat$temp_arr <- 1/((dat$temp + 273.15) * 8.617332e-05)

dat$cmax <- 0.82444 + 0.45867*dat$mass + 0.11493*dat$temp -0.00502*dat$temp*dat$temp

dat$log_met <- 1.59 + 0.77*dat$mass -0.61*dat$temp_arr + 0.017*dat$temp_arr*dat$mass

dat$met <- exp(dat$log_met)

dat$met_stand <- dat$met/max(dat$met)
dat$cmax_stand <- dat$cmax/max(dat$cmax)

dat$diff <- dat$cmax_stand - dat$met_stand 

# C. PLOT ==========================================================================
dat2 <- dat %>% 
  select(mass, temp, met_stand, cmax_stand, diff) %>% 
  pivot_longer(cols = 3:5)

T_opt_s <- dat2 %>% 
  group_by(mass) %>% 
  filter(name == "diff") %>% 
  filter(value == max(value))

T_opt_s


dat2$name <- factor(dat2$name, levels = c("cmax_stand", "met_stand", "diff"))

pal <- RColorBrewer::brewer.pal("Dark2", n = 3)

ggplot(
  dat2, aes(temp, value, color = factor(name), linetype = factor(mass), alpha = factor(name))
  #dat2, aes(temp, value, color = factor(mass), linetype = factor(name), alpha = factor(mass))
  ) +
  # facet_wrap(~ name) +
  #facet_wrap(~ mass) +
  geom_line(size = 1) +
  coord_cartesian(ylim = c(0, 1.1),
                  xlim = c(-5, 20), 
                  expand = 0) + 
  theme_classic(base_size = 14) + 
  # scale_color_brewer(palette = "Dark2",
  #                    labels = c("Consumption", "Metabolism", "Consumption-Meteabolism")
  #                    ) +
  scale_color_manual(values = c(pal[1], pal[3], pal[2]),
                     labels = c("Consumption", "Metabolism", "Consumption-Metabolism")) +
  scale_linetype_manual(values = c(1, 2), 
                        labels = c("Small", "Large")) +
  scale_alpha_manual(values = c(0.75, 0.75, 1)) +
  labs(x = expression(paste("Rescaled temperature [", degree*C, "]")),
       y = "Rescaled rates",
       color = "Rate",
       linetype = "Rescaled mass"
       # color = "Mass",
       # linetype = "Rate"
       ) +
  geom_segment(data = T_opt_s, aes(x = temp, xend = temp, y = c(T_opt_s$value[1], T_opt_s$value[2]), yend = 0),
               size = 1, arrow = arrow(length = unit(0.35, "cm")), show.legend = FALSE) +
  guides(alpha = FALSE,
         linetype = guide_legend(override.aes = list(size = 0.5))) +
  theme(aspect.ratio = 3/4,
        #legend.text = element_text(size = 8),
        legend.title = element_text(size = 12)) +
  NULL

#ggsave("figures/concept.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


