#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.06: Max Lindmark
#
# Code to fit and evaluate non-linear (polynomial) models to consumption rates
# 
# A. Load libraries
#
# B. Read data
#
# C. Define models
#
# D. Model selection
#
# E. Model validation
# 
# F. Model fit
# 
# G. Plot predictions
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

# Temp-range used for prediction
temp_pred = seq(from = min(con$temp_norm_ct), 
                to = max(con$temp_norm_ct),
                length.out = 50)

# Plotting is the easiest way to see which masses on normal scale correspond to which 
# standardized and mean centered masses. For prediction, I will use -2, 0 and 2
plot(con$mass_norm ~ con$log_mass_norm_ct)


# C. DEFINE MODEL =================================================================
#** Without Mass-Temperature Interaction ================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred[k] <- b0 + b1*0 + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k]
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.00001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  b0 ~ dnorm(0, 5)
  b1 ~ dnorm(0, 5)
  b2 ~ dnorm(0, 5)
  b3 ~ dnorm(0, 5)
  b4 ~ dnorm(0, 5)
  sigma ~ dunif(0, 10) 
  tau <- 1/sigma^2
  
  }", fill = TRUE, file = "R/analysis/non-linear model/non-linear model.txt")

nl_model = "R/analysis/non-linear model/non-linear model.txt"


# D. FIT MODELS ====================================================================
t <- data.frame()
tt <- data.frame()
df <- data.frame()
jdat <- list()
pred_dat <- list()

for(i in unique(con$species_ab)){

  df <- filter(con, species_ab == i)
  
  jdat = list(
    #y = log(df$y_norm), 
    y = df$y_norm,
    #y = log(df$y),
    #y = df$y,
    #y = df$y - mean(df$y), 
    n_obs = length(df$y), 
    #mass = df$log_mass_norm_ct,
    mass = df$log_mass_ct,
    temp = df$temp_norm_ct,
    #temp = df$temp_norm,
    #temp = df$temp_c,
    #temp_pred = temp_pred,
    #temp_pred = seq(min(con$temp_c), max(con$temp_c), 1),
    temp_pred = seq(min(df$temp_norm_ct), max(df$temp_norm_ct), 0.1))
  
  jm = jags.model(nl_model,
                  data = jdat, 
                  n.adapt = 5000,
                  n.chains = 3) 
  
  burn.in = 10000
  
  update(jm, n.iter = burn.in) 
  
  # Plot model validation
  cs <- coda.samples(jm,
                     variable.names = c("b0", "b1", "b2", "b3",
                                        "mean_y", "mean_y_sim", "p_mean"), 
                     n.iter = 10000, 
                     thin = 5)
  
  p1 <- ggs(cs) %>% 
    filter(Parameter %in% c("b0", "b1", "b2", "b3", "b4")) %>% 
    ggs_density(.) + 
    facet_wrap(~ Parameter, ncol = 1, scales = "free") +
    theme_classic(base_size = 11) + 
    geom_density(alpha = 0.05) +
    scale_color_brewer(palette = "Dark2") + 
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Value", y = "Density") +
    guides(fill = FALSE, color = FALSE) +
    ggtitle(i) +
    NULL
  
  p2 <- ggs(cs) %>% 
    filter(Parameter %in% c("b0", "b1", "b2", "b3", "b4")) %>% 
    ggs_traceplot(.) +
    facet_wrap(~ Parameter, ncol = 1, scales = "free") +
    theme_classic(base_size = 11) + 
    geom_line(alpha = 0.3) +
    scale_color_brewer(palette = "Dark2") + 
    labs(x = "Iteration", y = "Value", color = "Chain #") +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(axis.text.x = element_text(size = 6)) +
    NULL
  
  p3 <- p1 + p2
  
  name <- paste("figures/supp/non-linear/", i, sep = "", "_validation.pdf")
  ggsave(name, plot = p3, scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
  
  # Plot fit (mean)
  cs_fit_df <- data.frame(as.matrix(cs))
  n_bins <- round(1 + 3.2*log(nrow(cs_fit_df)))
  
  p4 <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
    geom_histogram(bins = n_bins) +
    geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
               linetype = 2, size = 0.4) +
    theme_classic(base_size = 11) +
    annotate("text", -Inf, Inf, label = "A", size = 4, 
             fontface = "bold", hjust = -0.5, vjust = 1.3) +
    annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
             label = paste("P =", round(mean(cs_fit_df$p_mean), digits = 3))) +
    labs(x = "Mean simulated consumption", y = "count") +
    theme(aspect.ratio = 1) +
    coord_cartesian(expand = 0) + 
    ggtitle(i) +
    NULL
  name <- paste("figures/supp/non-linear/", i, sep = "", "_fit.pdf")
  ggsave(name, plot = p4, scale = 1, width = 16, height = 16, units = "cm", dpi = 300)
  
  
  js = jags.samples(jm,
                    variable.names = "pred", 
                    n.iter = 10000, 
                    thin = 5)
  
  t <- summary(js$pred, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat
  
  # Create data frame for the predictions
  pred_dat[[i]] <- data.frame(lwr_95 = t[1, ],
                              lwr_80 = t[2, ],
                              median = t[3, ],
                              upr_80 = t[4, ],
                              upr_95 = t[5, ],
                              mass = 0,
                              #temp = temp_pred,
                              #temp = seq(min(df$temp_c), max(df$temp_c), 0.1),
                              #temp = seq(min(df$temp_norm), max(df$temp_norm), 0.1),
                              temp = seq(min(df$temp_norm_ct), max(df$temp_norm_ct), 0.1),
                              species_ab = i)  
}

# Create list from data-frame
pred_dat_df <- dplyr::bind_rows(pred_dat)


# E. PLOT MODELS ===================================================================
# ** Plot predictions ==============================================================

# First, expand the Dark2 palette in RColorBrewer to accomodate 9 levels
colourCount = length(unique(con$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)


# Plot "raw" predictions
pred_dat_df %>% 
  ggplot(., aes(temp, median, color = factor(species_ab))) +
  geom_line() +
  scale_color_manual(values = pal) +
  theme_classic(base_size = 14) + 
  labs(x = "Rescaled temperature",
       y = "Predicted consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4) +
  NULL

# Now do facet by species to check fit individually...
pred_dat_df %>% 
  ggplot(., aes(temp, median)) +
  facet_wrap(~ factor(species_ab), scales = "free") +
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
#ggsave("figures/supp/non-linear/nl_model_species_fit.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# ** Plot standardized predictions, similar to Englund et al 2011 ==================
pred_dat_df <- dplyr::bind_rows(pred_dat)

pred_dat_df <- pred_dat_df %>% 
  dplyr::group_by(species_ab) %>% 
  dplyr::mutate(median_stand = median / max(median))
  
# Find T_opt within species
T_opt_s <- pred_dat_df %>% 
  dplyr::group_by(species_ab) %>% 
  dplyr::filter(median_stand == 1) %>% 
  dplyr::mutate(opt_temp = temp)

# Now get overall mean T_opt
T_opt <- pred_dat_df %>% 
  dplyr::filter(median_stand == 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(mean_temp = mean(temp),
                   stdev_temp = sd(temp))

# Summarize max rate within species
sum <- pred_dat_df %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(species_ab) %>% 
  dplyr::summarise(max_y_pred = max(median)) %>% 
  dplyr::ungroup()

# Filter raw data to plot with predicted
dat <- con %>% 
  select(species_ab, y_norm, temp_norm) %>% 
  rename("temp_dat" = "temp_norm")

# Add max rate by species into long data frame
dat <- left_join(dat, sum)

# Add T_opt by species into long data frame
dat <- left_join(dat, 
                 T_opt_s %>% select(opt_temp, species_ab))

# Standardize data with respect to max predicted and temp with respect to opt within and grand opt
dat <- dat %>% 
  mutate(y_stand_pred = y_norm / max_y_pred,
         y_stand_dat = y_norm / max(y_norm),
         t_stand = temp_dat - opt_temp + T_opt$mean_temp)

# Now we need to standardize predictions as well
# Add in species T_opt in prediction-data
pred_dat_df <- left_join(pred_dat_df,
                         T_opt_s %>% select(opt_temp, species_ab))

# Add in max prediction by species
pred_dat_df <- left_join(pred_dat_df,
                         sum)

# Mutate standardized columns
pred_dat_df <- pred_dat_df %>% 
  dplyr::mutate(pred_stand = median/max(median),
                lwr_95_stand = lwr_95/max(median),
                upr_95_stand = upr_95/max(median),
                lwr_80_stand = lwr_80/max(median),
                upr_80_stand = upr_80/max(median),
                t_stand = temp - opt_temp + T_opt$mean_temp)

# Plot standardized data as in Englund, prediction and credible intervals
pred_dat_df %>% 
  ggplot(., aes(t_stand, median_stand, color = factor(species_ab))) +
  geom_ribbon(data = pred_dat_df, inherit.aes = FALSE,
              aes(x = t_stand, fill = factor(species_ab), ymax = upr_80_stand, ymin = lwr_80_stand),
              size = 1, alpha = 0.15) +
  geom_vline(xintercept = T_opt$mean_temp, linetype = 2, color = "gray20", size = 0.7) + 
  geom_vline(xintercept = T_opt$mean_temp - T_opt$stdev_temp, 
             linetype = 3, color = "gray20", size = 0.7) + 
  geom_vline(xintercept = T_opt$mean_temp + T_opt$stdev_temp, 
             linetype = 3, color = "gray20", size = 0.7) + 
  geom_point(data = dat,  aes(t_stand, y_stand_pred, color = factor(species_ab)),
             size = 3, alpha = 0.6, shape = 21) +
  geom_line(size = 1, alpha = 0.8) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  guides(fill = FALSE) +
  theme_classic(base_size = 14) + 
  labs(x = "Rescaled temperature",
       y = "Rescaled consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4,
        legend.position = c(0.13, 0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12)) +
  NULL
#ggsave("figures/supp/non-linear/nl_model_stand.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)


# Plot non-standardized data and prediction (as it's already divided by the mean)
pred_dat_df %>% 
  ggplot(., aes(temp, median, color = factor(species_ab))) +
  geom_ribbon(data = pred_dat_df, inherit.aes = FALSE,
              aes(x = temp, fill = factor(species_ab), ymax = upr_80, ymin = lwr_80), 
              size = 1, alpha = 0.15) +
  geom_vline(xintercept = T_opt$mean_temp, linetype = 2, color = "gray20", size = 0.7) + 
  geom_vline(xintercept = T_opt$mean_temp - T_opt$stdev_temp, 
             linetype = 3, color = "gray20", size = 0.7) + 
  geom_vline(xintercept = T_opt$mean_temp + T_opt$stdev_temp, 
             linetype = 3, color = "gray20", size = 0.7) + 
  geom_point(data = dat,  aes(temp_dat, y_norm, color = factor(species_ab)),
             size = 3, alpha = 0.6, shape = 21) +
  geom_line(size = 1, alpha = 0.8) +
  geom_segment(data = T_opt_s, aes(x = opt_temp, xend = opt_temp, 
                                   y = min(pred_dat_df$lwr_80), yend = -0.4),
               size = 1, arrow = arrow(length = unit(0.35, "cm")), show.legend = FALSE) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim = c(min(pred_dat_df$lwr_80), max(dat$y_norm*1.03)),
                  xlim = c(min(dat$temp_dat*1.03), max(dat$temp_dat*1.03)),
                  expand = 0) +
  scale_fill_manual(values = pal) +
  guides(fill = FALSE)  +
  theme_classic(base_size = 14) + 
  labs(x = "Rescaled temperature",
       y = "Rescaled consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4,
        legend.position = c(0.13, 0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12)) +
  NULL
#ggsave("figures/nl_model.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)

