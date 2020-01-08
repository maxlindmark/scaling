#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.06: Max Lindmark
#
# - Code to fit non-linear (polynomial) models to consumption 
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
  dplyr::ungroup() %>% 
  dplyr::mutate(temp_norm_ct = temp_norm - mean(temp_norm))

# Which species have data above optimum?
spec <- unique(filter(con, above_optimum == "Y"))$species
  
con %>% 
  filter(species %in% spec) %>% 
  #ggplot(., aes(temp_norm_ct, y_norm, color = species)) +
  ggplot(., aes(temp_norm_ct, y, color = species)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  stat_smooth(se = FALSE) +
  labs(x = "Standardized temperature", y = "Standardized consumption") +
  scale_color_viridis(discrete = TRUE, option = "magma") +
  NULL

con %>% 
  filter(species %in% spec) %>% 
  #ggplot(., aes(temp_norm_ct, y_norm, color = log_mass_norm)) +
  ggplot(., aes(temp_norm_ct, y, color = log_mass_norm)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  stat_smooth(se = FALSE) +
  labs(x = "Standardized temperature", y = "Standardized consumption") +
  scale_color_viridis(option = "magma") +
  facet_wrap(~ species, scales = "free") +
  NULL

con <- con %>% filter(species %in% spec)

# Temp-range used for prediction
temp_pred = seq(from = min(con$temp_norm_ct), 
                to = max(con$temp_norm_ct),
                length.out = 50)


# # Data list for JAGS model
# con_data = list(
#   # y = con$y_norm,
#   # y = con$y,
#   y = log(con$y_norm), 
#   #y = log(con$y), 
#   n_obs = length(con$y_norm), 
#   mass = con$log_mass_norm_ct,
#   temp = con$temp_norm_ct,
#   temp_pred = temp_pred,
#   species_n = as.numeric(con$species)
# )


# Plotting is the easiest way to see which masses on normal scale correspond to which 
# standardized and mean centered masses. For prediction, I will use -2, 0 and 2
plot(con$mass_norm ~ con$log_mass_norm_ct)


# C. DEFINE MODELS =================================================================
#** Mass-Temperature Interaction ===================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i] + b4*temp[i]*mass[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred[k] <- b0 + b1*0 + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*0
    
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
  
  }", fill = TRUE, file = "R/analysis/non-linear model/interaction.txt")

model_inter = "R/analysis/non-linear model/interaction.txt"


#** No Mass-Temperature Interaction ================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1*mass[i] + b2*temp[i] + b3*temp[i]*temp[i] + b4*temp[i]*mass[i]
  
    # Add log likelihood computation for each observation
    pd[i] <- dnorm(y[i], mu[i], tau)
  
    # Calculates the log PPD
    log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Predictions
  for(k in 1:length(temp_pred)){
      
    pred[k] <- b0 + b1*0 + b2*temp_pred[k] + b3*temp_pred[k]*temp_pred[k] + b4*temp_pred[k]*0
    
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
  
  }", fill = TRUE, file = "R/analysis/non-linear model/no_interaction.txt")

model_no_inter = "R/analysis/non-linear model/no_interaction.txt"


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
    #y = df$y - mean(df$y), 
    n_obs = length(df$y), 
    mass = df$log_mass_norm_ct,
    #temp = df$temp_norm_ct,
    temp = df$temp_norm,
    #temp = df$temp_c,
    #temp_pred = temp_pred,
    #temp_pred = seq(min(con$temp_c), max(con$temp_c), 1),
    temp_pred = seq(min(con$temp_norm), max(con$temp_norm), 1))
  
  jm = jags.model(model_no_inter,
                  data = jdat, 
                  n.adapt = 5000, # 5000, 
                  n.chains = 3) #3)
  
  burn.in = 10000
  
  update(jm, n.iter = burn.in) 
  
  cs <- coda.samples(jm,
                     variable.names = c("b0", "b1", "b2", "b3", "b4"), 
                     n.iter = 10000, 
                     thin = 5)
  
  p1 <- ggs(cs) %>% 
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
  
  name <- paste("figures/supp/non-linear/", i, sep = "", ".pdf")
  ggsave(name, plot = p3, scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
  
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
                              #temp = seq(min(con$temp_c), max(con$temp_c), 1),
                              temp = seq(min(con$temp_norm), max(con$temp_norm), 1),
                              species = i)  
}

# Create list from data-frame
pred_dat_df <- dplyr::bind_rows(pred_dat)

str(pred_dat_df)

# Standardize predictions, as in Englund et al 2011
pred_dat_df <- pred_dat_df %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(median_stand = median/max(median))

# Plot "raw" predictions
pred_dat_df %>% 
  filter(median > 0) %>% 
  ggplot(., aes(temp, median, color = factor(species))) +
  geom_line() +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 14) + 
  labs(x = "Standardized temperature",
       y = "Standardized consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4) +
  NULL

# Plot standardized predictions
pred_dat_df %>% 
  filter(median > 0) %>% # try with and without this...
  ggplot(., aes(temp, median_stand, color = factor(species))) +
  geom_line() +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 14) + 
  labs(x = "Standardized temperature",
       y = "Standardized consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4) +
  NULL


# Plot standardized predictions with standardized data
sub <- con %>% 
  select(species_ab, y_norm, temp_norm) %>% 
  rename("species" = "species_ab",
         "temp_dat" = "temp_norm")

# First, summarize max rate within species
sum <- pred_dat_df %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(max_y_pred = max(median))

sum

sub2 <- left_join(sub, sum)

# Now standardize y wtr max predicted
sub2 <- sub2 %>% mutate(y_stand = y_norm / max_y_pred)

# Now get overall mean T_opt
sum_tot <- pred_dat_df %>% 
  dplyr::filter(median_stand == 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::summarise(mean_temp = mean(temp),
                   stdev_temp = sd(temp))


# Plot normalized data, prediction and credible intervals
pred_dat_df %>% 
  filter(median > 0) %>% 
  ggplot(., aes(temp, median_stand, color = factor(species))) +
  geom_ribbon(data = filter(pred_dat_df, median > 0), inherit.aes = FALSE,
              aes(x = temp, fill = factor(species), ymax = upr_80, ymin = lwr_80), 
              size = 1, alpha = 0.15) +
  geom_vline(xintercept = sum_tot$mean_temp, linetype = 2, color = "gray20", size = 0.7) + 
  geom_point(data = filter(sub2, y_stand > 0),  aes(temp_dat, y_stand, color = factor(species)),
             size = 3, alpha = 0.8, shape = 21) +
  geom_line(size = 1, alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = FALSE) +
  theme_classic(base_size = 14) + 
  labs(x = "Standardized temperature",
       y = "Standardized consumption rate",
       color = "Species") +
  coord_cartesian(ylim = c(0, 6)) +
  theme(aspect.ratio = 3/4,
        legend.position = c(0.13, 0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12)) +
  NULL
#ggsave("figures/nl_model.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)

# Plot normalized data, prediction : NO credible intervals
pred_dat_df %>% 
  filter(median > 0) %>% 
  ggplot(., aes(temp, median_stand, color = factor(species))) +
  geom_vline(xintercept = sum_tot$mean_temp, linetype = 2, color = "gray20", size = 0.7) + 
  geom_point(data = filter(sub2, y_stand > 0),  aes(temp_dat, y_stand, color = factor(species)),
             size = 3, alpha = 0.8, shape = 21) +
  geom_line(size = 1, alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = FALSE) +
  theme_classic(base_size = 14) + 
  labs(x = "Standardized temperature",
       y = "Standardized consumption rate",
       color = "Species") +
  theme(aspect.ratio = 3/4,
        legend.position = c(0.13, 0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12)) +
  NULL
#ggsave("figures/nl_model.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)



# Think it could be correct and that the "problem" is that I simply predict for a given 
# body mass of 0


##------ CONTINUE HERE.... DID I REALLY STANDARDIZE PREDICTIONS CORRECTLY? check left 
# join with miniale example


# Then do model fit and within loop
# then repeat for T_opt model. Look at the posterior of the interaction term.











# pred_dat_df %>% filter(median > 0) %>% ggplot(., aes(temp, median, color = factor(species))) +
#   facet_wrap(~ species) +
#   scale_color_brewer(palette = "Set1") +
#   geom_line() +
#   NULL

  # geom_point(data = met, aes(temp_norm_ct, log(y_norm)), size = 2.8, shape = 21,
  #            alpha = 0.2, color = "white", fill = "grey40") +
  # geom_ribbon(data = m_pdat, aes(x = temp, ymin = lwr_95, ymax = upr_95, fill = factor(mass)), 
  #             size = 0.6, alpha = 0.25, inherit.aes = FALSE) +
  # geom_ribbon(data = m_pdat, aes(x = temp, ymin = lwr_80, ymax = upr_80, fill = factor(mass)), 
  #             size = 0.6, alpha = 0.4, inherit.aes = FALSE) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  geom_line(size = 0.6, alpha = 1) +
  theme_classic(base_size = 14) + 
  labs(x = "",
       y = "ln(standardized\nmetabolic rate)") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  theme(aspect.ratio = 3/4) +
  guides(fill = FALSE, color = FALSE) +
  xlim(x_min, x_max) +
  geom_segment(aes(x = x_s_met, xend = x_s_met, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[1]) +
  geom_segment(aes(x = x_m_met, xend = x_m_met, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[2]) +
  geom_segment(aes(x = x_l_met, xend = x_l_met, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[3]) +
  NULL





# D. MODEL SELECTION ===============================================================
#** Metabolism =====================================================================
#**** With interaction =============================================================
jm_met_inter = jags.model(model_inter,
                          data = met_data, 
                          n.adapt = 5000, 
                          n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm_met_inter, n.iter = burn.in) 

#dic.samples(jm_con_inter, n.iter = 2500, type = "pD")

# Monitor the likelihood to calculate WAIC
zj_met_inter = jags.samples(jm_met_inter, 
                            variable.names = c("pd", "log_pd"), 
                            n.iter = 10000, 
                            thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd_met_inter <- -2*sum(log(summary(zj_met_inter$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC_met_inter <- sum((summary(zj_met_inter$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC_met_inter <- lppd_met_inter + 2*pd.WAIC_met_inter

# > WAIC_met_inter
# [1] 7465.851

# MEAN-standardizing instead of MAX
# > WAIC_met_inter
# [1] 6850.653

#**** Without interaction =============================================================
jm_met = jags.model(model,
                    data = met_data, 
                    n.adapt = 5000, 
                    n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm_met, n.iter = burn.in) 

#dic.samples(jm_con, n.iter = 2500, type = "pD")

# Monitor the likelihood to calculate WAIC
zj_met = jags.samples(jm_met, 
                      variable.names = c("pd", "log_pd"), 
                      n.iter = 10000, 
                      thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd_met <- -2*sum(log(summary(zj_met$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC_met <- sum((summary(zj_met$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC_met <- lppd_met + 2*pd.WAIC_met

# > WAIC_met
# [1] 7513.53

# MEAN-standardizing instead of MAX
# > WAIC_met
# [1] 6882.178


#** Consumption ====================================================================
#**** With interaction =============================================================
jm_con_inter = jags.model(model_inter,
                          data = con_data, 
                          n.adapt = 5000, 
                          n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm_con_inter, n.iter = burn.in) 

#dic.samples(jm_con_inter, n.iter = 2500, type = "pD")

# Monitor the likelihood to calculate WAIC
zj_con_inter = jags.samples(jm_con_inter, 
                            variable.names = c("pd", "log_pd"), 
                            n.iter = 10000, 
                            thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd_con_inter <- -2*sum(log(summary(zj_con_inter$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC_con_inter <- sum((summary(zj_con_inter$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC_con_inter <- lppd_con_inter + 2*pd.WAIC_con_inter
# > WAIC_con_inter
# [1] 1465.439

# MEAN-standardizing instead of MAX
# > WAIC_con_inter
# [1] 1399.706


#**** Without interaction =============================================================
jm_con = jags.model(model,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm_con, n.iter = burn.in) 

#dic.samples(jm_con, n.iter = 2500, type = "pD")

# Monitor the likelihood to calculate WAIC
zj_con = jags.samples(jm_con, 
                      variable.names = c("pd", "log_pd"), 
                      n.iter = 10000, 
                      thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd_con <- -2*sum(log(summary(zj_con$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC_con <- sum((summary(zj_con$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC_con <- lppd_con + 2*pd.WAIC_con
# > WAIC_con
# [1] 1466.829

# MEAN-standardizing instead of MAX
# > WAIC_con
# [1] 1399.519


# E. MODEL VALIDATION ==============================================================
#** Metabolism =====================================================================
# Sample from the posterior ========================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# CODA - Nice for getting the raw posteriors
cs_met <- coda.samples(#jm_met,
                       jm_met_inter,
                       variable.names = c("b0", "b1", "b2", "b3", "sigma"), 
                       n.iter = samples, 
                       thin = n.thin)

summary(cs_met)

# Evaluate convergence =============================================================
# Convert to ggplottable data frame
cs_met_df <- ggs(cs_met)

# Plot posterior densities of parameters and predictions
p3 <- cs_met_df %>% 
  filter(Parameter %in% c("b0", "b1", "b2", "b3", "sigma")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 1, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density") +
  guides(fill = FALSE, color = FALSE) +
  ggtitle("Metabolism") +
  NULL

# Traceplot for evaluating chain convergence
p4 <- cs_met_df %>% 
  filter(Parameter %in% c("b0", "b1", "b2", "b3", "b4", "sigma")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 1, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL

p3 + p4
#ggsave("figures/supp/nl_model_validation_met.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


#** Consumption ====================================================================
# Sample from the posterior ========================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# CODA - Nice for getting the raw posteriors
cs_con <- coda.samples(jm_con_inter,
                       variable.names = c("b0", "b1", "b2", "b3", "b4", "sigma", "mu_b0"), 
                       n.iter = samples, 
                       thin = n.thin)

summary(cs_con)

# Evaluate convergence =============================================================
# Convert to ggplottable data frame
cs_con_df <- ggs(cs_con)

# Plot posterior densities of parameters and predictions
p5 <- cs_con_df %>% 
  filter(Parameter %in% c("b0", "b1", "b2", "b3", "b4", "sigma")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 1, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density") +
  guides(fill = FALSE, color = FALSE) +
  ggtitle("Maximum consumption") +
  NULL

# Traceplot for evaluating chain convergence
p6 <- cs_con_df %>% 
  filter(Parameter %in% c("b0", "b1", "b2", "b3", "b4", "sigma")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 1, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL

p5 + p6
#ggsave("figures/supp/nl_model_validation_con.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


# F. MODEL FIT =====================================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# Plot mean of simulated data vs mean of observed data
# First convert your matrix 
cs_fit_met = coda.samples(jm_met_inter,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

cs_fit_con = coda.samples(jm_con_inter,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

# Make sure cs now samples simulated data and the mean of the data
cs_fit_df_met <- data.frame(as.matrix(cs_fit_met))
cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df_met)))

# Metabolism
p7 <- ggplot(cs_fit_df_met, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_met$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df_met$p_mean), digits = 3))) +
  labs(x = "Mean simulated metabolism", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  #ggtitle("Metabolism") +
  NULL

p8 <- ggplot(cs_fit_df_met, aes(cv_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_met$cv_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3, 
           label = paste("P =", round(mean(cs_fit_df_met$p_cv), digits = 3))) +
  labs(x = "cv simulated metabolism", y = "") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) +
  NULL

# Consumption
p9 <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
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

p10 <- ggplot(cs_fit_df_con, aes(cv_y_sim)) + 
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

(p7 + p8) / (p9 + p10)

p7 + p9
#ggsave("figures/supp/nl_cv_mean_fit.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# G. PLOT PREDICTIONS ==============================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

#** Metabolism =====================================================================
js_met = jags.samples(jm_met_inter,
                      #jm_met, 
                      variable.names = c("pred_large", "pred_medium", "pred_small"), 
                      n.iter = samples, 
                      thin = n.thin)

# Generate medians and quantiles that can be used for storing info, plotting etc.
# Large size:
m_pred_large <- summary(js_met$pred_large, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
m_pred_large_df <- data.frame(lwr_95 = m_pred_large[1, ],
                              lwr_80 = m_pred_large[2, ],
                              median = m_pred_large[3, ],
                              upr_80 = m_pred_large[4, ],
                              upr_95 = m_pred_large[5, ],
                              mass = 4,
                              temp = temp_pred_met)

# Medium size:
m_pred_medium <- summary(js_met$pred_medium, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
m_pred_medium_df <- data.frame(lwr_95 = m_pred_medium[1, ],
                               lwr_80 = m_pred_medium[2, ],
                               median = m_pred_medium[3, ],
                               upr_80 = m_pred_medium[4, ],
                               upr_95 = m_pred_medium[5, ],
                               mass = 0,
                               temp = temp_pred_met)

# Small size:
m_pred_small <- summary(js_met$pred_small, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
m_pred_small_df <- data.frame(lwr_95 = m_pred_small[1, ],
                              lwr_80 = m_pred_small[2, ],
                              median = m_pred_small[3, ],
                              upr_80 = m_pred_small[4, ],
                              upr_95 = m_pred_small[5, ],
                              mass = -4,
                              temp = temp_pred_met)

plot(m_pred_large_df$median ~ m_pred_large_df$temp)


#** Consumption ====================================================================
js_con = jags.samples(#jm_con_inter,
                      jm_con, 
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
                              mass = 4,
                              temp = temp_pred_con)

# Medium size:
c_pred_medium <- summary(js_con$pred_medium, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_medium_df <- data.frame(lwr_95 = c_pred_medium[1, ],
                               lwr_80 = c_pred_medium[2, ],
                               median = c_pred_medium[3, ],
                               upr_80 = c_pred_medium[4, ],
                               upr_95 = c_pred_medium[5, ],
                               mass = 0,
                               temp = temp_pred_con)

# Small size:
c_pred_small <- summary(js_con$pred_small, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_small_df <- data.frame(lwr_95 = c_pred_small[1, ],
                              lwr_80 = c_pred_small[2, ],
                              median = c_pred_small[3, ],
                              upr_80 = c_pred_small[4, ],
                              upr_95 = c_pred_small[5, ],
                              mass = -4,
                              temp = temp_pred_con)

plot(c_pred_large_df$median ~ c_pred_large_df$temp)


#** Plot ===========================================================================
# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
pal <- brewer.pal("Dark2", n = 5)

x_min <- min(met$temp_norm_ct, con$temp_norm_ct)
x_max <- max(met$temp_norm_ct, con$temp_norm_ct)

#**** Metabolism ===================================================================
m_pdat <- rbind(m_pred_large_df, m_pred_medium_df, m_pred_small_df)

# Calcualte stuf for the arrows in the below figure
max_rates_met <- m_pdat %>% 
  group_by(factor(mass)) %>% 
  filter(median == max(median)) %>% 
  arrange(mass)

y_end <- min(log(met$y_norm))
x_s_met <- max_rates_met$temp[1]
x_m_met <- max_rates_met$temp[2]
x_l_met <- max_rates_met$temp[3]

p11 <- ggplot(m_pdat, aes(temp, median, color = factor(mass))) +
  geom_point(data = met, aes(temp_norm_ct, log(y_norm)), size = 2.8, shape = 21,
             alpha = 0.2, color = "white", fill = "grey40") +
  # geom_point(data = met, aes(temp_norm_ct, log(y)), size = 2.8, shape = 21, 
  #            alpha = 0.2, color = "white", fill = "grey40") +
  geom_ribbon(data = m_pdat, aes(x = temp, ymin = lwr_95, ymax = upr_95, fill = factor(mass)), 
              size = 0.6, alpha = 0.25, inherit.aes = FALSE) +
  geom_ribbon(data = m_pdat, aes(x = temp, ymin = lwr_80, ymax = upr_80, fill = factor(mass)), 
              size = 0.6, alpha = 0.4, inherit.aes = FALSE) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  geom_line(size = 0.6, alpha = 1) +
  theme_classic(base_size = 14) + 
  labs(x = "",
       y = "ln(standardized\nmetabolic rate)") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  theme(aspect.ratio = 3/4) +
  guides(fill = FALSE, color = FALSE) +
  xlim(x_min, x_max) +
  geom_segment(aes(x = x_s_met, xend = x_s_met, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[1]) +
  geom_segment(aes(x = x_m_met, xend = x_m_met, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[2]) +
  geom_segment(aes(x = x_l_met, xend = x_l_met, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[3]) +
  NULL


#**** Consumption ==================================================================
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

p12 <- ggplot(c_pdat, aes(temp, median, color = factor(mass))) +
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
  geom_segment(aes(x = x_s_con, xend = x_s_con, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[1]) +
  geom_segment(aes(x = x_m_con, xend = x_m_con, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[2]) +
  geom_segment(aes(x = x_l_con, xend = x_l_con, y = -5, yend = y_end), arrow = arrow(length = unit(0.3, "cm")), col = pal[3]) +
  NULL

#**** Together ==================================================================
p11 / p12
#ggsave("figures/nl_model.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# NOTES:
# LOG Y NORM: OPTIMUM DECLINES WITH SIZE
# Y NORM: OPTIMUM INCREASES WITH SIZE
# Y: NO OPTIMUM! HERE WE DON*T LOOK AT RELATIVE CONSUMPTION AND THIS WILL AFFECT THINGS IF SPECIES HAVE DIFFERENT THERMAL RANGES AND SIZE COMBINATIONS
# LOG Y: ALL HAVE THE SAME OPTIMUM (OR ACTUALLY NO OPTIMUM EVIDENT)


test <- filter(con, temp_norm_ct > 12) 
data.frame(test)

test <- filter(met, temp_norm_ct > 16) 
data.frame(test)


