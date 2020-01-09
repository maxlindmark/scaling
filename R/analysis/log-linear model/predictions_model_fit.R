#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of maximum consumption rate as a function of 
# temperature with different group-effects and compare DIC and WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit models
# 
# E. Evaluate model fit
#
# D. Plot predicted mass-scaling slopes in different temperatures
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
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file(s)
met <- read.csv("data/met_analysis.csv")
con <- read.csv("data/con_analysis.csv")

# Filter data points at below optimum temperatures
met <- met %>% filter(above_optimum == "N")
con <- con %>% filter(above_optimum == "N")

# Data in list-format for JAGS
# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species))
con$species_n <- as.numeric(as.factor(con$species))

# Mass-range used for prediction
mass_pred_con = seq(from = min(con$log_mass_norm_ct), 
                    to = max(con$log_mass_norm_ct),
                    length.out = 100)

# Data list for consumption model
con_data = list(
  y = log(con$y), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_norm_ct,
  temp = con$temp_norm_arr_ct,
  mass_pred_con = mass_pred_con
)

# Mass-range used for prediction
mass_pred_met = seq(from = min(met$log_mass_norm_ct), 
                    to = max(met$log_mass_norm_ct),
                    length.out = 100)

# Data list for metabolism model
met_data = list(
  y = log(met$y), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_norm_ct,
  temp = met$temp_norm_arr_ct,
  mass_pred_met = mass_pred_met
)

# Check which temperatures to use
ggplot(met, aes(temp_norm_arr_ct, temp_norm)) + 
  geom_point()


# C. FIT MODELS ====================================================================
# Refit chosen models from the model selection part
#**** Metabolism (M1) ==============================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +                # varying intercept 
      b1[species_n[i]]*mass[i] +        # varying mass-exponent
      b2[species_n[i]]*temp[i] +        # varying activation energy
      b3[species_n[i]]*mass[i]*temp[i]  # varying M*T interaction
  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
    b3[j] ~ dnorm(mu_b3, tau_b3)
  }
  
  # Predictions
  for(k in 1:length(mass_pred_met)){
      
    pred_warm[k] <- mu_b0 + mu_b1*mass_pred_met[k] + mu_b2*-1 + mu_b3*mass_pred_met[k]*-1
    pred_cold[k] <- mu_b0 + mu_b1*mass_pred_met[k] + mu_b2*1 + mu_b3*mass_pred_met[k]*1
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # varying mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying activation energy
  mu_b3 ~ dnorm(0, 0.5)      # varying interaction
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  sigma_b3 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  tau_b3 <- 1/sigma_b3^2
  
  }", fill = TRUE, file = "R/analysis/log-linear model/models/m1_metabolism_pred.txt")

met_model = "R/analysis/log-linear model/models/m1_metabolism_pred.txt"

jm_met = jags.model(met_model,
                    data = met_data, 
                    n.adapt = 5000, 
                    n.chains = 3)


#**** Consumption (M5) =============================================================
# Predictions
cat(
  "model{
  
    for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +           # varying intercept 
      b1[species_n[i]]*mass[i] +   # varying mass-exponent
      b2[species_n[i]]*temp[i]     # varying activation energy
  
  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
    b2[j] ~ dnorm(mu_b2, tau_b2)
  }

    for(k in 1:length(mass_pred_con)){
  
        pred_warm[k] <- mu_b0 + mu_b1*mass_pred_con[k] + mu_b2*-1
        pred_cold[k] <- mu_b0 + mu_b1*mass_pred_con[k] + mu_b2*1
  
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # varying intercept
  mu_b1 ~ dnorm(-0.25, 0.5)  # varying mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # varying activation energy
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  sigma_b2 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  tau_b2 <- 1/sigma_b2^2
  
  }", fill = TRUE, file = "R/analysis/log-linear model/models/m5_consumption_pred.txt")

model_con = "R/analysis/log-linear model/models/m5_consumption_pred.txt"

jm_con = jags.model(model_con,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3)


# D. EVALUATE MODEL FIT ============================================================
# Plot mean of simulated data vs mean of observed data
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# First convert your matrix 
cs_fit_met = coda.samples(jm_met,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

cs_fit_con = coda.samples(jm_con,
                          variable.names = c("mean_y",
                                             "mean_y_sim", 
                                             "p_mean",
                                             "cv_y",
                                             "cv_y_sim",
                                             "p_cv"), 
                          n.iter = samples, 
                          thin = n.thin)

# Convert to data frames
cs_fit_df_met <- data.frame(as.matrix(cs_fit_met))
cs_fit_df_con <- data.frame(as.matrix(cs_fit_con))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df_met)))

# Metabolism
p1 <- ggplot(cs_fit_df_met, aes(mean_y_sim)) + 
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
  ggtitle("Metabolism") +
  NULL

p2 <- ggplot(cs_fit_df_met, aes(cv_y_sim)) + 
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
p3 <- ggplot(cs_fit_df_con, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df_con$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "C", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df_con$p_mean), digits = 3))) +
  labs(x = "Mean simulated consumption", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  ggtitle("Consumption") +
  NULL

p4 <- ggplot(cs_fit_df_con, aes(cv_y_sim)) + 
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

(p1 + p2) / (p3 + p4)
#ggsave("figures/supp/cv_mean_fit_met_con.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# E. PLOT PREDICTIONS ==============================================================
# JAGS - Nice for summaries and predictions
# Extract the prediction at each x including credible interaval
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# Metabolism
js_met = jags.samples(jm_met, 
                  variable.names = c("pred_warm", "pred_cold"), 
                  n.iter = samples, 
                  thin = n.thin)

# Generate medians and quantiles that can be used for storing info, plotting etc.
# Warm temp:
m_pred_warm <- summary(js_met$pred_warm, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
m_pred_warm_df <- data.frame(lwr_95 = m_pred_warm[1, ],
                             lwr_80 = m_pred_warm[2, ],
                             median = m_pred_warm[3, ],
                             upr_80 = m_pred_warm[4, ],
                             upr_95 = m_pred_warm[5, ],
                             mass = mass_pred_met,
                             temp = -1.5)

# Cold temp:
m_pred_cold <- summary(js_met$pred_cold, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
m_pred_cold_df <- data.frame(lwr_95 = m_pred_cold[1, ],
                           lwr_80 = m_pred_cold[2, ],
                           median = m_pred_cold[3, ],
                           upr_80 = m_pred_cold[4, ],
                           upr_95 = m_pred_cold[5, ],
                           mass = mass_pred_met,
                           temp = 1.5)

# Consumption
js_con = jags.samples(jm_con, 
                      variable.names = c("pred_warm", "pred_cold"), 
                      n.iter = samples, 
                      thin = n.thin)

# Generate medians and quantiles that can be used for storing info, plotting etc.
# Warm temp:
c_pred_warm <- summary(js_con$pred_warm, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_warm_df <- data.frame(lwr_95 = c_pred_warm[1, ],
                             lwr_80 = c_pred_warm[2, ],
                             median = c_pred_warm[3, ],
                             upr_80 = c_pred_warm[4, ],
                             upr_95 = c_pred_warm[5, ],
                             mass = mass_pred_con,
                             temp = -1)

# Cold temp:
c_pred_cold <- summary(js_con$pred_cold, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
c_pred_cold_df <- data.frame(lwr_95 = c_pred_cold[1, ],
                             lwr_80 = c_pred_cold[2, ],
                             median = c_pred_cold[3, ],
                             upr_80 = c_pred_cold[4, ],
                             upr_95 = c_pred_cold[5, ],
                             mass = mass_pred_con,
                             temp = 1)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
#pal <- viridis(option = "magma", n = 10)[c(2, 6)]
pal <- brewer.pal("Dark2", n = 5)[c(1,3)]

m_pdat <- rbind(m_pred_cold_df, m_pred_warm_df)
c_pdat <- rbind(c_pred_cold_df, c_pred_warm_df)

p5 <- ggplot(m_pdat, aes(mass, median, color = factor(temp))) +
  geom_point(data = met, aes(log_mass_norm_ct, log(y)), size = 2.8, shape = 21, 
             alpha = 0.2, color = "white", fill = "grey40") +
  geom_ribbon(data = m_pdat, aes(x = mass, ymin = lwr_95, ymax = upr_95, fill = factor(temp)), 
              size = 0.6, alpha = 0.25, inherit.aes = FALSE)+
  geom_ribbon(data = m_pdat, aes(x = mass, ymin = lwr_80, ymax = upr_80, fill = factor(temp)), 
              size = 0.6, alpha = 0.4, inherit.aes = FALSE) +
  # scale_color_manual(values = pal) +
  # scale_fill_manual(values = pal) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_line(size = 0.6, alpha = 1) +
  theme_classic(base_size = 13) + 
  labs(x = "ln(standardized mass)",
       y = "ln(metabolic rate)") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  theme(aspect.ratio = 3/4) +
  guides(fill = FALSE, color = FALSE) +
  NULL

p5

p6 <- ggplot(c_pdat, aes(mass, median, color = factor(temp))) +
  geom_point(data = con, aes(log_mass_norm_ct, log(y)), size = 2.8, shape = 21, 
             alpha = 0.2, color = "white", fill = "grey40") +
  geom_ribbon(data = c_pdat, aes(x = mass, ymin = lwr_95, ymax = upr_95, fill = factor(temp)), 
              size = 0.6, alpha = 0.25, inherit.aes = FALSE)+
  geom_ribbon(data = c_pdat, aes(x = mass, ymin = lwr_80, ymax = upr_80, fill = factor(temp)), 
              size = 0.6, alpha = 0.4, inherit.aes = FALSE) +
  # scale_color_manual(values = pal) +
  # scale_fill_manual(values = pal) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_line(size = 0.6, alpha = 1) +
  theme_classic(base_size = 13) + 
  labs(x = "ln(standardized mass)",
       y = "ln(maximum consumption rate)",
       color = "Standardized Arrhenius Temperature") +
  guides(fill = FALSE) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  theme(aspect.ratio = 3/4,
        legend.position = "bottom") +
  NULL

p6

p5 / p6
#ggsave("figures/pred_warm_cold_metcon.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)

