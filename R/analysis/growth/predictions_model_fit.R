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
dat <- read.csv("data/growth_analysis.csv")

str(dat)

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_optimum == "N")

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Mass-range used for prediction
mass_pred = seq(from = min(dat$log_mass_norm_ct), 
                to = max(dat$log_mass_norm_ct),
                length.out = 100)

# Filter only positive growth rates
dat <- dat %>% filter(., G > 0)

# Data in list-format for JAGS
data = list(
  y = log(dat$G), 
  n_obs = length(dat$G), 
  species_n = dat$species_n,
  mass = dat$log_mass_norm_ct,
  temp = dat$temp_norm_arr_ct,
  mass_pred = mass_pred
)

summary(lm(log(G) ~ temp_norm_arr_ct * log_mass_norm_ct, data = dat))

# C. FIT MODELS ====================================================================
# Refit chosen models from the model selection part
# Need to modify them to track new variables for assessing model fit
#**** Growth (M1) ==============================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +               # varying intercept 
      b1[species_n[i]]*mass[i] +       # varying mass-exponent
      b2[species_n[i]]*temp[i] +       # varying activation energy
      b3[species_n[i]]*mass[i]*temp[i] # varying M*T interaction
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
  for(k in 1:length(mass_pred)){
      
    pred_warm[k] <- mu_b0 + mu_b1*mass_pred[k] + mu_b2*-1 + mu_b3*mass_pred[k]*-1
    pred_cold[k] <- mu_b0 + mu_b1*mass_pred[k] + mu_b2*1 + mu_b3*mass_pred[k]*1
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  mu_b0 ~ dnorm(0, 0.5)      # global mean
  mu_b1 ~ dnorm(-0.25, 0.5)  # global mass-exponent
  mu_b2 ~ dnorm(-0.6, 0.5)   # global activation energy
  mu_b3 ~ dnorm(0, 0.5)      # global interaction
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
  
  }", fill = TRUE, file = "R/analysis/growth/models/m1_growth_pred.txt")

gro_model = "R/analysis/growth/models/m1_growth_pred.txt"

jm_gro = jags.model(gro_model,
                    data = data, 
                    n.adapt = 5000, 
                    n.chains = 3)


# D. EVALUATE MODEL FIT ============================================================
# Plot mean of simulated data vs mean of observed data
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# First convert your matrix 
cs_fit = coda.samples(jm_gro,
                      variable.names = c("mean_y",
                                         "mean_y_sim", 
                                         "p_mean",
                                         "cv_y",
                                         "cv_y_sim",
                                         "p_cv"), 
                      n.iter = samples, 
                      thin = n.thin)

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df)))

# Growth
p1 <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df$p_mean), digits = 3))) +
  labs(x = "Mean simulated growth", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  #ggtitle("Growth rate") +
  NULL

p2 <- ggplot(cs_fit_df, aes(cv_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df$cv_y, color = "white", 
             linetype = 2, size = 0.4) +
  theme_classic(base_size = 11) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3, 
           label = paste("P =", round(mean(cs_fit_df$p_cv), digits = 3))) +
  labs(x = "cv simulated growth", y = "") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) +
  NULL


p1 + p2
#ggsave("figures/supp/cv_mean_fit_gro.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# E. PLOT PREDICTIONS ==============================================================
# JAGS - Nice for summaries and predictions
# Extract the prediction at each x including credible interaval
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# Metabolism
js_gro = jags.samples(jm_gro, 
                      variable.names = c("pred_warm", "pred_cold"), 
                      n.iter = samples, 
                      thin = n.thin)

# Generate medians and quantiles that can be used for storing info, plotting etc.
# Warm temp:
pred_warm <- summary(js_gro$pred_warm, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
pred_warm_df <- data.frame(lwr_95 = pred_warm[1, ],
                           lwr_80 = pred_warm[2, ],
                           median = pred_warm[3, ],
                           upr_80 = pred_warm[4, ],
                           upr_95 = pred_warm[5, ],
                           mass = mass_pred,
                           temp = -1)

# Cold temp:
pred_cold <- summary(js_gro$pred_cold, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
pred_cold_df <- data.frame(lwr_95 = pred_cold[1, ],
                           lwr_80 = pred_cold[2, ],
                           median = pred_cold[3, ],
                           upr_80 = pred_cold[4, ],
                           upr_95 = pred_cold[5, ],
                           mass = mass_pred,
                           temp = 1)


# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
#pal <- viridis(option = "magma", n = 10)[c(2, 6)]
pal <- brewer.pal("Set1", n = 5)

pdat <- rbind(pred_cold_df, pred_warm_df)

p3 <- ggplot(pdat, aes(mass, median, color = factor(temp), fill = factor(temp))) +
  geom_ribbon(data = filter(pdat, temp == 1), aes(x = mass, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.2, inherit.aes = FALSE, fill = pal[2]) +
  geom_ribbon(data = filter(pdat, temp == 1), aes(x = mass, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = pal[2]) +
  geom_ribbon(data = filter(pdat, temp == -1), aes(x = mass, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.2, inherit.aes = FALSE, fill = pal[1]) +
  geom_ribbon(data = filter(pdat, temp == -1), aes(x = mass, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = pal[1]) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(data = dat, aes(log_mass_norm_ct, log(G)),
             size = 2.8, shape = 21, alpha = 0.8, color = "white", fill = "grey40") +
  scale_color_manual(values = pal) +
  theme_classic(base_size = 11) + 
  labs(x = "ln(standardized mass)",
       y = "ln(growth rate [%/day])",
       color = "Temperature\n(centered\nArrhenius)") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

p3

# Add posterior distributions of parameters
cs = coda.samples(jm_gro,
                   variable.names = c(
                     "mu_b1",
                     "mu_b2",
                     "mu_b3"), 
                   n.iter = samples, 
                   thin = n.thin)


cs_df <- ggs(cs)

# Posterior of parameters
color_scheme_set("gray")
sum_dat <- data.frame(summary(cs)[1])

# Mass-coefficient
p4 <- cs %>% 
  mcmc_dens(pars = "mu_b1") +
  theme_classic(base_size = 11) + 
  geom_vline(xintercept = sum_dat[1, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass-exponent") +
  NULL

# Temperature-coefficient
p5 <- cs %>% 
  mcmc_dens(pars = "mu_b2") +
  theme_classic(base_size = 11) + 
  geom_vline(xintercept = sum_dat[2, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "C", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Activation energy") +
  NULL

# Mass-temperature interaction
p6 <- cs %>% 
  mcmc_dens(pars = "mu_b3") +
  theme_classic(base_size = 11) + 
  geom_vline(xintercept = 0, color = "red", size = 0.6, linetype = 2) +
  geom_vline(xintercept = sum_dat[3, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "D", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "M*T interaction") +
  NULL

# Plot all together
p3 / (p4 + p5 + p6) + plot_layout(ncol = 1, heights = c(2.5, 1, 1))

#ggsave("figures/pred_warm_cold_gro.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)
