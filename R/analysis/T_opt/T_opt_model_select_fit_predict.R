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
# D. Evaluate model fit
#
# E. Evaluate model convergence
#
# F. Plot predicted mass-scaling slopes in different temperatures
#
# G. Plot temperature data
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
dat <- read.csv("data/topt_analysis.csv")

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

# Data in list-format for JAGS
data = list(
  y = dat$opt_temp_c_ct, 
  n_obs = length(dat$opt_temp_c_ct), 
  species_n = dat$species_n,
  mass = dat$log_mass_norm_ct,
  mass_pred = mass_pred
)

summary(lm(opt_temp_c_ct ~ log_mass_norm_ct, data = dat))

# Check distribution of relative size range
p1 <- ggplot(dat, aes(mass_norm)) +
  geom_histogram()

p2 <- ggplot(dat, aes(mass)) +
  geom_histogram()

p1/p2

library(tidylog)
t <- filter(dat, mass_norm < 0.0100001)

length(dat$opt_temp_c_ct)

ggplot(dat, aes(log_mass_norm_ct, opt_temp_c_ct)) + geom_point()
ggplot(dat, aes((mass_norm - mean(mass_norm)), opt_temp_c_ct)) + geom_point()
ggplot(dat, aes(mass_norm, opt_temp_c_ct)) + geom_point()
ggplot(dat, aes(mass_norm, opt_temp_c_ct)) + geom_point() + coord_cartesian(xlim = c(0, 0.024))

# C. MODEL SELECTION ===============================================================
#**** Random intercept & slope (M1) ================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +          # varying intercept 
      b1[species_n[i]]*mass[i]    # varying mass-exponent

  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  
  # Predictions
  for(k in 1:length(mass_pred)){
      
    pred[k] <- mu_b0 + mu_b1*mass_pred[k]
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  mu_b0 ~ dnorm(0, 5)    
  mu_b1 ~ dnorm(0, 5)    
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  sigma_b1 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  tau_b1 <- 1/sigma_b1^2
  
  }", fill = TRUE, file = "R/analysis/T_opt/models/m1_opt_pred.txt")

opt_model = "R/analysis/T_opt/models/m1_opt_pred.txt"

jm1 = jags.model(opt_model,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm1, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj1 = jags.samples(jm1, 
                   variable.names = c("pd", "log_pd"), 
                   n.iter = 10000, 
                   thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd1 <- -2*sum(log(summary(zj1$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC1 <- sum((summary(zj1$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC1 <- lppd1 + 2*pd.WAIC1

c(pd.WAIC1, WAIC1)
waic_m1 <- WAIC1
waic_m1


#**** Random intercept (M2) ================================================
cat(
  "model{
  
  for(i in 1:n_obs){
    # Simulate for comparison with data
    y_sim[i] ~ dnorm(mu[i], tau)
    
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- 
      b0[species_n[i]] +   # varying intercept 
      b1*mass[i]           # non-varying mass-exponent

  # Add log likelihood computation for each observation
  pd[i] <- dnorm(y[i], mu[i], tau)
  
  # Calculates the log PPD
  log_pd[i] <- log(dnorm(y[i], mu[i], tau))
  }
  
  # Second level (species-level effects)
  for(j in 1:max(species_n)){
    b0[j] ~ dnorm(mu_b0, tau_b0)
  }
  
  # Predictions
  for(k in 1:length(mass_pred)){
      
    pred[k] <- mu_b0 + b1*mass_pred[k]
    
  } 

  # Model fit
  mean_y <- mean(y[])
  mean_y_sim <- mean(y_sim[])
  p_mean <- step(mean_y_sim - mean_y) # Proportion of data above and below 
  
  cv_y <- sd(y[])/mean(y[])
  cv_y_sim <- sd(y_sim[])/max(0.0000001, mean(y_sim[])) # Not to divide by 0
  p_cv <- step(cv_y_sim - cv_y)

  #-- Priors	
  mu_b0 ~ dnorm(0, 5)    
  b1 ~ dnorm(0, 5)       
  sigma ~ dunif(0, 10) 
  sigma_b0 ~ dunif(0, 10)
  tau <- 1/sigma^2
  tau_b0 <- 1/sigma_b0^2
  
  }", fill = TRUE, file = "R/analysis/T_opt/models/m2_opt_pred.txt")

opt_model = "R/analysis/T_opt/models/m2_opt_pred.txt"

jm2 = jags.model(opt_model,
                 data = data, 
                 n.adapt = 5000, 
                 n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm2, n.iter = burn.in) 

# Monitor the likelihood to calculate WAIC
zj2 = jags.samples(jm2, 
                  variable.names = c("pd", "log_pd"), 
                  n.iter = 10000, 
                  thin = 1)

# Calculate model fit by summing over the log of means of the posterior distribution of 
# the PPD and multiply by -2 (i.e. negative log likelihood).
lppd2 <- -2*sum(log(summary(zj2$pd, mean)$stat))

# Calculate penalty (i.e. the number of parameters) as the variance of
# the log of PPD. Do this by squaring the standard deviation.
pd.WAIC2 <- sum((summary(zj2$log_pd, sd)$stat)^2) # Penalty

# WAIC = model fit + 2*penalty
WAIC2 <- lppd2 + 2*pd.WAIC2

c(pd.WAIC2, WAIC2)
waic_m2 <- WAIC2

# Compare both WAIC
waic_m1
waic_m2



# D. EVALUATE MODEL FIT ============================================================
# Plot mean of simulated data vs mean of observed data
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# Get main parameters
cs_fit_par = coda.samples(jm2,
                          variable.names = c("mu_b0", "b1"), 
                          n.iter = samples, 
                          thin = n.thin)

summary(cs_fit_par)

# First convert your matrix 
cs_fit = coda.samples(jm2,
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
  # annotate("text", -Inf, Inf, label = "A", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df$p_mean), digits = 3))) +
  labs(x = "Mean simulated growth", y = "count") +
  theme(aspect.ratio = 1) +
  coord_cartesian(expand = 0) + 
  ggtitle("Optimum growth temperature") +
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


#p1 + p2
p1
#ggsave("figures/supp/cv_mean_fit_T_opt.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm", dpi = 300)


# E. MODEL VALIDATION ==============================================================
#** Sample from the posterior ======================================================
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm2,
                   variable.names = c("b1", "b0", "mu_b0", "sigma_b0", "sigma"), 
                   n.iter = samples, 
                   thin = n.thin)

summary(cs) # Get the mean estimate and SE and 95% CIs

cs_df <- data.frame(summary(cs)[1])
cs_df$Parameter <- row.names(cs_df)


#** Evaluate convergence ===========================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

# Plot posterior densities of species intercepts
p1 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
p2 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
p1+p2
#ggsave("figures/supp/model_validation_t_opt_intercepts.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# Plot posterior densities of params
p1 <- cs_df %>% 
  filter(Parameter %in% c("b1", "mu_b0", "sigma_b0", "sigma")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL

# Traceplot for evaluating chain convergence
p2 <- cs_df %>% 
  filter(Parameter %in% c("b1", "mu_b0", "sigma_b0", "sigma")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  theme_classic(base_size = 11) + 
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
p1+p2
#ggsave("figures/supp/model_validation_t_opt.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)

# Rhat
cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.008) +
  theme_classic(base_size = 11) +
  geom_point(size = 2) +
  theme(aspect.ratio = 1)+
  NULL
#ggsave("figures/supp/rhat_t_opt.pdf", plot = last_plot(), scale = 1, width = 14, height = 14, units = "cm", dpi = 300)



# F. PLOT PREDICTIONS ==============================================================
# JAGS - Nice for summaries and predictions
# Extract the prediction at each x including credible interaval
samples = 10000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

js = jags.samples(jm2, 
                  variable.names = c("pred"), 
                  n.iter = samples, 
                  thin = n.thin)

# Generate medians and quantiles that can be used for storing info, plotting etc.
# Warm temp:
pred <- summary(js$pred, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
pred_df <- data.frame(lwr_95 = pred[1, ],
                      lwr_80 = pred[2, ],
                      median = pred[3, ],
                      upr_80 = pred[4, ],
                      upr_95 = pred[5, ],
                      mass = mass_pred)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
colourCount = length(unique(dat$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p3 <- ggplot(pred_df, aes(mass, median)) +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = mass, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 1, alpha = 1, col = "black") +
  geom_point(data = dat, aes(log_mass_norm_ct, opt_temp_c_ct, fill = species, size = mass),
             #size = 3.5, 
             shape = 21, 
             alpha = 0.8, 
             color = "white"
             ) +
  theme_classic(base_size = 13) + 
  scale_fill_manual(values = pal) +
  scale_size(range = c(2, 8), breaks = c(0, 1, 10, 100, 1000)) +
  theme(aspect.ratio = 4/5,
        #legend.position = c(0.11, 0.22),
        legend.title = element_text(size = 10)) +
  guides(fill = FALSE,
         size = guide_legend(override.aes = list(fill = "black",
                                                 color = "black"))) +
  labs(x = "ln(rescaled mass ct)",
       y = "Rescaled optimum growth temperature",
       size = "Mass [g]") +
  # annotate("text", -Inf, Inf, label = "A", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

p3

#ggsave("figures/topt_scatter.pdf", plot = last_plot(), scale = 1, width = 15, height = 15, units = "cm", dpi = 300)

# Add posterior distributions of parameters
cs = coda.samples(jm2,
                  variable.names = "b1", 
                  n.iter = samples, 
                  thin = n.thin)


cs_df <- ggs(cs)

# Posterior of parameters
color_scheme_set("gray")
sum_dat <- data.frame(summary(cs)[1])

# Mass-coefficient
p4 <- cs %>% 
  mcmc_dens(pars = "b1") +
  theme_classic(base_size = 14) + 
  geom_vline(xintercept = sum_dat[1, 1], color = "white", size = 0.6, linetype = 2) +
  geom_vline(xintercept = 0, color = "red", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "B", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "T_opt mass-coefficient") +
  NULL
p4

# Calculate the proportion of the posterior that is less than zero
js = jags.samples(jm2, 
                  variable.names = c("b1"), 
                  n.iter = samples, 
                  thin = n.thin)

ecdf(js$b1)(0) # We are 99% certain the slope is smaller than 0

#ggsave("figures/supp/T_opt_mass_posterior.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)




# G. PLOT TEMPERATURE DATA =========================================================
# Plot optimum growth compared to experimental and environmental temperature
# Read in growth data as that has experimental temperature
dat2 <- read.csv("data/growth_analysis.csv")

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat2$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat2$species )

dat2$species_ab <- paste(sp1, sp2, sep = ".")

dat2$env_temp_max <- as.numeric(as.character(dat2$env_temp_max))
dat2$env_temp_min <- as.numeric(as.character(dat2$env_temp_min))

specs <- sort(unique(dat2$species_ab))
dat2$env_temp_min

# Separate temperature sources
dat2$temp_source <- 1

dat2$env_temp_min[is.na(dat2$env_temp_min)] <- -9

dat2$temp_source <- ifelse(dat2$env_temp_min == -9,
                           2,
                           dat2$temp_source)

unique(filter(dat2, env_temp_min == -9)$species_ab)

# Beacuse I didn't put the range of the "preffered" temperature, I put them here now
dat2$env_temp_min <- ifelse(dat2$species_ab == "P.yokohamae",
                            3,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "P.yokohamae",
                            24,
                            dat2$env_temp_max)
dat2$env_temp_min <- ifelse(dat2$species_ab == "C.lumpus",
                            0.6,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "C.lumpus",
                            11.4,
                            dat2$env_temp_max)
dat2$env_temp_min <- ifelse(dat2$species_ab == "P.olivaceus",
                            8.6,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "P.olivaceus",
                            25,
                            dat2$env_temp_max)
dat2$env_temp_min <- ifelse(dat2$species_ab == "H.hippoglossus",
                            0.4,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "H.hippoglossus",
                            7.9,
                            dat2$env_temp_max)
dat2$env_temp_min <- ifelse(dat2$species_ab == "S.maximus",
                            5.9,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "S.maximus",
                            11.9,
                            dat2$env_temp_max)
dat2$env_temp_min <- ifelse(dat2$species_ab == "A.minor",
                            0.6,
                            dat2$env_temp_min)
dat2$env_temp_max <- ifelse(dat2$species_ab == "A.minor",
                            7.6,
                            dat2$env_temp_max)


# Prepare experimental data
sub <- dat2 %>% 
  select(temp_c, median_temp, species_ab, log_mass, temp_source) %>% 
  gather(source, temp, 1:2) %>%
  ungroup()

# Raw T_opt data
# CHECK THIS
sub2 <- dat %>% 
  mutate(source = "Optimum (data)",
         temp = opt_temp_c,
         log_mass = log(mass)) %>% 
  select(temp, species_ab, source, log_mass)

head(sub2)
head(sub)

sub <- bind_rows(sub2, sub)

sub$source <- ifelse(sub$source == "temp_c",
                     "Experimental\ntemperature",
                     sub$source)

sub$source <- ifelse(sub$source == "median_temp",
                     "Mid-point env. temperature",
                     sub$source)

  
# For ordering overplotting data points in plot...
sub$source2 <- 1
sub$source2 <- ifelse(sub$source == "Optimum (data)", 3, sub$source2)
sub$source2 <- ifelse(sub$source == "Mid-point env. temperature", 2, sub$source2)

str(dat2)

# need to find another way of keeong variables other than calculating mean...
dat3 <- dat2 %>% 
  group_by(species_ab) %>% 
  tidyr::drop_na(env_temp_max) %>% 
  tidyr::drop_na(env_temp_min) %>% 
  summarise(upper = mean(env_temp_max), 
            lower = mean(env_temp_min),
            temp = mean(temp_c),
            temp_source = mean(temp_source)) %>% 
  arrange(species_ab)

dat3

# Filter T_opt data
sub3 <- dat %>%
  group_by(species_ab) %>%
  summarize(temp_opt   = mean(opt_temp_c),
            source     = "Optimum (data)") %>%
  ungroup()

dat4 <- full_join(sub3, dat3) %>% drop_na() %>% arrange()

dat4

sub %>% filter(source == "Optimum (data)") %>% arrange(temp)

# Need to find a common way of sorting... currently I have a mix of data in separate files..

unique(dat4$species_ab)
dat4$sort <- 1
dat4$sort <- ifelse(dat4$species_ab == "A.minor", 2, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "H.hippoglossus", 3, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "C.lumpus", 4, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "S.salar", 5, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "G.morhua", 6, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "S.alpinus", 7, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "S.maximus", 8, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "P.yokohamae", 9, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "P.olivaceus", 10, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "P.fulvidraco", 11, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "L.calcarifer", 12, dat4$sort)
dat4$sort <- ifelse(dat4$species_ab == "R.canadum", 13, dat4$sort)

sub$sort <- 1
sub$sort <- ifelse(sub$species_ab == "A.minor", 2, sub$sort)
sub$sort <- ifelse(sub$species_ab == "H.hippoglossus", 3, sub$sort)
sub$sort <- ifelse(sub$species_ab == "C.lumpus", 4, sub$sort)
sub$sort <- ifelse(sub$species_ab == "S.salar", 5, sub$sort)
sub$sort <- ifelse(sub$species_ab == "G.morhua", 6, sub$sort)
sub$sort <- ifelse(sub$species_ab == "S.alpinus", 7, sub$sort)
sub$sort <- ifelse(sub$species_ab == "S.maximus", 8, sub$sort)
sub$sort <- ifelse(sub$species_ab == "P.yokohamae", 9, sub$sort)
sub$sort <- ifelse(sub$species_ab == "P.olivaceus", 10, sub$sort)
sub$sort <- ifelse(sub$species_ab == "P.fulvidraco", 11, sub$sort)
sub$sort <- ifelse(sub$species_ab == "L.calcarifer", 12, sub$sort)
sub$sort <- ifelse(sub$species_ab == "R.canadum", 13, sub$sort)

# Plot
sub6 <- sub %>% arrange(source2)

pal2 <- RColorBrewer::brewer.pal("Dark2", n = 5)



######## LESS DIFERRENCE BETWEEN n +  DO OPEN VS DOLID CIRCLES FOR GREEN TEMP SOURCE



ggplot() +
  geom_point(data = sub6, aes(x = reorder(species_ab, sort), y = temp, 
                              fill = source, alpha = source, 
                              #shape = factor(temp_source)
                              ),
             size = 3.5, 
             shape = 21, 
             color = "white", 
             position = position_dodge(width = 1)) +
  geom_point(data = filter(sub6, temp_source == 1 & source == "Mid-point env. temperature"), 
             aes(x = reorder(species_ab, sort), y = temp),
             size = 3.5, shape = 25, fill = pal2[1], color = "white", position = position_dodge(width = 1)) +
  geom_errorbar(data = dat4, aes(x = reorder(species_ab, sort), 
                                 ymin = lower, ymax = upper, linetype = factor(temp_source)),
                width = 0.5, color = pal2[1], shape = 23) +
  geom_count(data = filter(sub6, source == "Optimum (data)"),
             aes(x = reorder(species_ab, sort), y = temp),
             alpha = 0.8, shape = 21, color = pal2[2], fill = pal2[2]) +
  scale_fill_manual(values = c("grey75", pal2[1], pal2[2], pal2[2])) +
  scale_alpha_manual(values = c(0.8, 0.8, 0)) +
  scale_linetype_manual(values = c(1,1)) +
  scale_size_area(max_size = 4) +
  #scale_shape_manual(values = c(23, 21)) +
  theme_classic(base_size = 16) +
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,1),
                                                 color = "white")),
         linetype = FALSE,
         shape = FALSE
         ) +
  xlab("") + 
  scale_size_area(max_size = 5) +
  ylab(expression(paste("Temperature [", degree*C, "]"))) + 
  coord_flip() +
  theme(
        #legend.position = c(0.85, 0.26),
        axis.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 10),
        aspect.ratio = 6/7
        ) +
  NULL 

#ggsave("figures/opt_env_exp.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm", dpi = 300)




# I need to go from long to wide in sub, and then do group by and calculate - by species - the mean distance from environmental midpoint

sub3 <- dat %>%
  group_by(species_ab) %>%
  summarize(temp = mean(opt_temp_c)) %>%
  arrange(species_ab) %>% 
  ungroup()

sub3

sub4 <- dat2 %>%
  group_by(species_ab) %>%
  summarize(median_temp = mean(median_temp)) %>%
  arrange(species_ab) %>% 
  ungroup()

sub4

sub3$median_temp <- sub4$median_temp

sub3$mean_t_diff <- sub3$temp - sub3$median_temp

mean(sub3$mean_t_diff)
sd(sub3$mean_t_diff)

# head(sub)
# 
# sub_wide <- sub %>% 
#   ungroup() %>% 
#   filter(source %in% c("Optimum (data)", "Mid-point env. temperature")) %>% 
#   select(temp, species_ab, source) %>% 
#   pivot_wider(values_from = temp, names_from = source)
# 
# sub_wide
