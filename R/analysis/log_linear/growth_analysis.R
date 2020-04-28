#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of growth rate as a function of 
# temperature and mass
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit models
#
# D. Model validation (convergence)
# 
# E. Evaluate model fit
#
# F. Plot predictions
# 
# G. Additional calculations on the posterior
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
library(MCMCvis)
library(scales)

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file
dat <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

str(dat)

# Filter data points at below optimum temperatures
dat <- dat %>% filter(above_optimum == "N")

# Rename species factor for JAGS (must be numbered 1:n)
str(dat)
dat$species_n <- as.numeric(as.factor(dat$species))

# Print numeric factor level vs species name
sort(unique(dat$species_ab))
unique(dat$species_n)

# Mean-center predictor variables
dat$log_mass_ct <- dat$log_mass - mean(dat$log_mass)
dat$temp_arr_ct <- dat$temp_arr - mean(dat$temp_arr)

# Use only positive growth rates
dat <- dat %>% filter(y > 0)

# Masses for prediction
mass_pred <-  seq(from = min(dat$log_mass_ct), 
                  to = max(dat$log_mass_ct),
                  length.out = 100)

# Temperature for prediction
temp_pred <- 0 # This means we use mean temperature as it is centered 

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Data in list-format for JAGS
data = list(
  y = log(dat$y), 
  n_obs = length(dat$y), 
  species_n = dat$species_n,
  mass = dat$log_mass_ct,
  temp = dat$temp_arr_ct,
  mass_pred = mass_pred,
  temp_pred = temp_pred
)


# C. FIT MODELS ====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

# Select model with lowest WAIC (see grow_model_selection.R)
model = "R/analysis/JAGS_models/log_linear/selected_models/m1_pred_fit.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    mu_b3 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    sigma_b3 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    mu_b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    sigma_b3 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    mu_b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    sigma_b3 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm = jags.model(model,
                data = data, 
                n.adapt = 5000, 
                n.chains = 3,
                inits = inits)

update(jm, n.iter = n.iter) 


# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm,
                   variable.names = c("b0", "b1", "b2", "b3", 
                                      "mu_b0", "mu_b1", "mu_b2", "mu_b3",
                                      "sigma_b0", "sigma_b1", "sigma_b2", "sigma_b3",
                                      "sigma"), 
                   n.iter = n.iter, 
                   thin = thin)

#** Evaluate convergence ===========================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

#**** Species intercepts ===========================================================
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p2 <- cs_df %>% 
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", 
                          "b0[8]", "b0[9]", "b0[10]", "b0[11]", "b0[12]", "b0[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord1 + pWord2
ggsave("figures/supp/log_linear_model/growth/validation_gro_intercepts.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species mass-effects =========================================================
# Plot posterior densities of species mass-effects
p3 <- cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))


# Traceplot for evaluating chain convergence
p4 <- cs_df %>% 
  filter(Parameter %in% c("b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", "b1[6]", "b1[7]", 
                          "b1[8]", "b1[9]", "b1[10]", "b1[11]", "b1[12]", "b1[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord3 + pWord4
ggsave("figures/supp/log_linear_model/growth/validation_gro_mass.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species temperature-effects ==================================================
# Plot posterior densities of temperature-effects
p5 <- cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p6 <- cs_df %>% 
  filter(Parameter %in% c("b2[1]", "b2[2]", "b2[3]", "b2[4]", "b2[5]", "b2[6]", "b2[7]", 
                          "b2[8]", "b2[9]", "b2[10]", "b2[11]", "b2[12]", "b2[13]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord5 + pWord6
ggsave("figures/supp/log_linear_model/growth/validation_gro_temp.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities of group-level means and standard deviations
p7 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "mu_b3",
                          "sigma_b0", "sigma_b1", "sigma_b2", "sigma_b3")) %>% 
  ggs_density(.) + 
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_density(alpha = 0.05) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Value", y = "Density", fill = "Chain #") +
  guides(color = FALSE, fill = FALSE) +
  NULL
pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))

# Traceplot for evaluating chain convergence
p8 <- cs_df %>% 
  filter(Parameter %in% c("mu_b0", "mu_b1", "mu_b2", "mu_b3",
                          "sigma_b0", "sigma_b1", "sigma_b2", "sigma_b3")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(size = 6)) +
  NULL
pWord8 <- p8 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord7 + pWord8
ggsave("figures/supp/log_linear_model/growth/validation_gro.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p9 <- cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.003) +
  geom_point(size = 2) +
  theme(aspect.ratio = 1)+
  NULL
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
ggsave("figures/supp/log_linear_model/growth/validation_rhat_gro.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(0, 0.04)   # mean of all species intercept
# mu_b1 ~ dnorm(-0.25, 1)  # mean of all species mass-exponent
# mu_b2 ~ dnorm(-0.6, 1)   # mean of all species temperature coefficient
# mu_b3 ~ dnorm(0, 1)      # mean of all species interaction

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# sigma = 1
# tau <- 1/sigma^2

# Define priors for plot
tau <- 1
tau_int <- 0.04

mu_b0 <- rnorm(25000, 0, sqrt(1/tau_int))
mu_b1 <- rnorm(25000, -0.25, sqrt(1/tau))
mu_b2 <- rnorm(25000, -0.6, sqrt(1/tau))
mu_b3 <- rnorm(25000, 0, sqrt(1/tau))

PR <- as.matrix(cbind(mu_b0, mu_b1, mu_b2, mu_b3))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/log_linear_model/growth/validation_prior_post_growth.png", 
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs,
          params = c("mu_b0", "mu_b1", "mu_b2", "mu_b3"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE,
          type = "density")   # removes the trace plot

dev.off()


# E. EVALUATE MODEL FIT ============================================================
# CODA - Nice for getting the raw posteriors
cs_fit = coda.samples(jm,
                      variable.names = c("mean_y",
                                         "mean_y_sim", 
                                         "p_mean",
                                         "cv_y",
                                         "cv_y_sim",
                                         "p_cv"), 
                      n.iter = n.iter, 
                      thin = thin)

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

# Plot mean y and mean cv at each iteration and compare to data
# General formula for number of bins..
n_bins <- round(1 + 3.2*log(nrow(cs_fit_df)))

# Growth
p10 <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3,
           label = paste("P =", round(mean(cs_fit_df$p_mean), digits = 3))) +
  labs(x = "Mean simulated growth", y = "count") +
  coord_cartesian(expand = 0) + 
  NULL
pWord10 <- p10 + theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)


p11 <- ggplot(cs_fit_df, aes(cv_y_sim)) + 
  geom_histogram(bins = n_bins) +
  geom_vline(xintercept = cs_fit_df$cv_y, color = "white", 
             linetype = 2, size = 0.4) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, size = 4, hjust = -0.2, vjust = 3.3, 
           label = paste("P =", round(mean(cs_fit_df$p_cv), digits = 3))) +
  labs(x = "cv simulated growth", y = "") +
  coord_cartesian(expand = 0) +
  NULL

pWord11 <- p11 + theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)
pWord10 + pWord11
ggsave("figures/supp/log_linear_model/growth/fit_gro_mean_cv.png", width = 6.5, height = 6.5, dpi = 600)


# F. PLOT PREDICTIONS ==============================================================
# jags.samples - Nice for summaries and predictions
# Extract the prediction at each x including credible interaval
js = jags.samples(jm, 
                  variable.names = c("pred"), 
                  n.iter = n.iter, 
                  thin = thin)

# Save quantiles
pred <- summary(js$pred, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# Create data frame for the predictions
pred_df <- data.frame(lwr_95 = pred[1, ],
                      lwr_80 = pred[2, ],
                      median = pred[3, ],
                      upr_80 = pred[4, ],
                      upr_95 = pred[5, ],
                      mass = mass_pred)

# This is the mean temperature for the predictions:
# round(mean(dat$temp_arr))  
# round(mean(dat$temp_c)) 

# Plot with mass on x and logarithmic axis instead...
# First add back the mean
pred_df$mass_non_ct <- pred_df$mass + mean(dat$log_mass)

# Then exponentiate
pred_df$mass_g <- exp(pred_df$mass_non_ct)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
# Expand color palette
colourCount = length(unique(dat$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

p12 <- ggplot(pred_df, aes(mass_g, median)) +
  geom_ribbon(data = pred_df, aes(x = mass_g, ymin = lwr_95, ymax = upr_95), 
              size = 2, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = mass_g, ymin = lwr_80, ymax = upr_80), 
              size = 2, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.8, alpha = 0.8) +
  geom_point(data = dat, aes(mass_g, log(dat$y), fill = species_ab),
             size = 2.8, shape = 21, alpha = 0.8, color = "white") +
  scale_fill_manual(values = pal) +
  scale_x_continuous(trans = scales::log_trans(),
                     labels = scales::number_format(accuracy = 1)) +
  guides(fill = FALSE) +
  labs(x = "mass [g]",
       y = "ln(growth rate [%/day])") +
  annotate("text", 0.2, 3.2, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  # annotate("text", -Inf, Inf, label = "A", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) + # This solution doesn't play with the ln axis...
  NULL
pWord12 <- p12 + theme_classic() + theme(text = element_text(size = 12))

# Add posterior distributions of parameters
cs = coda.samples(jm,
                  variable.names = c(
                    "mu_b0",
                    "mu_b1",
                    "mu_b2",
                    "mu_b3"), 
                  n.iter = n.iter, 
                  thin = thin)

cs_df <- ggs(cs)

# Posterior of parameters
color_scheme_set("gray")
sum_dat <- data.frame(summary(cs)[1])

# Mass-coefficient
p13 <- cs %>% 
  mcmc_dens(pars = "mu_b1") +
  geom_vline(xintercept = sum_dat[2, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass-exponent") +
  NULL
pWord13 <- p13 + theme_classic() + theme(text = element_text(size = 12))

# Temperature-coefficient
p14 <- cs %>% 
  mcmc_dens(pars = "mu_b2") +
  geom_vline(xintercept = sum_dat[3, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "C", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Temperature coefficient") +
  NULL
pWord14 <- p14 + theme_classic() + theme(text = element_text(size = 12))


# Mass-temperature interaction
p15 <- cs %>%
  mcmc_dens(pars = "mu_b3") +
  geom_vline(xintercept = 0, color = "red", size = 0.6, linetype = 1) +
  geom_vline(xintercept = sum_dat[4, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-0.2, 0.15)) +
  annotate("text", -Inf, Inf, label = "D", size = 4,
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "M*T interaction") +
  NULL
pWord15 <- p15 + theme_classic() + theme(text = element_text(size = 12))

# Plot all together
pWord12 / (pWord13 + pWord14 + pWord15) + plot_layout(heights = c(2.5, 1, 1, 1))

ggsave("figures/pred_gro.png", width = 6.5, height = 6.5, dpi = 600)


# G. ADDITINAL CALCULATIONS ON THE POSTERIOR =======================================
# Calculate the proportion of the posterior of activation energy that is less than zero
js = jags.samples(jm, 
                  variable.names = c("mu_b2", "mu_b3"), 
                  n.iter = n.iter, 
                  thin = thin)
 
ecdf(js$mu_b2)(-0.65) 
# [1] 0.8321111

ecdf(js$mu_b3)(0) # how much is below?

# How much does the mass exponent decline per change in unit T?
#summary(cs)

# Coefficient is 0.016
# head(dat)
# dat$b_a <- 0.016 * dat$temp_norm_arr_ct
# 
# summary(lm(b_a ~ temp_norm_arr_ct, data = dat))

# Now fit the same exponents to C
#summary(lm(b_a ~ temp_norm, data = dat))

