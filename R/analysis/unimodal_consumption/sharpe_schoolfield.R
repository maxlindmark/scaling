#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit mixed-effects Sharpe-Schoolfield model to unimodal consumption data
#   following the parameterization of Padfield et al 2020.
# 
# - Not all variables calculated in the beginning are currently used, some where
#   calculated for data exploration
# 
# A. Load libraries
#
# B. Read data
#
# C. Fit model
# 
# D. Model validation (convergence, fit, residuals)
#
# E. Plot predictions
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
library(tidylog)
library(rTPC)
library(nls.multstart)
library(broom)

# sessionInfo()
# other attached packages:
# [1] broom_0.7.10        nls.multstart_1.2.0 rTPC_1.0.0          tidylog_1.0.2      
# [5] scales_1.1.1        MCMCvis_0.14.0      bayesplot_1.7.2     patchwork_1.1.1    
# [9] viridis_0.5.1       viridisLite_0.4.0   magrittr_2.0.1      readxl_1.3.1       
# [13] RCurl_1.98-1.5      ggmcmc_1.4.1        ggplot2_3.3.5       tidyr_1.1.4        
# [17] dplyr_1.0.7         RColorBrewer_1.1-2  rjags_4-10          coda_0.19-4   


# B. READ IN DATA ==================================================================
# Read in data file
con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

nrow(con)

# Which species have data above optimum?
spec <- unique(filter(con, above_peak_temp == "Y"))$species

con <- con %>% filter(species %in% spec) %>% droplevels()

nrow(con)

# Center temperature relative to mean in environment by species
con$temp_env_ct <- con$temp_c - con$median_temp

# Mass-normalize consumption (whole organism exponent based on log-linear consumption model)
con$y_spec <- con$y / con$mass_g^0.63

# Express consumption rate as fraction of mean within species
con <- con %>%
  group_by(species_ab) %>%
  mutate(y_mean_species = mean(y_spec),
         y_ct = y_spec / y_mean_species) %>%
  ungroup()

# Standardize to modeled peak temperature within each species
# Loop through each species and fit a model with a quadratic temperature term
# Then filter the temperature at peak rate, which we will subtract from y further down
tmp_dat <- data.frame()
datalist <- list()

for(i in unique(con$species_ab)) {
  
  tmp_dat <- con %>% filter(species_ab == i)
  
  # Get starting values
  start_vals <- get_start_vals(con$temp_c, con$y_ct, model_name = 'sharpeschoolhigh_1981')
  
  # Get limits
  low_lims <- get_lower_lims(con$temp_c, con$y_ct, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(con$temp_c, con$y_ct, model_name = 'sharpeschoolhigh_1981')
  
  # Fit Sharpe-Schoolfield
  fit <- nls_multstart(y_ct ~ sharpeschoolhigh_1981(temp = temp_c, r_tref, e, eh, th, tref = 0),
                       data = tmp_dat,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  # Predict
  new_data <- data.frame(temp_c = seq(min(tmp_dat$temp_c), max(tmp_dat$temp_c), length.out = 50))
  preds <- augment(fit, newdata = new_data)
  
  # Extract T_opt and add in species name
  topt <- calc_params(fit) %>%
    mutate_all(round, 2) %>% 
    select(topt) %>% 
    rename("peak_temp_c_model" = "topt") %>% 
    mutate(species_ab = i)

  # Plot data and model fit
  ggplot(tmp_dat, aes(temp_c, y_ct)) +
    geom_point() +
    geom_line(aes(temp_c, .fitted), preds, col = 'blue') +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature [ºC]',
         y = 'Maximum consumption rate') %>% 
    ggtitle(i) 
    
  ggsave(paste("figures/supp/unimodal_consumption/single_species_sharpe_schoolfield_models/", i, ".png", sep = ""))
  
  datalist[[i]] <- topt
  
}

est <- dplyr::bind_rows(datalist)

# Add in modeled peak temperature
colnames(con)

con <- left_join(con, est)

head(con)

# Center temperature now with modeled peak temperature within species
con$peak_temp_c_model_ct <- con$temp_c - con$peak_temp_c_model

# Filter species without a peak temperature
con <- con %>% filter(!species_ab == "P.yokohamae")

nrow(con)

# Plot data
ggplot(con, aes(peak_temp_c_model_ct, y_ct)) +
  geom_point() +
  theme_bw(base_size = 12) +
  stat_smooth() +
  labs(x = 'Temperature_centered (ºC)',
       y = 'Consumption rate')

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Rename species factor for JAGS (must be numbered 1:n)
con$species_n <- as.numeric(as.factor(con$species_ab))

# Arrange by species name so there's no confusion between the index and the name
con <- con %>% arrange(species_n)

# Data in list-format for JAGS
data = list(
  y = con$y_ct,
  temp = con$peak_temp_c_model_ct, 
  n_obs = length(con$y_ct), 
  species_n = as.numeric(as.factor(con$species_ab)),
  temp_pred = seq(from = min(con$peak_temp_c_model_ct),
                  max(con$peak_temp_c_model_ct),
                  length.out = 100)
)


# C. FIT MODEL =====================================================================
# Some settings:
burn.in <- 15000 # Length of burn-in
n.iter <- 15000  # Number of samples
thin <- 5        # Save every 5th sample

model = "JAGS_models/unimodal_consumption/sharpe_school.txt"

# Manually set initial values, because otherwise all the chains get the same
inits = list(
  list(
    mu_b0 = 0.1,
    sigma_b0 = 0.1,
    mu_E = 0.1,
    Eh = 0.1,
    Th = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2 # This is to reproduce the same samples
  ),
  list(
    mu_b0 = 1,
    sigma_b0 = 1,
    mu_E = 1,
    Eh = 1,
    Th = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    sigma_b0 = 2,
    mu_E = 2,
    Eh = 2,
    Th = 2,
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
                   variable.names = c("mu_b0", "sigma_b0", "mu_E", "Eh", "Th", "E", "b0"), 
                   n.iter = n.iter, 
                   thin = thin)

summary(cs)

# 2. Quantiles for each variable:
#   
#           2.5%   25%    50%    75%    97.5%
# E[1]      0.8534 0.9628 1.0250 1.0885 1.2123
# E[2]      0.5252 0.6162 0.6708 0.7279 0.8462
# E[3]      0.5742 0.7250 0.8158 0.9062 1.0967
# E[4]      0.4592 0.5516 0.6022 0.6560 0.7665
# E[5]      0.3814 0.5340 0.6115 0.6923 0.8533
# E[6]      0.6439 0.7968 0.8838 0.9715 1.1376
# E[7]      0.6441 0.7814 0.8537 0.9302 1.0871
# E[8]      0.5581 0.7123 0.7978 0.8853 1.0604
# E[9]      0.3370 0.3951 0.4272 0.4605 0.5327
# E[10]     0.4639 0.5743 0.6390 0.7019 0.8294
# Eh        1.6753 1.8099 1.8849 1.9578 2.0987
# Th       -0.8571 0.1837 0.7463 1.2769 2.3665
# b0[1]     0.4756 0.5516 0.5923 0.6338 0.7179
# b0[2]     0.9363 1.0478 1.1072 1.1665 1.2809
# b0[3]     0.4407 0.5710 0.6419 0.7165 0.8627
# b0[4]     0.7769 0.8535 0.8945 0.9363 1.0172
# b0[5]     0.5651 0.6753 0.7352 0.8020 0.9353
# b0[6]     0.3615 0.4538 0.5114 0.5736 0.7015
# b0[7]     0.4848 0.5998 0.6658 0.7318 0.8585
# b0[8]     0.4223 0.5392 0.6026 0.6716 0.8108
# b0[9]     1.1216 1.2027 1.2473 1.2922 1.3844
# b0[10]    0.6886 0.7893 0.8441 0.9015 1.0085
# mu_E      0.5353 0.6591 0.7275 0.7962 0.9408
# mu_b0     0.5769 0.7227 0.7880 0.8521 0.9908
# sigma_b0  0.1658 0.2301 0.2776 0.3373 0.5227


# Evaluate convergence =============================================================
# Convert to ggplottable data frame
cs_df <- ggs(cs)

#**** Species activation energies ==================================================
# Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p1 <- cs_df %>% 
  filter(Parameter %in% c("E[1]", "E[2]", "E[3]", "E[4]", "E[5]", "E[6]", "E[7]", "E[8]",
                          "E[9]", "E[10]", "E[11]")) %>% 
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
  filter(Parameter %in% c("E[1]", "E[2]", "E[3]", "E[4]", "E[5]", "E[6]", "E[7]", "E[8]",
                          "E[9]", "E[10]", "E[11]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord1 + pWord2
ggsave("figures/supp/unimodal_consumption/validation_SharpeSchool_E.png", width = 6.5, height = 6.5, dpi = 600)


#**** Species b0 ===================================================================
#Plot posterior densities of species intercepts
unique(cs_df$Parameter)

p3 <- cs_df %>%
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", "b0[8]",
                          "b0[9]", "b0[10]", "b0[11]")) %>% 
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
  filter(Parameter %in% c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b0[5]", "b0[6]", "b0[7]", "b0[8]",
                          "b0[9]", "b0[10]", "b0[11]")) %>% 
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord3 + pWord4
ggsave("figures/supp/unimodal_consumption/validation_SharpeSchool_b0.png", width = 6.5, height = 6.5, dpi = 600)


#**** Group-level means ============================================================
# Plot posterior densities
unique(cs_df$Parameter)

p5 <- cs_df %>% 
  filter(Parameter %in% c("b1", "mu_b0", "sigma_b0", "mu_E", "Eh", "Th")) %>%
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
  filter(Parameter %in% c("b1", "mu_b0", "sigma_b0", "mu_E", "Eh", "Th")) %>%
  ggs_traceplot(.) +
  facet_wrap(~ Parameter, ncol = 2, scales = "free") +
  geom_line(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Iteration", y = "Value", color = "Chain #") +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  NULL

pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord5 + pWord6
ggsave("figures/supp/unimodal_consumption/validation_SharpeSchool.png", width = 6.5, height = 6.5, dpi = 600)


#**** Chain convergencve (Rhat) ====================================================
p7 <- cs_df %>% 
  ggs_Rhat(.) + 
  xlab("R_hat") +
  xlim(0.999, 1.01) +
  geom_point(size = 2) +
  NULL

pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord7
ggsave("figures/supp/unimodal_consumption/validation_rhat_SharpeSchool.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# b ~ dnorm(0.6, 1)
# mu_b0 ~ dnorm(1, 1)
# mu_E ~ dnorm(0.5, 4)
# sigma ~ dunif(0, 3) 
# sigma_b0 ~ dunif(0, 3) 
# sigma_b ~ dunif(0, 3) 
# sigma_E ~ dunif(0, 3)
# Eh ~ dnorm(2, 0.25)
# Th ~ dnorm(5, 0.25)

# Remember: distributions in JAGS have arguments mean and precision (inverse of variance)
# tau = 1/variance
# from sigma to tau 
# sigma = 1
# tau <- 1/sigma^2
# sigma^2 <- 1/tau
# sigma <- sqrt(1/tau)

# Define priors for plot
# set.seed(42)

mu_b0 <- rnorm(25000, mean = 1, sd = sqrt(1/(1)))
mu_E <- rnorm(25000, mean = 0.5, sd = sqrt(1/(4)))
Eh <- rnorm(25000, mean = 2, sd = sqrt(1/(0.25)))
Th <- rnorm(25000, mean = 5, sd = sqrt(1/(0.25)))

PR <- as.matrix(cbind(mu_b0, mu_E, Eh, Th))

# This is not a ggplot...
png(file = "/Users/maxlindmark/Desktop/R_STUDIO_PROJECTS/scaling/figures/supp/unimodal_consumption/validation_prior_post_SharpeSchool.png",
    units = "px", width = 1800, height = 1800, res = 300)

MCMCtrace(cs,
          params = c("mu_b0", "mu_E", "Eh", "Th"),
          ISB = FALSE,
          priors = PR,
          pdf = FALSE,
          Rhat = FALSE,
          n.eff = TRUE,
          type = "density",
          sz_txt = 1) # removes the trace plot

dev.off()


# Evaluate model fit & residuals =================================================
# https://rpubs.com/Niko/332320

# Extract generated data and data
cs_fit = coda.samples(jm, n.iter = n.iter, thin = thin,
                      variable.names = c("mean_y", "mean_y_sim", "p_mean"))

# Convert to data frames
cs_fit_df <- data.frame(as.matrix(cs_fit))

#-- Model fit
p_fit <- ggplot(cs_fit_df, aes(mean_y_sim)) + 
  coord_cartesian(expand = 0) +
  geom_histogram(bins = round(1 + 3.2*log(nrow(cs_fit_df)))) +
  geom_vline(xintercept = cs_fit_df$mean_y, color = "white", 
             linetype = 2, size = 0.4) +
  labs(x = "mean simulated data", y = "count") +
  theme_classic() +
  annotate("text", -Inf, Inf, label = round(mean(cs_fit_df$p_mean), digits = 3), 
           size = 3, hjust = -0.5, vjust = 1.3) +
  theme(text = element_text(size = 12), aspect.ratio = 1) +
  NULL


#-- Posterior predictive distributions
# https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/#posterior-predictive-analysis

# Extract posteriors for each data point for calculation of residuals
y_sim <- coda.samples(jm, variable.names = c("y_sim"), n.iter = n.iter, thin = thin)

# Tidy-up
df_y_sim <- ggs(y_sim)

pal <- brewer.pal(n = 3, name = "Dark2")

pp <- ggplot() +
  geom_density(data = df_y_sim, aes(value, fill = 'Posterior\nPredictive'), alpha = 0.6) +
  geom_density(data = con, aes(y_ct, fill = 'Observed'), alpha = 0.6) +
  scale_fill_manual(values = pal[c(3,2)]) +
  coord_cartesian(expand = 0) +
  theme_classic() +
  theme(text = element_text(size = 12), aspect.ratio = 1,
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.2, 0.95),
        legend.key.size = unit(0.3, "cm"))


#-- Residual vs fitted
df_y_sim <- df_y_sim %>%
  ungroup() %>%
  group_by(Parameter) %>%
  summarize(median = median(value)) %>% 
  rename("yhat" = "median") %>% 
  mutate(y = con$y_ct,
         resid = y - yhat)

p_resid <- ggplot(df_y_sim, aes(yhat, resid)) +
  geom_point(fill = "black", color = "white", shape = 21) + 
  theme_classic() +
  theme(text = element_text(size = 12)) 

(p_fit | pp) / p_resid + plot_annotation(tag_levels = 'A')

ggsave("figures/supp/unimodal_consumption/fit_pp_resid_SharpeSchool.png", width = 7, height = 7, dpi = 600)


# E. PLOT PREDICTIONS ==============================================================
#** Fits and data ==================================================================
# jags.samples - Nice for summaries and predictions
# Extract the prediction at each x including credible interval
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
                      temp = data$temp_pred)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
# Extend color palette
colourCount = length(unique(data$species))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
pal <- getPalette(colourCount)

min(pred_df$median)

p9 <- ggplot(pred_df, aes(temp, median)) +
  geom_ribbon(data = pred_df, aes(x = temp, ymin = lwr_95, ymax = upr_95),
             size = 0.3, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = temp, ymin = lwr_80, ymax = upr_80), 
              size = 0.3, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.3, alpha = 1, col = "black") +
  geom_point(data = con, aes(peak_temp_c_model_ct, y_ct, fill = species_ab),
             size = 1.2, alpha = 0.8, shape = 21, color = "white", stroke = 0.2) +
  scale_fill_manual(values = pal, name = "Species") +
  annotate("text", 5, 2.5, label = paste("n=", nrow(con), sep = ""), size = 2) +
  labs(x = "Rescaled temperature",
       y = "Rescaled consumption rate")
  
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 8),
                                       aspect.ratio = 1,
                                       legend.spacing.x = unit(0, 'cm'),
                                       legend.box.spacing = unit(0, 'cm'),
                                       legend.key.size = unit(0.1, 'cm'), 
                                       legend.title = element_text(size = 4),
                                       legend.position = c(0.15, 0.8),
                                       legend.text = element_text(size = 4, face = "italic"))
pWord9

ggsave("figures/SharpeSchool_scatter.png", width = 6, height = 6, dpi = 600, unit = "cm")

# Large version
p10 <- ggplot(pred_df, aes(temp, median)) +
  geom_ribbon(data = pred_df, aes(x = temp, ymin = lwr_95, ymax = upr_95),
              size = 0.3, alpha = 0.25, inherit.aes = FALSE, fill = "grey45") +
  geom_ribbon(data = pred_df, aes(x = temp, ymin = lwr_80, ymax = upr_80), 
              size = 0.3, alpha = 0.35, inherit.aes = FALSE, fill = "grey35") +
  geom_line(size = 0.8, alpha = 1, col = "black") +
  geom_point(data = con, aes(peak_temp_c_model_ct, y_ct, fill = species_ab),
             size = 2, alpha = 0.8, shape = 21, color = "white", stroke = 0.2) +
  scale_fill_manual(values = pal, name = "Species") +
  annotate("text", 5, 2.5, label = paste("n=", nrow(con), sep = ""), size = 3) +
  labs(x = "Rescaled temperature",
       y = "Rescaled consumption rate")

pWord10 <- p10 + theme_classic() + theme(text = element_text(size = 12),
                                         aspect.ratio = 1,
                                         legend.spacing.x = unit(0, 'cm'),
                                         legend.box.spacing = unit(0, 'cm'),
                                         legend.key.size = unit(0.1, 'cm'), 
                                         legend.title = element_text(size = 7),
                                         legend.position = c(0.15, 0.8),
                                         legend.text = element_text(size = 7, face = "italic"))
pWord10

ggsave("figures/SharpeSchool_scatter_v2.png", width = 10, height = 10, dpi = 600, unit = "cm")

