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

# sessionInfo()
# other attached packages:
# [1] tidylog_1.0.2      scales_1.1.1       MCMCvis_0.14.0     bayesplot_1.7.2    patchwork_1.0.1    viridis_0.5.1      viridisLite_0.3.0 
# [8] magrittr_2.0.1     readxl_1.3.1       RCurl_1.98-1.2     ggmcmc_1.4.1       ggplot2_3.3.2      tidyr_1.1.2        dplyr_1.0.2       
# [15] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-4  


# B. READ IN DATA ==================================================================
# Read in data file
con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Which species have data above optimum?
spec <- unique(filter(con, above_peak_temp == "Y"))$species

con <- con %>% filter(species %in% spec) %>% droplevels()

# Center temperature relative to mean in environment by species
con$temp_env_ct <- con$temp_c - con$median_temp

# Mean-center mass
con$mass_g_ct <- con$mass_g - mean(con$mass_g)

# Add mass-specific consumption
con$y_spec <- con$y / con$mass_g^0.62

# Express consumption rate as fraction of mean within species
con <- con %>%
  group_by(species_ab) %>%
  mutate(y_mean_species = mean(y_spec),
         y_ct = y_spec / y_mean_species) %>% 
  ungroup()

# Center consumption relative to the one where it is maximized and express consumption as % of max
# con <- con %>%
#   group_by(species) %>%
#   mutate(con_percent = y_spec/max(y_spec)) %>%
#   ungroup()

# Now summarize this data
# sum <- con %>%
#   group_by(species) %>%
#   filter(con_percent == 1) %>%
#   mutate(peak_temp = temp_c) %>% 
#   as.data.frame() %>% 
#   select(peak_temp, common_name)

# Now do a left_join
# con <- left_join(con, sum, by = "common_name") %>% as.data.frame()

# Standardize temperature
# con <- con %>% mutate(temp_ct = temp_c - peak_temp)

# Standardize to modeled peak temperature
# Loop through each species and fit a model with a quadratic temperature term
# Then filter the temperature at peak rate
tmp_dat <- data.frame()
datalist <- list()

for(i in unique(con$species_ab)) {
  
  tmp_dat <- con %>% filter(species_ab == i)
  
  tmp_dat$temp_c_sq <- tmp_dat$temp_c*tmp_dat$temp_c
  
  fit <- lm(y_ct ~ temp_c + temp_c_sq, data = tmp_dat)
  
  summary(fit)
  
  # Predict
  new_data <- data.frame(temp_c = seq(min(tmp_dat$temp_c), max(tmp_dat$temp_c),
                                      length.out = 50))
  
  new_data$temp_c_sq <- new_data$temp_c*new_data$temp_c
  
  new_data$preds <- predict(fit, newdata = new_data)
  
  topt <- new_data %>%
    filter(preds == max(preds)) %>%
    select(temp_c) %>% 
    rename(peak_temp_c_model = temp_c)
  
  topt$species_ab <- i
  
  ggplot(tmp_dat, aes(temp_c, y_ct)) +
    geom_point() +
    geom_line(data = new_data, aes(temp_c, preds), col = 'blue') +
    theme_bw(base_size = 12) +
    geom_vline(xintercept = topt$peak_temp_c_model) +
    labs(x = 'Temperature_centered (ºC)',
         y = 'Consumption rate') + 
    ggtitle(i)
  
  ggsave(paste("figures/supp/unimodal_consumption/single_species_SH_models/", i, ".png", sep = ""))
  
  datalist[[i]] <- topt
  
}

est <- dplyr::bind_rows(datalist)

est

# Add in modeled peak temperature
colnames(con)

con <- left_join(con, est)

head(con)

# Center temperature now with modeled temp
con$peak_temp_c_model_ct <- con$temp_c - con$peak_temp_c_model

# Filter species without a peak temperature
con <- con %>% filter(!species_ab == "P.yokohamae")

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
                inits = inits
)

update(jm, n.iter = n.iter)


# D. MODEL VALIDATION ==============================================================
# CODA - Nice for getting the raw posteriors
cs <- coda.samples(jm,
                   variable.names = c("mu_b0", "sigma_b0", "mu_E", "Eh", "Th",
                                      "E", "b0"), 
                   n.iter = n.iter, 
                   thin = thin)

summary(cs)

# 2. Quantiles for each variable:
#   
#          2.5%    25%    50%    75%  97.5%
# E[1]     0.5452 0.6096 0.6467 0.6861 0.7664
# E[2]     0.3319 0.3865 0.4170 0.4506 0.5212
# E[3]     0.4356 0.5639 0.6360 0.7116 0.8665
# E[4]     0.5499 0.6413 0.6912 0.7439 0.8517
# E[5]     0.3656 0.5001 0.5723 0.6494 0.7953
# E[6]     0.3632 0.4798 0.5406 0.6026 0.7243
# E[7]     0.4378 0.5393 0.5960 0.6549 0.7755
# E[8]     0.4115 0.5394 0.6109 0.6846 0.8470
# E[9]     0.2938 0.3460 0.3742 0.4035 0.4640
# E[10]    0.5237 0.6345 0.6964 0.7587 0.8892
# b0[1]    0.6476 0.7155 0.7519 0.7865 0.8547
# b0[2]    1.0107 1.0995 1.1448 1.1896 1.2772
# b0[3]    0.4048 0.5205 0.5809 0.6485 0.7825
# b0[4]    0.4544 0.5258 0.5647 0.6040 0.6781
# b0[5]    0.3978 0.4956 0.5541 0.6165 0.7409
# b0[6]    0.4611 0.5493 0.5995 0.6525 0.7651
# b0[7]    0.5248 0.6289 0.6863 0.7425 0.8529
# b0[8]    0.3780 0.4875 0.5439 0.6063 0.7252
# b0[9]    0.9177 0.9841 1.0187 1.0541 1.1207
# b0[10]   0.4161 0.5064 0.5568 0.6070 0.7077

#          2.5%    25%    50%    75%  97.5%
# Eh       2.1712 2.4675 2.6374 2.8278 3.2211
# Th       2.8158 3.6412 4.0288 4.3708 4.9798
# mu_E     0.4454 0.5308 0.5770 0.6254 0.7359
# mu_b0    0.5204 0.6463 0.7029 0.7611 0.8853
# sigma_b0 0.1498 0.2046 0.2463 0.3012 0.4666


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
  filter(Parameter %in% c("mu_b0", "sigma_b0", "mu_E", "Eh", "Th")) %>%
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
  filter(Parameter %in% c("mu_b0", "sigma_b0", "mu_E", "Eh", "Th")) %>%
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
  xlim(0.999, 1.003) +
  geom_point(size = 2) +
  NULL

pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 10),
                                       axis.text = element_text(size = 5))
pWord7
ggsave("figures/supp/unimodal_consumption/validation_rhat_SharpeSchool.png", width = 6.5, height = 6.5, dpi = 600)


#**** Prior vs posterior ===========================================================
# https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Priors from JAGS
# mu_b0 ~ dnorm(1, 1)
# mu_E ~ dnorm(0.5, 4)
# sigma_b0 ~ dunif(0, 3) 
# sigma_E ~ dunif(0, 3)
# sigma ~ dunif(0, 3) 
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
  annotate("text", 5, 2.8, label = paste("n=", nrow(con), sep = ""), size = 2) +
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
  annotate("text", 5, 2.8, label = paste("n=", nrow(con), sep = ""), size = 2) +
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

