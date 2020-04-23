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
# C. Plot parameter estimates
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
library(tidybayes)

# other attached packages:
# [1] bayesplot_1.7.1    patchwork_0.0.1    viridis_0.5.1      viridisLite_0.3.0  magrittr_1.5       readxl_1.3.1      
# [7] RCurl_1.95-4.12    bitops_1.0-6       ggmcmc_1.3         ggplot2_3.2.1      tidyr_1.0.0        dplyr_0.8.3       
# [13] RColorBrewer_1.1-2 rjags_4-10         coda_0.19-3    


# B. READ IN DATA ==================================================================
# Read in your data file
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Filter data points at below optimum temperatures
met <- met %>% filter(above_optimum == "N")
con <- con %>% filter(above_optimum == "N")

# Rename species factor for JAGS (must be numbered 1:n)
met$species_n <- as.numeric(as.factor(met$species))
con$species_n <- as.numeric(as.factor(con$species))

# Mean-center predictor variables
met$log_mass_ct <- met$log_mass - mean(met$log_mass)
met$temp_arr_ct <- met$temp_arr - mean(met$temp_arr)
con$log_mass_ct <- con$log_mass - mean(con$log_mass)
con$temp_arr_ct <- con$temp_arr - mean(con$temp_arr)

# Use mass-specific values
met$y_spec <- met$y / met$mass_g
con$y_spec <- con$y / con$mass_g

# Prepare data for JAGS
met_data = NULL # Clear any old data lists that might confuse things
con_data = NULL

# Data in list-format for JAGS
met_data = list(
  y = log(met$y_spec), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_ct,
  temp = met$temp_arr_ct
)

# Data in list-format for JAGS
con_data = list(
  y = log(con$y_spec), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_ct,
  temp = con$temp_arr_ct
)


# C. PLOT PARAMETER ESTIMATES ======================================================
# Refit chosen models from the model selection part

# Maximum consumption rate
con_model = "R/analysis/JAGS_models/log_linear/m5.txt"

# Manually set initial values, because otherwise all the chains get the same
# NOTE I don't do it for all parameters...
inits_con = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))

jm_con = jags.model(con_model,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3,
                    inits = inits_con)

burn.in = 15000 # Length of burn-in

update(jm_con, n.iter = burn.in) 


# Metabolic rate
met_model = "R/analysis/JAGS_models/log_linear/m2.txt"

# Manually set initial values, because otherwise all the chains get the same
inits_met = list(
  list(
    mu_b0 = 0.1,
    mu_b1 = 0.1,
    mu_b2 = 0.1,
    b3 = 1,
    sigma = 0.1,
    sigma_b0 = 0.1,
    sigma_b1 = 0.1,
    sigma_b2 = 0.1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 1,
    mu_b1 = 1,
    mu_b2 = 1,
    b3 = 1,
    sigma = 1,
    sigma_b0 = 1,
    sigma_b1 = 1,
    sigma_b2 = 1,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ),
  list(
    mu_b0 = 2,
    mu_b1 = 2,
    mu_b2 = 2,
    b3 = 2,
    sigma = 2,
    sigma_b0 = 2,
    sigma_b1 = 2,
    sigma_b2 = 2,
    .RNG.name = "base::Super-Duper", .RNG.seed = 2
  ))


jm_met = jags.model(met_model,
                    data = met_data, 
                    n.adapt = 5000, 
                    n.chains = 3,
                    inits = inits_met)

burn.in = 15000 # Length of burn-in

update(jm_met, n.iter = burn.in) 

# Sample from the posteriors =======================================================
samples = 15000 # How many samples to take from the posterior
n.thin = 5 # Thinning?

# CODA - Nice for getting the raw posteriors
# Maximum consumption rate
cs_con <- coda.samples(jm_con,
                       variable.names = c("b0", "b1", "b2", 
                                          "mu_b0", "mu_b1", "mu_b2", 
                                          "sigma_b0", "sigma_b1", "sigma_b2",
                                          "sigma"), 
                       n.iter = samples, 
                       thin = n.thin)

summary(cs_con) # Get the mean estimate and SE and 95% CIs

# Metabolic rate - NOTE I'm also sampling b3 here! That parameter is not in consumption
cs_met <- coda.samples(jm_met,
                       variable.names = c("b0", "b1", "b2", "b3",
                                          "mu_b0", "mu_b1", "mu_b2", 
                                          "sigma_b0", "sigma_b1", "sigma_b2",
                                          "sigma"),
                       n.iter = samples, 
                       thin = n.thin)

summary(cs_met)

# Convert to ggplottable data frame


#** Plot species-predictions =======================================================
# Maxiumum consumption
con_df <- data.frame(summary(cs_con)[2])
con_df$Parameter <- rownames(con_df)
con_df$Parameter_sub <- factor(substring(con_df$Parameter, 1, 2))

std_con <- data.frame(summary(cs_con)[1])
std_con$Parameter <- rownames(std_con)

#** Mass exponent
con_b <- con_df %>% filter(Parameter_sub == "b1")
con_b$Species <- unique(con$species_ab)
#con_b$Species <- unique(con_data$species)
con_b$Rate <- "Maximum consumption rate"
con_b$Parameter_mte <- "Mass exponent"
con_b$pred <- filter(con_df, Parameter == "mu_b1")$quantiles.50.
con_b$pred_sd <- filter(std_con, Parameter == "mu_b1")$statistics.SD

#** Activation energy
con_e <- con_df %>% filter(Parameter_sub == "b2")
con_e$Species <- unique(con$species_ab)
#con_e$Species <- unique(con_data$species)
con_e$Rate <- "Maximum consumption rate"
con_e$Parameter_mte <- "Activation energy"
con_e$pred <- filter(con_df, Parameter == "mu_b2")$quantiles.50.
con_e$pred_sd <- filter(std_con, Parameter == "mu_b2")$statistics.SD

# Metabolism
met_df <- data.frame(summary(cs_met)[2])
met_df$Parameter <- rownames(met_df)
met_df$Parameter_sub <- factor(substring(met_df$Parameter, 1, 2))

std_met <- data.frame(summary(cs_met)[1])
std_met$Parameter <- rownames(std_met)

#** Mass exponent
met_b <- met_df %>% filter(Parameter_sub == "b1")
met_b$Species <- unique(met$species_ab)
#met_b$Species <- unique(met_data$species)
met_b$Rate <- "Metabolic rate"
met_b$Parameter_mte <- "Mass exponent"
met_b$pred <- filter(met_df, Parameter == "mu_b1")$quantiles.50.
met_b$pred_sd <- filter(std_met, Parameter == "mu_b1")$statistics.SD

#** Activation energy
met_e <- met_df %>% filter(Parameter_sub == "b2")
met_e$Species <- unique(met$species_ab)
#met_e$Species <- unique(met_data$species)
met_e$Rate <- "Metabolic rate"
met_e$Parameter_mte <- "Activation energy"
met_e$pred <- filter(met_df, Parameter == "mu_b2")$quantiles.50.
met_e$pred_sd <- filter(std_met, Parameter == "mu_b2")$statistics.SD

#** M*T interaction
met_c <- met_df %>% filter(Parameter == "b3")
met_c$Species <- NA
met_c$Rate <- "Metabolic rate"
met_c$Parameter_mte <- "M*T interaction"
met_c$pred <- filter(met_df, Parameter == "b3")$quantiles.50.
met_c$pred_sd <- NA

# Merge data frames
df <- rbind(con_b, con_e, met_b, met_e, met_c)

# Define color palettes
pal <- brewer.pal("Dark2", n = 5)[c(1,3)]

# Convert temperature ceofficient to activation energy by multiplying with -1
df$quantiles.2.5. <- ifelse(df$Parameter_mte == "Activation energy",
                            df$quantiles.2.5. * -1,
                            df$quantiles.2.5.)

df$quantiles.25. <- ifelse(df$Parameter_mte == "Activation energy",
                           df$quantiles.25. * -1,
                           df$quantiles.25.)

df$quantiles.50. <- ifelse(df$Parameter_mte == "Activation energy",
                           df$quantiles.50. * -1,
                           df$quantiles.50.)

df$quantiles.75. <- ifelse(df$Parameter_mte == "Activation energy",
                           df$quantiles.75. * -1,
                           df$quantiles.75.)

df$quantiles.97.5. <- ifelse(df$Parameter_mte == "Activation energy",
                             df$quantiles.97.5. * -1,
                             df$quantiles.97.5.)

df$pred <- ifelse(df$Parameter_mte == "Activation energy",
                  df$pred * -1,
                  df$pred)

df$pred_sd <- ifelse(df$Parameter_mte == "Activation energy",
                     df$pred_sd * -1,
                     df$pred_sd)

# Create data frame for rectangles
df_std <- df[!duplicated(df$pred_sd), ]
df_std$ymax <- df_std$pred + 2*df_std$pred_sd
df_std$ymin <- df_std$pred - 2*df_std$pred_sd

# Plot all species varying estimates and global mean
p1 <- df %>% 
  filter(Parameter_mte %in% c("Activation energy", "Mass exponent")) %>% 
  ggplot(., aes(Species, quantiles.50., color = Rate, shape = Rate)) +
  facet_grid(Rate ~ Parameter_mte, scales = "free") +
  guides(color = FALSE, fill = FALSE, shape = FALSE) +
  scale_color_manual(values = pal[1:2]) +
  scale_fill_manual(values = pal[1:2]) +
  scale_shape_manual(values = c(21, 24)) +
  geom_hline(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass exponent")), 
             aes(yintercept = pred, color = Rate),
             size = 0.6, alpha = 1, linetype = "dashed") +
  geom_rect(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass exponent")), 
            inherit.aes = FALSE, aes(ymin = ymin, ymax = ymax, fill = Rate), xmin = 0, xmax = 50, 
            alpha = 0.2) +
  coord_flip() +
  geom_errorbar(aes(Species, quantiles.50., color = Rate, 
                    ymin = quantiles.2.5., ymax = quantiles.97.5.),
                size = 1, width = 0, alpha = 0.4) +
  geom_errorbar(aes(Species, quantiles.50., color = Rate, 
                    ymin = quantiles.25., ymax = quantiles.75.), 
                size = 1.5, width = 0, alpha = 0.7) +
  geom_point(size = 1.5, fill = "white") +
  labs(x = "Species", y = "Prediction") + 
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  theme(aspect.ratio = 2/1,
        legend.position = "bottom", 
        legend.title = element_blank()) +
  NULL

pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 7))
ggsave("figures/species_b_ea.png", width = 6.5, height = 6.5, dpi = 600)

#ggsave("figures/species_b_ea.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


#** Plot global-predictions ========================================================
color_scheme_set("gray")

p2 <- cs_met %>% 
  mcmc_dens(pars = "mu_b1") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Mass exponent" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass exponent") +
  ggtitle("") +
  xlim(-0.6, -0.1) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Mass exponent" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))

p3 <- cs_met %>% 
  mcmc_dens(pars = "mu_b2") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Activation energy" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Temperture coefficient") +
  ggtitle("Metabolic rate") +
  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Activation energy" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL
pWord3 <- p3 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))

p4 <- cs_met %>% 
  mcmc_dens(pars = "b3") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "M*T interaction" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "M*T interaction") +
  ggtitle("") +
  geom_vline(xintercept = filter(df, Parameter_mte == "M*T interaction" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "red") +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))

p5 <- cs_con %>% 
  mcmc_dens(pars = "mu_b1") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Mass exponent" & Rate == "Maximum consumption rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass exponent") +
  ggtitle("") +
  xlim(-0.6, -0.1) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Mass exponent" & Rate == "Maximum consumption rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))

p6 <- cs_con %>% 
  mcmc_dens(pars = "mu_b2") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Activation energy" & Rate == "Maximum consumption rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Temperture coefficient") +
  ggtitle("Maximum consumption rate") +
  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Activation energy" & Rate == "Maximum consumption rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL
pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text = element_text(size = 12))

pWord2 + pWord3 + pWord4 + pWord5 + pWord6 + plot_layout(ncol = 3)

ggsave("figures/supp/log_linear_model/met_con/species_b_ea.png", width = 6.5, height = 6.5, dpi = 600)


#** Do calculations on the posterior... ============================================
# How much of the interaction coefficient overlaps 0?
# Metabolic rate
js = jags.samples(jm_met, 
                  variable.names = c("b3"), 
                  n.iter = samples, 
                  thin = n.thin)

1-ecdf(js$b3)(0) # We are % certain the slope is smaller than 0
# [1] 0.9973333

# How big is a c of 0.015 on celcius scale?
b <- 0.75
cc <- 0.015

test <- con

test$b_pred <- b + cc * test$temp_norm_arr_ct

par(mfrow = c(2,1))
plot(test$b_pred ~ test$temp_norm_arr_ct)
plot(test$b_pred ~ test$temp_norm, col = "red")

summary(lm(test$b_pred ~ test$temp_norm))$coefficients[2]

# Not that big of a difference in the exponents.
