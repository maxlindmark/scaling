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
# Read in your data file(s)
met <- read.csv("data/met_analysis.csv")
con <- read.csv("data/con_analysis.csv")

# Filter data points at below optimum temperatures
met <- met %>% filter(above_optimum == "N")
con <- con %>% filter(above_optimum == "N")

# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(con$species, 1, 1)
sp2 <- gsub( ".*\\s", "", con$species )
con$species_ab <- paste(sp1, sp2, sep = ".")

sp1 <- substring(met$species, 1, 1)
sp2 <- gsub( ".*\\s", "", met$species )
met$species_ab <- paste(sp1, sp2, sep = ".")

# Prepare data for JAGS
data = NULL # Clear any old data lists that might confuse things

# Rename species factor for JAGS (must be numbered 1:n)
con$species_n <- as.numeric(as.factor(con$species))
met$species_n <- as.numeric(as.factor(met$species))

# Data in list-format for JAGS
con <- con %>% arrange(species_n)

con_data = list(
  y = log(con$y), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  #mass = con$log_mass_norm_ct,
  #mass = con$log_mass_ct,
  #temp = con$temp_norm_arr_ct
  mass = con$log_mass_ct,
  temp = con$temp_arr - mean(con$temp_arr)
)

met <- met %>% arrange(species_n)

met_data = list(
  y = log(met$y), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  #mass = met$log_mass_norm_ct,
  #mass = met$log_mass_ct,
  #temp = met$temp_norm_arr_ct
  mass = met$log_mass_ct,
  temp = met$temp_arr - mean(met$temp_arr)
)

# met_data$species_n
# met %>% filter(species == "Gadus morhua") # Cod is species no. 16
# ggplot(met, aes(species_n, species)) + geom_point()
# sort(unique(met_data$species_n))


# C. PLOT SPECIES VARYING COEFFICIENTS =============================================
# Refit chosen models from the model selection part

# Maximum consumption rate
con_model = "R/analysis/log-linear model/models/m5.txt"

jm_con = jags.model(con_model,
                    data = con_data, 
                    n.adapt = 5000, 
                    n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm_con, n.iter = burn.in) 


# Metabolic rate
met_model = "R/analysis/log-linear model/models/m2.txt"

jm_met = jags.model(met_model,
                    data = met_data, 
                    n.adapt = 5000, 
                    n.chains = 3)

burn.in = 10000 # Length of burn-in

update(jm_met, n.iter = burn.in) 

# Sample from the posteriors =======================================================
samples = 10000 # How many samples to take from the posterior
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

# Check overlap with zero
js = jags.samples(jm_con, 
                  variable.names = c("mu_b1"), 
                  n.iter = samples, 
                  thin = n.thin)

ecdf(js$mu_b1)(0.75) 
#0.9295

1- ecdf(js$b3)(0) 


# Metabolic rate
cs_met <- coda.samples(jm_met,
                       variable.names = c("b0", "b1", "b2", "b3",
                                          "mu_b0", "mu_b1", "mu_b2", 
                                          "sigma_b0", "sigma_b1", "sigma_b2",
                                          "sigma"),
                       n.iter = samples, 
                       thin = n.thin)

summary(cs_met)

# Check overlap with zero
js = jags.samples(jm_met, 
                  variable.names = c(
                                     #"mu_b1", 
                                     "b3"), 
                  n.iter = samples, 
                  thin = n.thin)

1-ecdf(js$mu_b1)(0.75) 
#0.9295

1-ecdf(js$b3)(0)


#** Plot posteriors ================================================================
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
con_b$Rate <- "Maximum Consumption"
con_b$Parameter_mte <- "Mass-exponent"
con_b$pred <- filter(con_df, Parameter == "mu_b1")$quantiles.50.
con_b$pred_sd <- filter(std_con, Parameter == "mu_b1")$statistics.SD

#** Activation energy
con_e <- con_df %>% filter(Parameter_sub == "b2")
con_e$Species <- unique(con$species_ab)
#con_e$Species <- unique(con_data$species)
con_e$Rate <- "Maximum Consumption"
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
met_b$Parameter_mte <- "Mass-exponent"
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
#pal <- brewer.pal("Dark2", n = 5)
#pal <- viridis(option = "magma", n = 10)[c(2, 6)]
pal <- brewer.pal("Dark2", n = 5)[c(1,3)]

# Create data frame for rectangles
df_std <- df[!duplicated(df$pred_sd), ]
df_std$ymax <- df_std$pred + 2*df_std$pred_sd
df_std$ymin <- df_std$pred - 2*df_std$pred_sd

# Plot all species varying estimates and global mean
df %>% 
  filter(Parameter_mte %in% c("Activation energy", "Mass-exponent")) %>% 
  ggplot(., aes(Species, quantiles.50., color = Rate, shape = Rate)) +
  facet_grid(~ Parameter_mte, scales = "free") +
  scale_color_manual(values = pal[1:2]) +
  scale_fill_manual(values = pal[1:2]) +
  scale_shape_manual(values = c(21, 24)) +
  geom_hline(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass-exponent")), 
             aes(yintercept = pred, color = Rate),
             size = 0.6, alpha = 1, linetype = "dashed") +
  geom_rect(data = filter(df_std, Parameter_mte %in% c("Activation energy", "Mass-exponent")), 
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
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  theme(aspect.ratio = 2/1,
        legend.position = "bottom", 
        legend.title = element_blank()) +
  NULL

#ggsave("figures/species_b_ea.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm", dpi = 300)


#** Plot global-predictions ========================================================
color_scheme_set("gray")

p1 <- cs_met %>% 
  mcmc_dens(pars = "mu_b1") +
  theme_classic(base_size = 11) + 
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Mass-exponent" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass-exponent") +
  ggtitle("") +
  xlim(0.48, 0.9) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Mass-exponent" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL

p2 <- cs_met %>% 
  mcmc_dens(pars = "mu_b2") +
  theme_classic(base_size = 11) + 
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Activation energy" & Rate == "Metabolic rate")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Activation energy") +
  ggtitle("Metabolic rate") +
  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Activation energy" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  NULL

p3 <- cs_met %>% 
  mcmc_dens(pars = "b3") +
  theme_classic(base_size = 11) + 
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "M*T interaction" & Rate == "Metabolic rate")$pred, 3)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "m*t interaction") +
  ggtitle("") +
#  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "M*T interaction" & Rate == "Metabolic rate")$pred, 
             linetype = "dashed", color = "white") +
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "red") +
  NULL

p4 <- cs_con %>% 
  mcmc_dens(pars = "mu_b1") +
  theme_classic(base_size = 11) + 
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Mass-exponent" & Rate == "Maximum Consumption")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass-exponent") +
  ggtitle("") +
  xlim(0.48, 0.9) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Mass-exponent" & Rate == "Maximum Consumption")$pred, 
             linetype = "dashed", color = "white") +
  NULL

p5 <- cs_con %>% 
  mcmc_dens(pars = "mu_b2") +
  theme_classic(base_size = 11) + 
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", -Inf, Inf, label = "C", size = 4, 
  #          fontface = "bold", hjust = -0.5, vjust = 1.3) +
  annotate("text", -Inf, Inf, label = round(filter(df, Parameter_mte == "Activation energy" & Rate == "Maximum Consumption")$pred, 2)[1], 
           size = 3, hjust = -0.5, vjust = 1.3) +
  labs(x = "Activation energy") +
  ggtitle("Maximum consumption rate") +
  xlim(-0.95, -0.4) +
  geom_vline(xintercept = filter(df, Parameter_mte == "Activation energy" & Rate == "Maximum Consumption")$pred, 
             linetype = "dashed", color = "white") +
  NULL

p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 3)

#ggsave("figures/supp/posterior_mte_parameters.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)


# How much of the interaction coefficient overlaps 0?
# Metabolic rate
js = jags.samples(jm_met, 
                  variable.names = c("b3"), 
                  n.iter = samples, 
                  thin = n.thin)

1-ecdf(js$b3)(0) # We are % certain the slope is smaller than 0
# [1] 0.9973333

# How big is a c of 0.014 on celcius scale?
b <- 0.75
cc <- 0.017

test <- con

test$b_pred <- b + cc * test$temp_norm_arr_ct

par(mfrow = c(2,1))
plot(test$b_pred ~ test$temp_norm_arr_ct)
plot(test$b_pred ~ test$temp_norm, col = "red")

summary(lm(test$b_pred ~ test$temp_norm))$coefficients[2]

# Not that big of a difference in the exponents.
