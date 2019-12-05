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
# C. Plot data
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
con_data = list(
  y = log(con$y), 
  n_obs = length(con$y), 
  species_n = con$species_n,
  mass = con$log_mass_norm_ct,
  temp = con$temp_norm_arr_ct
)


met_data = list(
  y = log(met$y), 
  n_obs = length(met$y), 
  species_n = met$species_n,
  mass = met$log_mass_norm_ct,
  temp = met$temp_norm_arr_ct
)

# C. PLOT DATA =====================================================================
# Plot response vs explanatory - color species
p1 <- ggplot(con, aes(log_mass_norm_ct, log(y), fill = species_ab, color = species_ab)) + 
  geom_point(size = 1.1, alpha = 0.7, shape = 21, stroke = 0.01) +
  #  stat_smooth(method = "lm", size = 2, alpha = 0.4, color = "black") +
  theme_classic(base_size = 11) +
  scale_fill_viridis(discrete = TRUE, option = "cividis") +
  scale_color_viridis(discrete = TRUE, option = "cividis") +
  guides(color = FALSE, fill = FALSE) +
  theme(aspect.ratio = 1) +
  labs(x = "ln(standardized mass)",
       y = "ln(maximun consumption rate)",
       color = "Standardized\nArrhenius\ntemperature") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

# Plot response vs explanatory - color species
p2 <- ggplot(met, aes(log_mass_norm_ct, log(y), fill = species_ab, color = species_ab)) + 
  geom_point(size = 1.1, alpha = 0.7, shape = 21, stroke = 0.01) +
  #  stat_smooth(method = "lm", size = 2, alpha = 0.4, color = "black") +
  theme_classic(base_size = 11) +
  scale_fill_viridis(discrete = TRUE, option = "cividis") +
  scale_color_viridis(discrete = TRUE, option = "cividis") +
  guides(color = FALSE, fill = FALSE) +
  theme(aspect.ratio = 1) +
  labs(x = "ln(standardized mass)",
       y = "ln(metabolic rate)",
       color = "Species") +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

# Plot response vs explanatory - color mass
p3 <- ggplot(con, aes(temp_norm_arr_ct, log(y), fill = log_mass_norm_ct, color = log_mass_norm_ct)) + 
  geom_point(size = 1.1, alpha = 0.7, shape = 21, stroke = 0.01) +
  #  stat_smooth(method = "lm", size = 2, alpha = 0.4, color = "black") +
  theme_classic(base_size = 11) +
  scale_fill_viridis(discrete = FALSE, option = "magma") +
  scale_color_viridis(discrete = FALSE, option = "magma") +
  theme(aspect.ratio = 1) +
  labs(x = "Standardized Arrhenius temperature",
       y = "ln(maximun consumption rate)",
       color = "ln(standardized mass)") +
  theme(legend.position = "bottom") +
  annotate("text", -Inf, Inf, label = "C", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

# Plot response vs explanatory - color mass
p4 <- ggplot(met, aes(temp_norm_arr_ct, log(y), fill = log_mass_norm_ct, color = log_mass_norm_ct)) + 
  geom_point(size = 1.2, alpha = 0.7, shape = 21, stroke = 0.01) +
  #  stat_smooth(method = "lm", size = 2, alpha = 0.4, color = "black") +
  theme_classic(base_size = 11) +
  scale_fill_viridis(discrete = FALSE, option = "magma") +
  scale_color_viridis(discrete = FALSE, option = "magma") +
  theme(aspect.ratio = 1) +
  labs(x = "Standardized Arrhenius temperature",
       y = "ln(metabolic rate)",
       color = "ln(standardized mass)") +
  theme(legend.position = "bottom") +
  annotate("text", -Inf, Inf, label = "D", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

(p1 + p2) / (p3 + p4)

#ggsave("figures/model_data.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)
