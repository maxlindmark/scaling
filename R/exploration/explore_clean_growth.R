#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore body growth rate data
# - For exploration of metadata - see 1_explore_Topt_growth.R
# 
# A. Load libraries & read data
#
# B. Explore growth data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES & READ DATA ====================================================
rm(list = ls())

# When doing a fresh start I need to check I'm in the right libpath:
# .libPaths() 
# .libPaths("C:/Program Files/R/R-3.5.0/library")

# Load libraries, install if needed
library(dplyr)
library(tidyr)
library(tidylog)
library(readxl)
library(RCurl)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(magrittr)
library(readr)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Print package versions
# print(sessionInfo())
# other attached packages:
# [1] patchwork_0.0.1    magrittr_1.5       RColorBrewer_1.1-2 viridis_0.5.1     
# [5] viridisLite_0.3.0  ggplot2_3.1.1      RCurl_1.95-4.12    bitops_1.0-6      
# [9] readxl_1.3.1       tidylog_0.1.0      tidyr_0.8.3        dplyr_0.8.1    

# ***  Will create a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))
dat <- read_excel("data/growth_data.xlsx")

glimpse(dat)

# Which cols to make numeric?
cols = c(1, 2, 3, 4, 5, 6, 16, 17, 18, 19, 10, 20, 21, 22)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)

colnames(dat)[1] <-  "G"


# B. CLEAN DATA ==================================================================
#** Predictors =====================================================================
dat$pref_temp_mid[is.na(dat$pref_temp_mid)] <- -9
dat$env_temp_mid[is.na(dat$env_temp_mid)] <- -9

# Some species with no temperature-information
# t <- data.frame(filter(dat, pref_temp_mid == -9 & env_temp_mid == -9))
# unique(t$common_name)

# Temperature-variable for analysis
dat$median_temp <- dat$env_temp_mid

# Take median of "preferred" temperature if environment temp is NA
# Replace NA with -9...
dat$median_temp <- ifelse(dat$median_temp == -9,
                          dat$pref_temp_mid,
                          dat$median_temp)

# Any NA's still?
dplyr::filter(dat, median_temp == -9)

# Convert median_temp to Arrhenius scale:
dat$median_temp_arr <- 1/((dat$median_temp + 273.15) * 8.617332e-05)

# Convert experimental to Arrhenius scale:
dat$temp_arr <- 1/((dat$temp_c + 273.15) * 8.617332e-05)

# Standardize temperatures to median-reference temperature on Arrhenius scale
dat$temp_norm_arr <- dat$temp_arr - dat$median_temp_arr

# Standardize temperatures to median-reference temperature on C scale
dat$temp_norm <- dat$temp_c - dat$median_temp

# Mean center predictor variables
dat$temp_norm_arr_ct <- dat$temp_norm_arr - mean(dat$temp_norm_arr)
dat$temp_norm_ct <- dat$temp_norm - mean(dat$temp_norm)

# We use geometric mean for size, and if this is not possible we'll go with mid point of size range. # Mass for analysis:
dat$mass <- dat$geom_mean_mass_g

# Replace NA with -9...
dat$mass[is.na(dat$mass)] <- -9

dat$mass <- ifelse(dat$mass == -9,
                   dat$size_group,
                   dat$mass)

# Check we only have real values:
dat$mass

# Calculate log mass
dat$log_mass <- log(dat$mass)

# Calculate mean-centered log mass
dat$log_mass_ct <- dat$log_mass - mean(dat$log_mass)

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass / dat$w_max_published_g

# Calculate log mass
dat$log_mass_norm <- log(dat$mass_norm)

# Mean center mass
dat$log_mass_norm_ct <- dat$log_mass_norm - mean(dat$log_mass_norm)

# Save data:
str(dat)
head(dat)

# dat %>%
#   select(G, mass, log_mass_ct, geom_mean_mass_g, temp_c, above_optimum, common_name, species,
#          median_temp, median_temp_arr, temp_arr, temp_norm_arr, temp_norm,
#          temp_norm_arr_ct, temp_norm, mass_norm, log_mass_norm, log_mass_norm_ct) %>%
#   write_csv(., "data/growth_analysis.csv", ";")


#** Growth rate ====================================================================
# Plot env. temperature (Fishbase) compared to experimental temperature range
pal <- brewer.pal("Dark2", n = 5)

dat %>% 
  select(temp_c, median_temp, species) %>% 
  gather(source, temp, 1:2) %>% 
  ggplot(., aes(x = reorder(species, temp), y = temp, color = source)) +
  geom_point(size = 2) +
  scale_color_manual(labels = c("Env. Temperature [C]", 
                                "Experimental\ntemperature"), 
                     values = pal[c(1,2)]) +
  theme_classic(base_size = 12) +
  guides(color = guide_legend(title = "")) +
  xlab("") + 
  ylab("Temperature [C]") + 
  coord_flip() +
  theme(legend.position = c(0.8, 0.2)) +
  NULL 
#ggsave("figures/supp/exp_env_temps_gro.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm", dpi = 300)

# Plot growth rate over temp within species, how many points beynd optimum?
ggplot(dat, aes(temp_norm_arr, G, color = factor(log10(mass_norm)))) + 
  geom_point(size = 2, alpha = 1) +
  geom_line() +
  facet_wrap(~common_name, scales = "free") +
  theme_bw(base_size = 15) +
  guides(color = FALSE) +
  scale_color_viridis(discrete = T, option = "cividis") +
  labs(x = "Relative temperature (environment", y = "Specific growth rate (%/day)") +
  #guides(color = F) +
  NULL

# Plot growth rate over temp within species, how many points beynd optimum?
# One species at the time
dat %>% filter(common_name == "Marbled flounder") %>% 
ggplot(., aes(temp_norm, G, color = factor(log10(mass_norm)))) + 
  geom_point(size = 2, alpha = 1) +
  geom_line() +
  #facet_wrap(~common_name, scales = "free") +
  theme_bw(base_size = 15) +
  guides(color = FALSE) +
  scale_color_viridis(discrete = T, option = "cividis") +
  labs(x = "Relative temperature (environment", y = "Specific growth rate (%/day)") +
  guides(color = F) +
  NULL
  
# G ~ temp, color = normalized mass
p1 <- ggplot(dat, aes(temp_norm, G, color = log10(mass_norm))) + 
  geom_point(size = 4, alpha = 1) +
  theme_bw(base_size = 15) +
  scale_color_viridis(discrete = FALSE, option = "cividis") +
  NULL

p2 <- ggplot(dat, aes(temp_norm, G, color = log10(mass))) + 
  geom_point(size = 4, alpha = 1) +
  theme_bw(base_size = 15) +
  scale_color_viridis(discrete = FALSE, option = "cividis") +
  NULL

p1 + p2

# ln(G) ~ ln(mass), color temp
ggplot(dat, aes(log(mass), log(G), color = temp_norm)) + 
  geom_point(size = 4, alpha = 1) +
  theme_bw(base_size = 15) +
  scale_color_viridis(discrete = FALSE, option = "cividis") +
  NULL

# Create discrete temperature-classes
dat <- dat %>% 
  dplyr::mutate(temp_cat = ntile(temp_norm, 20)) %>% 
  ungroup()

# Check it worked
ggplot(dat, aes(temp_cat, temp_norm, group = temp_cat)) + geom_boxplot()

# Number of data points within each temp-cat
ggplot(dat, aes(temp_cat)) + geom_bar()

# Check size-distribution in each temperature-class
ggplot(dat, aes(temp_cat, log10(mass_norm), group = temp_cat)) + geom_boxplot()

# Plot again and color by size-class
# ln(G) ~ ln(mass), color temp
ggplot(dat, aes(log(mass), log(G), color = factor(temp_cat), group = factor(temp_cat))) +   geom_point(size = 4, alpha = 0.7) +
  stat_smooth(method = "lm", se = FALSE, size = 2, alpha = 0.8) +
  theme_bw(base_size = 15) +
  scale_color_viridis(discrete = TRUE, option = "cividis") +
  NULL

# So, higher temperature is associated with higher slopes
summary(lm(log(dat$G) ~ log(dat$mass_norm) * dat$temp_norm))
summary(lm(log(dat$G) ~ log(dat$mass) * dat$temp_norm))
summary(lm(log(dat$G) ~ log(dat$mass) * dat$temp_cat))

