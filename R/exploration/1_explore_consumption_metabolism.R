#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.23: Max Lindmark
#
# - Explore consumption and metabolic data
# 
# A. Load libraries & read data
#
# B. Explore data
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
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Print package versions
# print(sessionInfo())
# other attached packages:
# magrittr_1.5       RColorBrewer_1.1-2 viridis_0.5.1      viridisLite_0.3.0  
# ggplot2_3.2.1      RCurl_1.95-4.12    bitops_1.0-6       readxl_1.3.1      
# tidylog_0.2.0      patchwork_0.0.1    tidyr_1.0.0        dplyr_0.8.3

# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

con <- read_excel("data/consumption_data.xlsx")
met <- read_excel("data/metabolism_data.xlsx")

# Which cols to make numeric?
glimpse(con)
glimpse(met)

con$type <- NA

cols = c(1,2,3,13,14,15,16,17,18)

con[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
met[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

con <- con %>% rename(y = consumption)
met <- met %>% rename(y = metabolic_rate)

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- data.frame(bind_rows(con, met))

head(dat)


# B. EXPLORE DATA ==================================================================
#** Normalize variables ============================================================
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate) + 
  coord_cartesian(expand = 0) +
  theme_classic() +
  NULL

dat$pref_temp_mid[is.na(dat$pref_temp_mid)] <- -9
dat$env_temp_mid[is.na(dat$env_temp_mid)] <- -9
t <- data.frame(filter(dat, pref_temp_mid == -9 & env_temp_mid == -9))

# Some species with no temperature-information
unique(t$common_name)

# For Stechlin cisco, well use the sympatric species Coregonus albula
dat$pref_temp_mid <- ifelse(dat$common_name == "Stechlin cisco", 
                            filter(dat, species == "Coregonus albula")$pref_temp_mid[1],
                            dat$pref_temp_mid)

filter(dat, common_name == "Stechlin cisco")$pref_temp_mid

# For Southern catfish, we'll assume it's the same as for the closely related species Amur catfish
# https://www.fishbase.se/summary/Silurus-asotus.html
dat$pref_temp_mid <- ifelse(dat$common_name == "Southern catfish", 
                            15,
                            dat$pref_temp_mid)

filter(dat, common_name == "Southern catfish")$pref_temp_mid

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

# Mean center predictor variable
dat$temp_norm_arr_ct <- dat$temp_norm_arr - mean(dat$temp_norm_arr)
dat$temp_norm_ct <- dat$temp_norm - mean(dat$temp_norm)

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass_g / dat$w_max_published_g
# ---- Stechlin cisco har larger size than max*

# Calculate log10 mass
dat$log_mass_norm <- log10(dat$mass_norm)

# Mean center mass
dat$log_mass_norm_ct <- dat$log_mass_norm - mean(dat$log_mass_norm)

# Create data with with species that have also temperature repliates within species (intra-specific analysis)
s_datc <- data.frame(
  dat %>% 
    dplyr::filter(rate == "consumption") %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(unique_t = n_distinct(temp_norm_arr_ct)) %>% 
    dplyr::filter(unique_t > 1)
)

str(s_datc)
ggplot(s_datc, aes(unique_t)) + geom_histogram()

s_datm <-  data.frame(
  dat %>% 
    dplyr::filter(rate == "metabolism") %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(unique_t = n_distinct(temp_norm_arr_ct)) %>% 
    dplyr::filter(unique_t > 1)
)

ggplot(s_datm, aes(unique_t)) + geom_histogram()

str(s_datm)


#** Plot general data ==============================================================
# Trophic level
ggplot(dat, aes(x = reorder(common_name, trophic_level), y = trophic_level)) +
  geom_point(stat = 'identity', size = 4) +
  scale_fill_manual(name = "trophic_level") + 
  theme_classic(base_size = 14) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("Trophic level") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# Temperature midpoint (Fishbase)
ggplot(dat, aes(x = reorder(common_name, env_temp_mid), 
                y = env_temp_mid)) +
  geom_point(stat = 'identity', size = 4) +
  scale_fill_manual(name = "env_temp_mid") + 
  theme_classic(base_size = 15) +
  guides(colour = T) +
  xlab("") + 
  ylab("Mid. Env. Temperature [C]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# Mid env. temperature (Fishbase) compared to experimental temperature range
pal <- rev(brewer.pal("Dark2", n = 5))
  
dat %>% 
  dplyr::group_by(common_name) %>% 
  ggplot(., aes(x = reorder(common_name, env_temp_mid), 
                y = env_temp_mid, color = "blue")) +
  geom_point(stat = 'identity', size = 1) +
  geom_errorbar(aes(reorder(common_name, env_temp_mid), 
                    ymin = env_temp_min, 
                    ymax = env_temp_max, color = "blue"), 
                size = 1,
                width = 0) +
  geom_point(aes(x = reorder(common_name, env_temp_mid), 
                 y = temp_c, color = "gray"), size = 1) +
  scale_color_manual(labels = c("Mid. Env. Temperature [C]", 
                                "Experimental\ntemperature"), 
                     values = c(pal)) +
  scale_fill_manual(name = "env_temp_mid") + 
  theme_classic(base_size = 12) +
  guides(color = guide_legend(title = "")) +
  xlab("") + 
  ylab("Temperature [C]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# ggsave("figures/supp/temperatures.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Max. published weight
ggplot(dat, aes(x = reorder(common_name, w_max_published_g), 
                y = log10(w_max_published_g))) +
  geom_point(stat = 'identity', size = 2) +
  scale_fill_manual(name = "w_max_published_g") + 
  theme_classic(base_size = 12) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("log10(max published weight) [g]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# ggsave("figures/supp/max_weight.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Phylogeny
dat %>% dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(order, fill = family)) +
  geom_bar() +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  #scale_fill_manual(values = pal) +  
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ rate) +
  NULL
# Doesn't seem like we can do much about phyologey here..

# ggsave("figures/supp/phylogeny.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Biogeography
dat %>% dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(biogeography, fill = biogeography)) +
  geom_bar() +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  #scale_fill_manual(values = mycolors) + 
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ rate) +
  guides(fill = FALSE) +
  NULL

# ggsave("figures/supp/biogeography.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Lifestyle
dat %>% 
  dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(habitat, fill  = lifestyle)) +
  geom_bar() +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  #scale_fill_manual(values = mycolors) + 
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ rate) +
  NULL

# ggsave("figures/supp/lifestyle.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")


#** Plot response variable =========================================================
str(dat)

# All data, by species
dat %>% 
  dplyr::filter(rate == "consumption") %>% 
  ggplot(., aes(temp_norm_ct, y)) + 
  geom_point(size = 2, alpha = 0.3) +
  theme_classic(base_size = 11) +
  facet_wrap(~species, scales = "free_y") +
  stat_smooth(se = F) +
  guides(color = F) +
  NULL

# Only intraspecific data, by species
s_datc %>% 
  ggplot(., aes(temp_norm_ct, y))+ 
  geom_point(size = 2, alpha = 0.3) +
  theme_classic(base_size = 11) +
  facet_wrap(~species, scales = "free_y") +
  stat_smooth(se = F) +
  guides(color = F) +
  NULL

s_datm %>% 
  ggplot(., aes(temp_norm_ct, y))+ 
  geom_point(size = 2, alpha = 0.3) +
  theme_classic(base_size = 11) +
  facet_wrap(~species, scales = "free_y") +
  stat_smooth(se = F) +
  guides(color = F) +
  NULL

# Add new column to see if measurements are below or within optimum of that size class and species

t <- c()

# Consumption
for(i in unique(s_datc$common_name)) {
  
  t <- filter(s_datc, common_name == i)
  t$round_mass_g <- round(t$mass_g, digits = 0)
  t$size_cl <- cut(t$round_mass_g, 10)
  title <- t$common_name[1]
  
  p <- t %>% 
    ggplot(., aes(temp_norm_ct, y, color = factor(size_cl)))+ 
    geom_jitter(size = 4, height = 0) +
    theme_classic(base_size = 11) + 
    scale_color_viridis(discrete = TRUE) +
    ggtitle(title)
    
  path <- paste("figures/supp/opt/con_", i, ".pdf", sep = "")
  ggsave(path, plot = p, scale = 1, width = 20, height = 20, units = "cm")
  
}  

# Metabolism
t <- c()

for(i in unique(s_datm$common_name)) {
  
  t <- filter(s_datm, common_name == i)
  t$round_mass_g <- round(t$mass_g, digits = 0)
  t$size_cl <- cut(t$round_mass_g, 10)
  title <- t$common_name[1]
  
  p <- t %>% 
    ggplot(., aes(temp_norm_ct, y, color = factor(size_cl)))+ 
    geom_jitter(size = 4, height = 0) +
    theme_classic(base_size = 11) + 
    scale_color_viridis(discrete = TRUE) +
    ggtitle(title)
  
  path <- paste("figures/supp/opt/met_", i, ".pdf", sep = "")
  ggsave(path, plot = p, scale = 1, width = 20, height = 20, units = "cm")
  
}  



# Group by normalized mass
# Consumption
p1 <- s_datc %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(temp_norm_ct, y_norm, color = log_mass_norm)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  #guides(color = FALSE) +
  scale_color_viridis() +
  labs(x = "Normalized temperature", y = "Normalized consumption") +
  stat_smooth(span = 1.0, se = F, size = 1, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 1.2, se = F, size = 1, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 0.8, se = F, size = 1, alpha = 0.4, geom = "line", color = pal[2]) +
  theme(aspect.ratio = 1) +
  NULL

# Metabolism
p2 <- s_datm %>% 
  dplyr::filter(unit == "mg O2/h") %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(temp_norm_ct, y_norm, color = log_mass_norm)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_viridis() +
  #guides(color = FALSE) +
  labs(x = "Normalized temperature", y = "Normalized metabolic rate") +
  stat_smooth(span = 1.0, se = F, size = 1, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 1.2, se = F, size = 1, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 0.8, se = F, size = 1, alpha = 0.4, geom = "line", color = pal[2]) +
  theme(aspect.ratio = 1) +
  NULL

p1 / p2

# ggsave("figures/supp/rates_temperature.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Now plot slope log(rate)~temp(arr), where the slope is the activation energy, 
# using sub-optimum temperatures! (<15 for consumption, based on previous plots). 
# Note I don't normalize here but plot the actual rates
# Consumption
p3 <- s_datc %>% 
  dplyr::filter(temp_norm_ct < 12) %>% 
  ggplot(., aes(temp_norm_arr_ct, log10(y), color = log_mass_norm_ct)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_viridis() +
  labs(x = "Inverse temperature [1/kT]", y = "log(consumption [g/day])") +
  stat_smooth(method = "lm", size = 1, alpha = 0.8, geom = "line", color = pal[2]) +
  NULL

# Metabolism
p4 <- s_datm %>% 
  ggplot(., aes(temp_norm_arr_ct, log10(y), color = log_mass_norm_ct)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_viridis() +
  labs(x = "Inverse temperature [1/kT]", y = "log(metabolic rate [mg O2/h])") +
  stat_smooth(method = "lm", size = 1, alpha = 0.8, geom = "line", color = pal[2]) +
  NULL

p3 / p4

# ggsave("figures/supp/rates_arr_temperature.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Now plot slope log(rate)~MASS, where the slope is the mass-exponent, 
# using sub-optimum temperatures! (<15 for consumption, based on previous plots). 
# Note I don't normalize here but plot the actual rates

# Consumption
p5 <- s_datc %>% 
  dplyr::filter(temp_norm_ct < 12) %>% 
  ggplot(., aes(log_mass_norm_ct, log10(y), color = temp_norm_ct)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_viridis() +
  labs(x = "Normalized mass", y = "log(consumption [g/day])", color = "Inverse temperature [1/kT]") +
  stat_smooth(method = "lm", size = 1, alpha = 0.8, geom = "line", color = pal[2]) +
  NULL

# Metabolism
p6 <- s_datm %>% 
  ggplot(., aes(log_mass_norm_ct, log10(y), color = temp_norm_ct)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_viridis() +
  labs(y = "log(metabolic rate [mg O2/h])", x = "Normalized mass", color = "Inverse temperature [1/kT]") +
  stat_smooth(method = "lm", size = 1, alpha = 0.8, geom = "line", color = pal[2]) +
  NULL

p5 / p6

# ggsave("figures/supp/rates_mass.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Seriously clustered data! color by species instead
# Consumption
p7 <- s_datc %>% 
  dplyr::filter(temp_norm_ct < 12) %>% 
  ggplot(., aes(log_mass_norm_ct, log10(y), color = species)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.6) +
  guides(color = FALSE) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Normalized mass", y = "log(consumption [g/day])", color = "Inverse temperature [1/kT]") +
  stat_smooth(method = "lm", size = 1, alpha = 0.8, geom = "line", color = pal[2]) +
  NULL

# Metabolism
p8 <- s_datm %>% 
  dplyr::filter(unit == "mg O2/h") %>% 
  ggplot(., aes(log_mass_norm_ct, log10(y), color = species)) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.6) +
  guides(color = FALSE) +
  scale_color_viridis(discrete = TRUE) +
  labs(y = "log(metabolic rate [mg O2/h])", x = "Normalized mass", color = "Inverse temperature [1/kT]") +
  stat_smooth(method = "lm", size = 1, alpha = 0.8, geom = "line", color = pal[2]) +
  NULL

p7 / p8

# ggsave("figures/supp/rates_species.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")

# Split the plot into discrete size classes. We normalize by ranking all normalized masses (which describes mass relative to max mass of that species) by splitting up the data in 8 classes, 8 being largest. We then normalize consumption within species, so that all consumption rates are relative to the maximum consumption rate within species, across all sizes and temperatures. We do this because the feeding rates differ across species (likely due to ecology and experimental setup).

# This is how ntile() works:
# k <- data.frame(value = head(unique(dat$mass_norm), 10) * 1000,
#                 species = rep(c("A", "B"), each = 5))
# 
# k %>% 
#   group_by(species) %>% 
#   dplyr::mutate(quartile = ntile(value, 4)) %>% 
#   arrange(species, quartile, value)

# Here I use only intra-specific data:
# Set number of ranks
nranks <- 5

# Group by normalized mass
# Consumption
p9 <- s_datc %>% 
  dplyr::mutate(quartile = ntile(mass_g, nranks)) %>%
  dplyr::group_by(quartile) %>% 
  dplyr::mutate(mean_mass_q = round(mean(mass_norm), digits = 2)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(temp_norm_ct, y_norm, color = factor(mean_mass_q))) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  labs(x = "Normalized temperature", y = "Normalized consumption") +
  scale_color_viridis(discrete = TRUE) +
  stat_smooth(span = 1.0, se = FALSE, size = 1, alpha = 0.8, geom = "line") +
  NULL

# Metabolism
p10 <- s_datm %>% 
  dplyr::mutate(quartile = ntile(mass_g, nranks)) %>%
  dplyr::group_by(quartile) %>% 
  dplyr::mutate(mean_mass_q = round(mean(mass_norm), digits = 2)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(temp_norm_ct, y_norm, color = factor(mean_mass_q))) + 
  theme_classic(base_size = 12) +
  geom_point(size = 1, alpha = 0.8) +
  labs(x = "Normalized temperature", y = "Normalized metabolism") +
  scale_color_viridis(discrete = TRUE) +
  stat_smooth(span = 1.0, se = F, size = 1, alpha = 0.8, geom = "line") +
  NULL

p9 / p10

# ggsave("figures/supp/rates_temp_discretemass.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")


#**** Summary ======================================================================
# Looks ok so far. I can show that there is a tendency in the data to have declining 
# optimum temperature for Cmax over ontogeny. I can also show that optimum curves 
# overall are more evident for Cmax than metabolism. This is ok to show as it is...
# ... But I will try and fit a polynomial like Garcia-Garcia to these data and see what 
# comes out. Then we can plot predicted lines for the average size in the size classes 
# I currently use for colors.
# Since everything is normalized... Can I compare Cmax and metabolism?
# Optimum "growth" is when the difference is maximized, so that should be OK as long 
# as I don't focus on the actual values! That gives me some data points I can compare
# with the T_opt figure.



