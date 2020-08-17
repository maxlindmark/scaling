#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore growth data
# 
# A. Load libraries & read data
#
# B. Explore & clean data
#
# C. Save data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES & READ DATA ====================================================
rm(list = ls())

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

dat <- read_excel("data/growth_data.xlsx")

glimpse(dat)

# Which cols to make numeric?
cols = c(19, 20, 21, 22)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)

dat$y <- dat$`growth_rate_%/day`

glimpse(dat)


# B. EXPLORE & CLEAN DATA ==========================================================
# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Create a single reference temperature for analysis. This is midpoint of environment (mainly),
# but sometimes midpoint of preferred (both from fishbase), and in two cases other literature
dat$pref_temp_mid[is.na(dat$pref_temp_mid)] <- -9
dat$env_temp_mid[is.na(dat$env_temp_mid)] <- -9

# New reference temperature (either mid of preference of environment temperature)
dat$median_temp <- dat$env_temp_mid

# Take median of "preferred" temperature if environment temp is NA
# Replace NA with -9...
dat$median_temp <- ifelse(dat$median_temp == -9,
                          dat$pref_temp_mid,
                          dat$median_temp)

# Bring back NA
dat$env_temp_mid <- ifelse(dat$env_temp_mid == -9,
                           NA,
                           dat$env_temp_mid)

dat$pref_temp_mid <- ifelse(dat$pref_temp_mid == -9,
                            NA,
                            dat$pref_temp_mid)

# Any NA's still?
dplyr::filter(dat, median_temp == -9)

# Inspect temperatures
ggplot(dat, aes(median_temp, fill = common_name)) +
  geom_histogram() + 
  coord_cartesian(expand = 0) +
  theme_classic() +
  NULL

# Convert experimental to Arrhenius scale:
dat$temp_arr <- 1/((dat$temp_c + 273.15) * 8.617332e-05)

# Standardize temperatures to median-reference temperature on C scale
dat$temp_norm <- dat$temp_c - dat$median_temp

# Use size_group if no geometric mean mass
dat$mass_g <- dat$geom_mean_mass_g

dat$mass_g[is.na(dat$mass_g)] <- -9

dat$mass_g <- ifelse(dat$mass_g == -9,
                     dat$size_group,
                     dat$mass_g)

# Calculate log mass
dat$log_mass <- log(dat$mass_g)

# Normalize mass with respect to max mass
dat$mass_norm <- dat$mass_g / dat$w_maturation_g

# Calculate log normalized mass
dat$log_mass_norm <- log(dat$mass_norm)


#** Plot general data ==============================================================
# Plot env. temperature (Fishbase) compared to experimental temperature range

# Trophic level
p1 <- ggplot(dat, aes(x = reorder(species, trophic_level), y = trophic_level)) +
  geom_point(stat = 'identity', size = 2) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("Trophic level") + 
  coord_flip() +
  NULL 
pWord1 <- p1 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text.y = element_text(size = 8, face = "italic"))

# Max. published weight
p2 <- ggplot(dat, aes(x = reorder(species, w_maturation_g), y = w_maturation_g)) +
  geom_point(stat = 'identity', size = 2) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("Mass at maturation [g]") + 
  scale_y_continuous(trans = 'log10') +
  coord_flip() +
  NULL 
pWord2 <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text.y = element_text(size = 8, face = "italic"))
pWord1 / pWord2
ggsave("figures/supp/data/growth_tl_mass.png", width = 6.5, height = 6.5, dpi = 600)


# Mid env. temperature (Fishbase) compared to experimental temperature range
pal <- brewer.pal("Dark2", n = 5)

p3 <- ggplot(dat) +
  coord_flip() +
  geom_point(aes(x = reorder(species, median_temp), 
                 y = median_temp, color = "Environment (median)"), size = 1.5, alpha = 0.6) +
  geom_point(aes(x = reorder(species, env_temp_mid), 
                 y = temp_c, color = "Experiment"), size = 1.5, alpha = 0.6) +
  geom_point(aes(x = reorder(species, env_temp_max), 
                 y = env_temp_max, color = "Environment (max)"), size = 1.5, alpha = 0.6) +
  geom_point(aes(x = reorder(species, env_temp_min), 
                 y = env_temp_min, color = "Environment (min)"), size = 1.5, alpha = 0.6) +
  scale_color_manual(values = rev(pal), name = "Temperature") +
  xlab("") + 
  guides(shape = F) +
  ylab(expression(paste("Temperature [", degree*C, "]"))) + 
  NULL 
pWord <- p3 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text = element_text(size = 8),
                                      axis.text.y = element_text(face = "italic"))
ggsave("figures/supp/data/growth_temperatures.png", width = 6.5, height = 6.5, dpi = 600)


# Phylogeny
p4 <- dat %>% dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(order, fill = family)) +
  geom_bar() +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  NULL
pWord4 <- p4 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
                                       legend.text = element_text(size = 6),
                                       legend.key.size = unit(0.4, "cm"))

# Lifestyle
p5 <- dat %>% 
  dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(habitat, fill  = lifestyle)) +
  geom_bar() +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  NULL
pWord5 <- p5 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text.x = element_text(angle = 30, hjust = 1),
                                       legend.text = element_text(size = 6),
                                       legend.key.size = unit(0.4, "cm"))
(pWord4 / pWord5) 
ggsave("figures/supp/data/growth_life_phylo.png", width = 6.5, height = 6.5, dpi = 600)


# Biogeography
p6 <- dat %>% dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(biogeography, fill = biogeography)) +
  geom_bar() +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  guides(fill = FALSE) +
  NULL
pWord6 <- p6 + theme_classic() + theme(text = element_text(size = 12),
                                       axis.text.x = element_text(angle = 0, hjust = 1))

# Test which sizes I use
p7 <- ggplot(dat, aes(mass_norm, fill = species)) + 
  geom_histogram() + 
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  coord_cartesian(expand = 0) + 
  labs(x = "Mass/Maturation mass") +
  guides(fill = FALSE) +
  NULL
pWord7 <- p7 + theme_classic() + theme(text = element_text(size = 12),
                                       legend.text = element_text(size = 8))
pWord6 / pWord7
ggsave("figures/supp/data/growth_bio_mass.png", width = 6.5, height = 6.5, dpi = 600)


#** Plot response variable =========================================================
# Plot growth rate over temp within species, how many points beynd optimum?
p8 <- ggplot(dat, aes(temp_c, y, 
                      group = factor(log10(mass_norm)),
                      color = log10(mass_norm))) + 
  geom_point(size = 2, alpha = 1) +
  geom_line() +
  facet_wrap(~common_name, scales = "free") +
  scale_color_viridis(option = "magma") +
  labs(x = expression(paste("Temperature [", degree*C, "]")), 
       y = "Specific growth rate (%/day)") +
  NULL
pWord <- p8 + theme_classic() + theme(text = element_text(size = 12),
                                      strip.text.x = element_text(size = 6.5))
ggsave("figures/supp/data/species_plots/growth/growth.png", width = 6.5, height = 6.5, dpi = 600)


#**** Plot all data combined =======================================================
dat$sample_size <- paste("n=", nrow(dat), sep = "")

# As a function of temperature
p9 <- dat %>% filter(y > 0) %>% 
  ggplot(., aes(x = temp_arr, y = y, fill = log10(mass_g))) + 
  geom_point(size = 3, alpha = 0.8, color = "white", shape = 21) +
  theme_classic(base_size = 11) +
  scale_y_log10() +
  labs(y = "specific growth rate (%/day)",
       x = "Arrhenius temperature [1/kT]") +
  scale_fill_viridis(option = "magma", name = "log10(mass)\n[g]") +
  geom_text(aes(x = Inf, y = Inf, label = sample_size, hjust = 1.05, vjust = 1.5), 
            color = "black") + 
  NULL
pWord9 <- p9 + theme_classic() + theme(text = element_text(size = 12))

# As a function of mass
p10 <- dat %>% filter(y > 0) %>% 
  ggplot(., aes(x = mass_g, y = y, fill = temp_arr)) + 
  geom_point(size = 3, alpha = 0.8, color = "white", shape = 21) +
  theme_classic(base_size = 11) +
  labs(y = "specific growth rate (%/day)",
       x = "mass [g]") +
  scale_fill_viridis(option = "magma", name = "Arrhenius\ntemperature\n[1/kT]") +
  scale_y_log10() +
  scale_x_log10() +
  geom_text(aes(x = Inf, y = Inf, label = sample_size, hjust = 1.05, vjust = 1.5), 
            color = "black") +
  NULL
pWord10 <- p10 + theme_classic() + theme(text = element_text(size = 12))
pWord9 / pWord10
ggsave("figures/supp/data/growth_rate_temp_mass.png", width = 6.5, height = 6.5, dpi = 600)


# C. SAVE DATA =====================================================================
glimpse(dat)
dat %>%
  select(y, `growth_rate_%/day`, geom_mean_mass_g, size_group, mass_g, log_mass, mass_norm, log_mass_norm, 
         temp_c, temp_arr, median_temp, above_peak_temp, common_name, species, species_ab, env_temp_min, env_temp_max) %>%
  write_csv(., "data/growth_analysis.csv", ";")
