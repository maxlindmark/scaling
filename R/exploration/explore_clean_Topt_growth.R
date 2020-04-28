#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore optimum growth temperature data. See "explore_clean_growth" for 
#   addition biological data 
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

dat <- read_excel("data/growth_data_Topt.xlsx")

glimpse(dat)

# Which cols to make numeric?
cols = c(15:21)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)

dat$y <- dat$`growth_rate_%/day`


# B. EXPLORE & CLEAN DATA ==================================================================
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

# Calculate mean optimum temperature within species
dat$mean_opt_temp_c <- ave(dat$opt_temp_c, dat$common_name)

# Center each size class' optimum temperature by mean optimum temperature for that species
dat$opt_temp_c_ct <- dat$opt_temp_c - dat$mean_opt_temp_c

# Use size_group if no geometric mean mass
dat$mass_g <- dat$geom_mean_mass_g

dat$mass_g[is.na(dat$mass_g)] <- -9

dat$mass_g <- ifelse(dat$mass_g == -9,
                     dat$size_group,
                     dat$mass_g)

# Calculate log mass
dat$log_mass <- log(dat$mass_g)

# Normalize mass with respect to max mass
dat$mass_norm_max <- dat$mass_g / dat$w_max_published_g

# Normalize mass with respect to mass at maturation
dat$mass_norm_mat <- dat$mass_g / dat$w_maturation_g



#** Plot general data ==============================================================
# Test which sizes I use
p1 <- ggplot(dat, aes(mass_norm_mat, fill = species)) + 
  geom_histogram() + 
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  coord_cartesian(expand = 0) + 
  labs(x = "Mass/Maturation mass") +
  guides(fill = FALSE) +
  NULL
pWord <- p1 + theme_classic() + theme(text = element_text(size = 12), aspect.ratio = 1)
ggsave("figures/supp/data/topt_size_range.png", width = 6.5, height = 6.5, dpi = 600)


#** Plot response variable =========================================================
# Optimum growth temperature as function of mass
ggplot(dat, aes(log(mass_norm_mat), opt_temp_c_ct, fill = species)) + 
  geom_point(size = 3, alpha = 1, color = "white", shape = 21) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_fill_viridis(discrete = T, option = "magma") +
  labs(x = "ln(mass) [g]", 
       y = expression(paste("Rescaled temperature [", degree*C, "]"))) +
  theme(aspect.ratio = 3/4,
        legend.text = element_text(size = 8)) +
  NULL

ggplot(dat, aes(log(mass_norm_max), opt_temp_c_ct, fill = species)) + 
  geom_point(size = 3, alpha = 1, color = "white", shape = 21) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_fill_viridis(discrete = T, option = "magma") +
  labs(x = "ln(mass) [g]", 
       y = expression(paste("Rescaled temperature [", degree*C, "]"))) +
  theme(aspect.ratio = 3/4,
        legend.text = element_text(size = 8)) +
  NULL


# C. SAVE DATA ==================================================================
glimpse(dat)
dat %>%
  select(y, `growth_rate_%/day`, geom_mean_mass_g, size_group, mass_g, mass_norm_max,
         mass_norm_mat, opt_temp_c, opt_temp_c_ct, median_temp, common_name, species, species_ab) %>%
  write_csv(., "data/topt_analysis.csv", ";")
