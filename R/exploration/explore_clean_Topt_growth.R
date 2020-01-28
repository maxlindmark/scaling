#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore body growth rate data
# 
# A. Load libraries & read data
#
# B. Explore growth data
#
# C. Explore growth data
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
# [1] patchwork_0.0.1    magrittr_1.5       RColorBrewer_1.1-2 viridis_0.5.1     
# [5] viridisLite_0.3.0  ggplot2_3.1.1      RCurl_1.95-4.12    bitops_1.0-6      
# [9] readxl_1.3.1       tidylog_0.1.0      tidyr_0.8.3        dplyr_0.8.1    


# ***  Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

dat <- read_excel("data/growth_data_Topt.xlsx")

glimpse(dat)

# Which cols to make numeric?
cols = c(1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)

colnames(dat)[1] <-  "G"

# B. EXPLORE DATA ==================================================================
#** General ========================================================================
# Trophic level
ggplot(dat, aes(x = reorder(common_name, trophic_level), y = trophic_level)) +
  geom_point(stat = 'identity', size=6) +
  scale_fill_manual(name = "trophic_level") + 
  theme_classic(base_size = 15) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("Trophic level") + 
  coord_flip() +
  NULL 

# Max. published weight
ggplot(dat, aes(x = reorder(common_name, w_max_published_g), 
                y = log10(w_max_published_g))) +
  geom_point(stat = 'identity', size=6) +
  scale_fill_manual(name = "w_max_published_g") + 
  theme_classic(base_size = 15) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("log10(max published weight) [g]") + 
  coord_flip() +
  NULL 

# Phylogeny
nb.cols <- length(unique(dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

dat %>% distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(order, fill = family)) +
  geom_bar() +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = mycolors) +
  NULL

# Biogeography
dat %>% distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(biogeography, fill = biogeography)) +
  geom_bar() +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = FALSE) +
  scale_fill_manual(values = mycolors) +
  NULL

# Lifestyle
dat %>% 
  distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(habitat, fill  = lifestyle)) +
  geom_bar() +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = mycolors) +
  NULL


# C. CLEAN DATA ==================================================================
# Calculate mean optimum temperature within species
dat$mean_opt_temp_c <- ave(dat$opt_temp_c, dat$common_name)

# Center each size class' optimum temperature by mean optimum temperature for that species
dat$opt_temp_c_ct <- dat$opt_temp_c - dat$mean_opt_temp_c

# We use geometric mean for size, and if this is not possible we'll go with mid point of size range. # Mass for analysis:
dat$mass <- dat$geom_mean_mass_g

# Replace NA with -9...
dat$mass[is.na(dat$mass)] <- -9

dat$mass <- ifelse(dat$mass == -9,
                   dat$size_group,
                   dat$mass)

# Check we only have real values:
unique(is.na(dat$mass))

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass / dat$w_max_published_g

# Calculate log mass
dat$log_mass_norm <- log(dat$mass_norm)

# Mean center mass
dat$log_mass_norm_ct <- dat$log_mass_norm - mean(dat$log_mass_norm)

# Save data:
str(dat)
head(dat)

# Plot
plot(opt_temp_c_ct ~ log_mass_norm_ct, data = dat)

# dat %>% 
#   select(G, geom_mean_mass_g, opt_temp_c, common_name, species,  
#          mean_opt_temp_c, opt_temp_c_ct, mass, mass_norm, log_mass_norm, log_mass_norm_ct) %>% 
#   write_csv(., "data/topt_analysis.csv", ";")


#** Extra ==========================================================================
# Growth rate at optimum over size 
ggplot(dat, aes(log10(mass_norm), G, 
                  color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm", se = FALSE, size = 2, alpha = 0.2)+ 
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

ggplot(dat, aes(log(mass), log(G), 
                  color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm", se = FALSE, size = 2, alpha = 0.2)+ 
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

summary(lm(log(s_dat$growth_rate) ~ log(s_dat$mass)))
