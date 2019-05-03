#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore consumption and metabolic scaling data
# 
# A. Load libraries & read data
#
# B. Explore data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES & READ DATA ========
rm(list = ls())

#==** Load packages ====
# Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "readxl",
          "RCurl",
          "ggplot2",
          "viridis",
          "plyr",
          "RColorBrewer")

# Install packages
#install.packages(pkgs)

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# x <- devtools::session_info(pkgs = pkgs)
# x <- as.data.frame(x$packages)
# x <- dplyr::filter(x, package %in% pkgs) %>% 
# dplyr::select(-`*`, -date, -source) %>% 
# dplyr::arrange(package)
# x

# package   version
# 1        dplyr   0.8.0.1
# 2      ggplot2     3.0.0
# 3         plyr     1.8.4
# 4 RColorBrewer     1.1-2
# 5        RCurl 1.95-4.10
# 6       readxl     1.3.1
# 7      tidylog     0.1.0
# 8        tidyr     0.8.3
# 9      viridis     0.5.1

#==** Read data ====
# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

con <- read_excel("data/consumption_scaling_data.xlsx")
met <- read_excel("data/metabolism_scaling_data.xlsx")

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- rbind(con, met)

glimpse(dat)

cols = c(1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)

#======== B. EXPLORE DATA ========
#==** Normalize variables ==== 
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate)

# Relative to temp in environment
dat$env_temp_mid_norm <- dat$temp_c - dat$env_temp_mid


#==** Plot general data ==== 
##-- Trophic level
ggplot(dat, aes(x = reorder(common_name, trophic_level), y = trophic_level)) +
  geom_point(stat = 'identity', size=6) +
  scale_fill_manual(name = "trophic_level") + 
  theme_classic(base_size = 15) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("Trophic level") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

##-- Temperature midpoint (Fishbase)
ggplot(dat, aes(x = reorder(common_name, env_temp_mid), 
                y = env_temp_mid)) +
  geom_point(stat = 'identity', size=5) +
  scale_fill_manual(name = "env_temp_mid") + 
  theme_classic(base_size = 15) +
  guides(colour = T) +
  xlab("") + 
  ylab("Mid. Env. Temperature [C]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

##-- Mid env. temperature (Fishbase) compared to experimental temperature range
dat %>% group_by(common_name) %>% 
  ggplot(., aes(x = reorder(common_name, env_temp_mid), 
                y = env_temp_mid, color = "blue")) +
  geom_point(stat = 'identity', size=5) +
  geom_errorbar(aes(reorder(common_name, env_temp_mid), 
                    ymin = env_temp_min, 
                    ymax = env_temp_max, color = "blue"), 
                size = 3,
                width = 0) +
  geom_point(aes(x = reorder(common_name, env_temp_mid), 
                 y = temp_c, color = "gray"), size = 3) +
  scale_color_manual(labels = c("Mid. Env. Temperature [C]", 
                                "Experimental\ntemperature"), 
                     values = c("#f1a340", "#998ec3")) +
  scale_fill_manual(name = "env_temp_mid") + 
  theme_classic(base_size = 18) +
  guides(color = guide_legend(title = "")) +
  xlab("") + 
  ylab("Temperature [C]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# Some species without environmental temperature. Can look for other source. Meanwhile, check if they have a preferred temperature
dat$t <- is.na(dat$env_temp_mid)
t <- subset(dat, t == TRUE)
data.frame(t)

# Utah chub, Stechlin cisco, Southern catfish have no temperature at all, will need to fix that

##-- Max. published weight
ggplot(dat, aes(x = reorder(common_name, w_max_published_g), 
                y = log10(w_max_published_g))) +
  geom_point(stat = 'identity', size=6) +
  scale_fill_manual(name = "w_max_published_g") + 
  theme_classic(base_size = 15) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("log10(max published weight) [g]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

##-- Phylogeny
nb.cols <- length(unique(dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

dat %>% distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(order, fill = family)) +
  geom_bar() +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = mycolors) +  
  facet_wrap(~ rate) +
  NULL
# Doesn't seem like we can do much about phyologey here..

##-- Biogeography
nb.cols <- length(unique(dat$biogeography))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

dat %>% distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(biogeography, fill = biogeography)) +
  geom_bar() +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = FALSE) +
  scale_fill_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

##-- Lifestyle
nb.cols <- length(unique(dat$lifestyle))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

dat %>% 
  distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(habitat, fill  = lifestyle)) +
  geom_bar() +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL


#==** Plot response variable ==== 
##-- All exponents
ggplot(dat, aes(temp_c, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

##-- Significant weight weight effects
ggplot(filter(dat, significant_size == "Y"), aes(temp_c, b)) + 
  ylim(0, 1.2) +
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  facet_wrap(~ rate) +
  NULL

##-- Color by species and add species-lines
nb.cols <- length(unique(dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(filter(dat, significant_size == "Y"), aes(temp_c, b, color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#### CONTINIUE FROM HERE ####





##-- Temperature centered to env. midpoint from fishbase
ggplot(filter(dat, significant_size == "Y"), aes(env_temp_mid_ct, b)) + 
  ylim(0, 1.2) +
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm") +
  theme_classic(base_size = 15) +
  NULL

# temperature centered to env. midpoint from fishbase, species grouping
ggplot(filter(dat, significant_size == "Y"), 
       aes(env_temp_mid_ct, b, color = species)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  ylim(0, 1.2) +
  stat_smooth(method = "lm", se = F, alpha = 0.2) +
  stat_smooth(aes(group = lm_all), method = "lm", 
              colour = "black", size = 3, alpha = 0.2) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# filter species with pref temp and plot b~pref_temp
ggplot(filter(dat, significant_size == "Y" & pref_temp_mid > 0), 
       aes(pref_temp_mid, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm") +
  theme_classic(base_size = 15) +
  NULL

# filter species with pref temp and plot b~pref_temp_centered
ggplot(filter(dat, significant_size == "Y" & pref_temp_mid > 0), 
       aes(pref_temp_mid_ct, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm") +
  theme_classic(base_size = 15) +
  NULL


##-- Intraspecific data:
s_dat <- dat %>% 
  filter(significant_size == "Y") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# plot all intra-specific data
ggplot(s_dat, aes(temp_c, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  NULL

# add in species-specific lines
colourCount <- length(unique(s_dat$common_name)) # number of levels

ggplot(s_dat, aes(temp_c, b, color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# plotting temp as difference between fishbase temp
ggplot(s_dat, aes(env_temp_mid_ct, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  NULL

summary(lm(s_dat$b ~ s_dat$env_temp_mid_ct))

# add in species colour
ggplot(s_dat, aes(env_temp_mid_ct, b, color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

ggplot(s_dat, aes(pref_temp_mid_ct, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  NULL

# plotting temp as difference between mean in experiment
ggplot(s_dat, aes(temp_c_ct, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  NULL

## CORRECT
ggplot(s_dat, aes(temp_c_ct2, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  NULL

## Not a good methods anyway because of the plot showing they experimental temepratures are differently spread around environmental temperature

# plot normalized constant a within species
s_dat %>% 
  group_by(species) %>%
  mutate(a_norm = (a - min(a))/(max(a) - min(a))) %>%   
  ggplot(., aes(env_temp_mid_ct, a_norm)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.5) +
  stat_smooth(span = 1, se = F, size = 3, col = "black") +
  #scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# species unique slopes
s_dat %>% 
  group_by(common_name) %>%
  filter(max(a) < 1000) %>% 
  mutate(a_norm = (a - min(a))/(max(a) - min(a))) %>%   
  ggplot(., aes(env_temp_mid_ct, a_norm, color = common_name)) + 
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~species, scales = "free") +
  theme_classic(base_size = 15) +
  stat_smooth(se = F, span = 1) +
  guides(color = F) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# b vs log max size
ggplot(dat, aes(log(w_max_published_g), b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# b vs trophic level
ggplot(dat, aes(trophic_level, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# b vs lifestyle
ggplot(dat, aes(lifestyle, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  geom_boxplot() +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# b vs biogeography
ggplot(dat, aes(biogeography, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  geom_boxplot() +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# b vs genus
ggplot(dat, aes(genus, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  geom_boxplot() +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# b vs difference between max size and max size in experiment
colourCount <- length(unique(dat$common_name)) # number of levels
getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))

dat %>% 
  mutate(size_ratio = max_mass_g/w_max_published_g) %>% 
  ggplot(., aes(size_ratio, b, color = common_name)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL




