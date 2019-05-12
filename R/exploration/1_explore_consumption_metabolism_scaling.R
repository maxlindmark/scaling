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
          "RColorBrewer",
          "magrittr")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
pkg_info(pkgs)

# package loadedversion
# 1         dplyr       0.8.0.1
# 2       ggplot2         3.1.1
# 3      magrittr           1.5
# 4          plyr         1.8.4
# 5  RColorBrewer         1.1-2
# 6         RCurl     1.95-4.12
# 7        readxl         1.3.1
# 8       tidylog         0.1.0
# 9         tidyr         0.8.3
# 10      viridis         0.5.1

#==** Read data ====
# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

con <- read_excel("data/consumption_scaling_data.xlsx")
met <- read_excel("data/metabolism_scaling_data.xlsx")

glimpse(dat)

cols = c(1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21)

con[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
met[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

con$rate <- "consumption"
met$rate <- "metabolism"

# Calculate mean b within species and rate
con <- con %>% 
  group_by(species, rate) %>% 
  mutate(mean_b = ave(b, common_name))

met <- met %>% 
  group_by(species, rate) %>% 
  mutate(mean_b = ave(b, common_name))

dat <- rbind(con, met)

dat$mean_b_ct <- dat$b - dat$mean_b

#======== B. EXPLORE DATA ========
#==** Normalize variables ==== 
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate)

# Relative to temp in environment
dat$env_temp_mid_norm <- dat$temp_c - dat$env_temp_mid

# Relative to "preferred" temp (not all species have this info)
dat$pref_temp_mid_norm <- dat$temp_c - dat$pref_temp_mid

#==** Plot general data ==== 
#-- Trophic level
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

#-- Temperature midpoint (Fishbase)
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

#-- Mid env. temperature (Fishbase) compared to experimental temperature range
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

#-- Max. published weight
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

#-- Phylogeny
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

#-- Biogeography
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

#-- Lifestyle
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
#==**** Interspecific ==== 
#-- All exponents
ggplot(dat, aes(temp_c, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#-- Significant weight weight effects
ggplot(filter(dat, significant_size == "Y"), aes(temp_c, b)) + 
  ylim(0, 1.2) +
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  facet_wrap(~ rate) +
  NULL

#-- Color by species and add species-lines
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

#-- Temperature centered to env. midpoint from fishbase
ggplot(filter(dat, significant_size == "Y"), 
       aes(env_temp_mid_norm, b)) + 
  geom_point(size = 6, fill = "black", 
             color = "white", shape = 21, alpha = 0.5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#-- Temperature centered to env. midpoint from fishbase, species lines
ggplot(filter(dat, significant_size == "Y"), 
       aes(env_temp_mid_norm, b, color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#-- Filter species with pref temp and plot b~pref_temp (Appendix analysis)
ggplot(filter(dat, significant_size == "Y" & pref_temp_mid > 0), 
       aes(pref_temp_mid, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm") +
  theme_classic(base_size = 15) +
  facet_wrap(~rate) +
  NULL

#-- Filter species with pref temp and plot b~pref_temp_centered
ggplot(filter(dat, significant_size == "Y" & pref_temp_mid > 0), 
       aes(pref_temp_mid_norm, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm") +
  theme_classic(base_size = 15) +
  facet_wrap(~rate) +
  NULL

#-- b vs log max size
ggplot(dat, aes(log(w_max_published_g), b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.7) + 
  facet_wrap(~rate) +
  NULL

#-- b vs trophic level
ggplot(dat, aes(trophic_level, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  facet_wrap(~rate) +
  NULL

#-- b vs lifestyle
ggplot(dat, aes(lifestyle, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  geom_boxplot() +
  facet_wrap(~rate) +
  NULL

#-- b vs biogeography
ggplot(dat, aes(biogeography, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  geom_boxplot() +
  facet_wrap(~rate) +
  NULL

#-- b vs genus
ggplot(dat, aes(genus, b)) + 
  theme_classic(base_size = 15) +
  geom_point(size = 5, alpha = 0.7) +
  geom_boxplot() +
  facet_wrap(~rate) +
  NULL


#==**** Intraspecific ==== 
s_dat <- dat %>% 
  filter(significant_size == "Y") %>% 
  group_by(common_name, rate) %>% 
  filter(n()>1)

data.frame(s_dat)

#-- Plot all intra-specific data
ggplot(s_dat, aes(temp_c, b)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  facet_wrap(~rate) +
  NULL

#-- Species specific lines
nb.cols <- length(unique(s_dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat, aes(temp_c, b, color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#-- Temperature centered to env. midpoint from fishbase
ggplot(s_dat, aes(env_temp_mid_norm, b)) + 
  geom_point(size = 6, fill = "black", 
             color = "white", shape = 21, alpha = 0.5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#-- Centered b across species
ggplot(s_dat, aes(env_temp_mid_norm, mean_b_ct)) + 
  geom_point(size = 6, fill = "black", 
             color = "white", shape = 21, alpha = 0.5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

#-- Temperature centered to env. midpoint from fishbase, species lines
ggplot(s_dat, aes(env_temp_mid_norm, b, color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL

s_dat %>% 
  filter(rate == "consumption") %>%
  ggplot(., aes(env_temp_mid_norm, b, color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F, size = 2) +
  scale_color_manual(values = mycolors) + 
  NULL


#-- Comparing Cmax and metabolism exponents
# Rainclouds with boxplots
# source code from github:
script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/raincloud_plot.R", ssl.verifypeer = FALSE)

eval(parse(text = script))

nb.cols <- 8

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat, aes(x = rate, y = b, fill = rate, colour = rate))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE, alpha = 0.8)+
  geom_point(position = position_jitter(width = .15), size = 3, alpha = 1)+
  geom_boxplot(aes(x = rate, y = b),
               outlier.shape = NA, alpha = 0.3, width = .2, color = "black", size = 1) +
  coord_flip() + 
  guides(fill = FALSE, colour = FALSE) +
  scale_color_manual(values = mycolors) + 
  scale_fill_manual(values = mycolors) + 
  theme_classic(base_size = 20) +
  labs(y = "Size-scaling exponent", x = "Rate") +
  NULL

summary(lm(b~rate, data=s_dat))

  



