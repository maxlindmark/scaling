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

# Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "readxl",
          "RCurl",
          "ggplot2",
          "viridis",
          "RColorBrewer",
          "magrittr")

# Install packages
#install.packages(pkgs)

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# pkg_info(pkgs)

# package loadedversion
# 1        dplyr         0.8.1
# 2      ggplot2         3.1.1
# 3     magrittr           1.5
# 4 RColorBrewer         1.1-2
# 5        RCurl     1.95-4.12
# 6       readxl         1.3.1
# 7      tidylog         0.1.0
# 8        tidyr         0.8.3
# 9      viridis         0.5.1

# Also add patchwork
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

con <- read_excel("data/consumption_data.xlsx")
met <- read_excel("data/metabolism_data.xlsx")

# Which cols to make numeric?
cols = c(1,2,3,12,13,14,15,16,17)

con[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
met[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

colnames(con)
head(con)

con <- con %>% rename(y = consumption)
met <- met %>% rename(y = metabolic_rate)

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- bind_rows(con, met)

head(dat)


# B. EXPLORE DATA ==================================================================
#** Normalize variables ============================================================
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate)

# Relative to temp in environment
dat$env_temp_mid_norm <- dat$temp_c - dat$env_temp_mid

# Relative to "preferred" temp (not all species have this info)
dat$pref_temp_mid_norm <- dat$temp_c - dat$pref_temp_mid

# Create inverse temp
dat$inv_temp <- 1/((dat$temp_c + 273.15) * 8.617332e-05)

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass_g / dat$w_max_published_g
# ---- Stechlin cisco har larger size than max*

# Create data with with species that have also temperature repliates within species (intra-specific analysis)
s_datc <- data.frame(
  dat %>% 
    filter(rate == "consumption") %>% 
    group_by(species) %>% 
    mutate(unique_t = n_distinct(env_temp_mid_norm)) %>% 
    filter(unique_t > 1)
)

s_datm <-  data.frame(
  dat %>% 
    filter(rate == "metabolism") %>% 
    group_by(species) %>% 
    mutate(unique_t = n_distinct(env_temp_mid_norm)) %>% 
    filter(unique_t > 1)
)


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
pal <- viridis(n = 5)
  
dat %>% group_by(common_name) %>% 
  ggplot(., aes(x = reorder(common_name, env_temp_mid), 
                y = env_temp_mid, color = "blue")) +
  geom_point(stat = 'identity', size = 5) +
  geom_errorbar(aes(reorder(common_name, env_temp_mid), 
                    ymin = env_temp_min, 
                    ymax = env_temp_max, color = "blue"), 
                size = 3,
                width = 0) +
  geom_point(aes(x = reorder(common_name, env_temp_mid), 
                 y = temp_c, color = "gray"), size = 3) +
  scale_color_manual(labels = c("Mid. Env. Temperature [C]", 
                                "Experimental\ntemperature"), 
                     values = pal[c(3, 1)]) + # c("#f1a340", "#998ec3")
  scale_fill_manual(name = "env_temp_mid") + 
  theme_classic(base_size = 18) +
  guides(color = guide_legend(title = "")) +
  xlab("") + 
  ylab("Temperature [C]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# Some species without environmental temperature. Can look for other source. 
# Meanwhile, check if they have a preferred temperature
# dat$t <- is.na(dat$env_temp_mid)
# t <- subset(dat, t == TRUE)
# unique(t$common_name)
# t$pref_temp_mid
# [1] "Japanese flounder"      "Orange-spotted grouper" "Stechlin cisco"        
# [4] "Southern catfish"       "Spotted grunter"   
# will need to fix that

# Max. published weight
ggplot(dat, aes(x = reorder(common_name, w_max_published_g), 
                y = log10(w_max_published_g))) +
  geom_point(stat = 'identity', size = 4) +
  scale_fill_manual(name = "w_max_published_g") + 
  theme_classic(base_size = 15) +
  guides(colour = FALSE) +
  xlab("") + 
  ylab("log10(max published weight) [g]") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 

# Phylogeny
nb.cols <- length(unique(dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

dat %>% distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(order, fill = family)) +
  geom_bar() +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = mycolors) +  
  facet_wrap(~ rate) +
  NULL
# Doesn't seem like we can do much about phyologey here..

# Biogeography
nb.cols <- length(unique(dat$biogeography))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

dat %>% distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(biogeography, fill = biogeography)) +
  geom_bar() +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  #scale_fill_manual(values = mycolors) + 
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ rate) +
  NULL

# Lifestyle
nb.cols <- length(unique(dat$lifestyle))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

dat %>% 
  distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(habitat, fill  = lifestyle)) +
  geom_bar() +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = mycolors) + 
  facet_wrap(~ rate) +
  NULL


#** Plot response variable =========================================================
# All data, by species
dat %>% 
  filter(rate == "consumption") %>% 
  ggplot(., aes(env_temp_mid_norm, y)) + 
  geom_point(size = 2, alpha = 0.7, color = pal[1]) +
  theme_classic(base_size = 11) +
  facet_wrap(~species, scales = "free_y") +
  stat_smooth(se = F) +
  guides(color = F) +
  NULL

# Only intraspecific data, by species
s_datc %>% 
  ggplot(., aes(env_temp_mid_norm, y))+ 
  geom_point(size = 2, alpha = 0.7, color = pal[1]) +
  theme_classic(base_size = 11) +
  facet_wrap(~species, scales = "free_y") +
  stat_smooth(se = F) +
  guides(color = F) +
  NULL

# This looks promising! For species with a large temperature range there is an optimum. Now, normalize all data and plot together. Or at least filter the most common unit. I could plot the actual rates, but it if I do that the optimum dissapears. I suspect this is driven by species having very different Cmax. Hence I normalize to max within species below.
ggplot(s_datc, aes(unit)) + geom_bar()
ggplot(s_datm, aes(unit)) + geom_bar()

pal <- viridis(n = 7)

# Group by normalized mass
# Consumption
p1 <- s_datc %>% 
  group_by(species) %>% 
  mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(env_temp_mid_norm, y_norm)) + 
  theme_classic(base_size = 14) +
  geom_point(size = 2, alpha = 0.4, color = pal[4]) +
  guides(color = FALSE) +
  labs(x = "Normalized temperature", y = "Normalized consumption") +
  stat_smooth(span = 1.0, se = F, size = 2, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 1.2, se = F, size = 2, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 0.8, se = F, size = 2, alpha = 0.4, geom = "line", color = pal[2]) +
  theme(aspect.ratio = 1) +
  coord_cartesian(xlim = c(min(c(min(s_datm$env_temp_mid_norm), min(s_datc$env_temp_mid_norm))),
                           max(c(max(s_datm$env_temp_mid_norm), max(s_datc$env_temp_mid_norm))))) +
  NULL

# Metabolism
p2 <- s_datm %>% 
  filter(unit == "mg O2/h") %>% 
  group_by(species) %>% 
  mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(env_temp_mid_norm, y_norm)) + 
  theme_classic(base_size = 14) +
  geom_point(size = 2, alpha = 0.4, color = pal[4]) +
  guides(color = FALSE) +
  labs(x = "Normalized temperature", y = "Normalized metabolic rate") +
  stat_smooth(span = 1.0, se = F, size = 2, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 1.2, se = F, size = 2, alpha = 0.4, geom = "line", color = pal[2]) +
  stat_smooth(span = 0.8, se = F, size = 2, alpha = 0.4, geom = "line", color = pal[2]) +
  theme(aspect.ratio = 1) +
  coord_cartesian(xlim = c(min(c(min(s_datm$env_temp_mid_norm), min(s_datc$env_temp_mid_norm))),
                           max(c(max(s_datm$env_temp_mid_norm), max(s_datc$env_temp_mid_norm))))) +
  NULL

p1 + p2

# ggsave("figs/opt_curves.pdf", plot = last_plot(), scale = 1, width = 20, height = 20, units = "cm")


# Now plot slope log(rate)~temp, where the slope is the activation energy, using sub-optimum temperatures! (<11, based on previous plots). Note I don't normalize here but plot the actual rates
# Consumption
p3 <- s_datc %>% 
  filter(env_temp_mid_norm < 12) %>% 
  ggplot(., aes(inv_temp, log(y))) + 
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.4, color = pal[2]) +
  labs(x = "Inverse temperature [1/kT]", y = "log(consumption [g/day])") +
  stat_smooth(method = "lm", size = 3, alpha = 0.4, geom = "line", color = pal[3]) +
  NULL

# Metabolism
p4 <- s_datm %>% 
  filter(unit == "mg O2/h" & env_temp_mid_norm < 12) %>% 
  ggplot(., aes(inv_temp, log(y))) + 
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.4, color = pal[2]) +
  labs(x = "Inverse temperature [1/kT]", y = "log(metabolic rate [mg O2/h])") +
  stat_smooth(method = "lm", size = 3, alpha = 0.4, geom = "line", color = pal[3]) +
  NULL

p3 + p4

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
pal <- viridis(n = nranks)
pal2 <- colorRampPalette(brewer.pal(8, "Dark2"))(nranks)

xmax <- max(c(max(s_datc$env_temp_mid_norm), max(s_datm$env_temp_mid_norm)))+0.1
xmin <- min(c(min(s_datc$env_temp_mid_norm), min(s_datm$env_temp_mid_norm)))+0.1

# Group by normalized mass
# Consumption
p3 <- s_datc %>% 
  mutate(quartile = ntile(mass_g, nranks)) %>% 
  group_by(species) %>% 
  mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(env_temp_mid_norm, y_norm, color = factor(quartile))) + 
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.4) +
  guides(color = FALSE) +
  labs(x = "Normalized temperature", y = "Normalized consumption") +
  scale_color_manual(values = rev(pal)) +
  #scale_color_manual(values = pal2) +
  stat_smooth(span = 1.0, se = F, size = 3, alpha = 1, geom = "line") +
  coord_cartesian(xlim = c(xmin, xmax)) +
  NULL

# Metabolism
p4 <- s_datm %>% 
  mutate(quartile = ntile(mass_g, nranks)) %>% 
  group_by(species) %>% 
  mutate(y_norm = y/max(y)) %>% 
  ggplot(., aes(env_temp_mid_norm, y_norm, color = factor(quartile))) + 
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.4) +
  labs(x = "Normalized temperature", y = "Normalized metabolism") +
  scale_color_manual(values = rev(pal)) +
  #scale_color_manual(values = pal2) +
  stat_smooth(span = 1.0, se = F, size = 3, alpha = 1, geom = "line") +
  coord_cartesian(xlim = c(xmin, xmax)) +
  NULL

p3 + p4


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



