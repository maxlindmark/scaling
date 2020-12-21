#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.23: Max Lindmark
#
# - Explore consumption and metabolic data
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
# magrittr_1.5       RColorBrewer_1.1-2 viridis_0.5.1      viridisLite_0.3.0  
# ggplot2_3.2.1      RCurl_1.95-4.12    bitops_1.0-6       readxl_1.3.1      
# tidylog_0.2.0      patchwork_0.0.1    tidyr_1.0.0        dplyr_0.8.3

con <- read_excel("data/consumption_data.xlsx")
met <- read_excel("data/metabolism_data.xlsx")

# Which cols to make numeric?
glimpse(con)
glimpse(met)

con$type <- NA

cols = c(15,16,17,18,19,20)

con[,cols] %<>% lapply(function(x) as.numeric(x))
met[,cols] %<>% lapply(function(x) as.numeric(x))

con <- con %>% rename(y = consumption)
met <- met %>% rename(y = metabolic_rate)

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- data.frame(bind_rows(con, met))

head(dat)
glimpse(dat)


# B. EXPLORE & CLEAN DATA ==========================================================
# Create abbreviated species name for plotting. Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Change settings for using scientific notation
options(scipen=999) 

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
ggplot(dat, aes(env_temp_mid, fill = common_name)) +
  geom_histogram() + 
  facet_wrap(~rate) + 
  coord_cartesian(expand = 0) +
  theme_classic() +
  NULL

# Convert experimental to Arrhenius scale:
dat$temp_arr <- 1/((dat$temp_c + 273.15) * 8.617332e-05)

# Standardize temperatures to median-reference temperature on C scale
dat$temp_norm <- dat$temp_c - dat$median_temp

# Calculate log mass
dat$log_mass <- log(dat$mass_g)

# Normalize mass with respect to max mass
dat$mass_norm <- dat$mass_g / dat$w_max_published_g

# Calculate log normalized mass
dat$log_mass_norm <- log(dat$mass_norm)


#** Plot general data ==============================================================
# Test which sizes I use
p1 <- ggplot(dat, aes(mass_norm, fill = species)) + 
  geom_histogram() + 
  facet_wrap(~ rate, nrow = 2, scales = "free_y") +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  coord_cartesian(expand = 0) + 
  labs(x = "Mass/Max mass") +
  guides(fill = FALSE) +
  NULL
pWord <- p1 + theme_classic() + theme(text = element_text(size = 12))
ggsave("figures/supp/data/meta_cons_size_range.png", width = 6.5, height = 6.5, dpi = 600)

# Plot size vs asymptotic size
# If we have a strong relationship between mass and max mass, then it could be argued that
# we also plot average interspecific relationships
ggplot(dat, aes(w_max_published_g, mass_g, color = species)) + 
  geom_point() + 
  facet_wrap(~ rate, nrow = 2, scales = "free_y") +
  scale_color_viridis(discrete = TRUE, option = "magma") +
  #coord_cartesian(expand = 0) + 
  #labs(x = "Mass/Max mass") +
  guides(fill = FALSE) +
  NULL


# Trophic level
p2 <- ggplot(dat, aes(x = reorder(species, trophic_level), y = trophic_level)) +
  geom_point(stat = 'identity', size = 2) +
  scale_fill_manual(name = "trophic_level") + 
  guides(colour = FALSE) +
  xlab("") + 
  ylab("Trophic level") + 
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 
pWord <- p2 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text.y = element_text(size = 8, face = "italic"))
ggsave("figures/supp/data/meta_cons_trophic_level.png", width = 6.5, height = 6.5, dpi = 600)


# Mid env. temperature (Fishbase) compared to experimental temperature range
pal <- brewer.pal("Dark2", n = 5)

pal[2] <- "grey"

p3 <- ggplot(dat) +
  coord_flip() +
  geom_jitter(aes(x = reorder(species, median_temp), 
                 y = temp_c, color = "Experiment"), size = 1, alpha = 0.6,
              width = 0.1, height = 0, shape = 21, fill = NA) +
  geom_point(aes(x = reorder(species, median_temp), 
                 y = median_temp, color = "Environment (median)"), size = 1, alpha = 0.6) +
  geom_point(aes(x = reorder(species, median_temp), 
                 y = env_temp_max, color = "Environment (max)"), size = 1, alpha = 0.6) +
  geom_point(aes(x = reorder(species, median_temp), 
                 y = env_temp_min, color = "Environment (min)"), size = 1, alpha = 0.6) +
  scale_color_manual(values = rev(pal), name = "Temperature") +
  xlab("") + 
  ylab(expression(paste("Temperature [", degree*C, "]"))) + 
  facet_wrap(~ rate) +
  NULL 
pWord <- p3 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text.y = element_text(size = 8, face = "italic"))
ggsave("figures/supp/data/meta_cons_temperatures.png", width = 6.5, height = 6.5, dpi = 600)


# Max. published weight
p4 <- ggplot(dat, aes(x = reorder(species, w_max_published_g), y = w_max_published_g)) +
  geom_point(stat = 'identity', size = 2) +
  scale_fill_manual(name = "w_max_published_g") + 
  guides(colour = FALSE) +
  xlab("") + 
  ylab("max published weight [g]") + 
  scale_y_continuous(trans = 'log10') +
  coord_flip() +
  facet_wrap(~ rate) +
  NULL 
pWord <- p4 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text.y = element_text(size = 8, face = "italic"))
ggsave("figures/supp/data/meta_cons_max_weight.png", width = 6.5, height = 6.5, dpi = 600)


# Phylogeny
p5 <- dat %>% dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(order, fill = family)) +
  geom_bar() +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  facet_wrap(~ rate) +
  NULL
pWord <- p5 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text.x = element_text(angle = 60, hjust = 1, size = 8,
                                                                 face = "italic"))
ggsave("figures/supp/data/meta_cons_phylogeny.png", width = 6.5, height = 6.5, dpi = 600)


# Biogeography
p6 <- dat %>% dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(biogeography, fill = biogeography)) +
  geom_bar() +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  facet_wrap(~ rate) +
  guides(fill = FALSE) +
  NULL
pWord <- p6 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("figures/supp/data/meta_cons_biogeography.png", width = 6.5, height = 6.5, dpi = 600)


# Lifestyle
p7 <- dat %>% 
  dplyr::distinct(common_name, .keep_all = TRUE) %>% 
  ggplot(., aes(habitat, fill  = lifestyle)) +
  geom_bar() +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  facet_wrap(~ rate) +
  NULL
pWord <- p7 + theme_classic() + theme(text = element_text(size = 12),
                                      axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("figures/supp/data/meta_cons_lifestyle.png", width = 6.5, height = 6.5, dpi = 600)


#** Plot response variable =========================================================
# Consumption
dat %>% filter(rate == "consumption") %>% 
  ggplot(., aes(temp_c, y)) + 
  geom_point(size = 2, alpha = 0.3) +
  facet_wrap(~species, scales = "free") +
  theme_classic(base_size = 11) +
  #stat_smooth(se = F) +
  guides(color = F) +
  NULL

# Metabolism
dat %>% filter(rate == "metabolism") %>% 
  ggplot(., aes(temp_c, y)) + 
  geom_point(size = 2, alpha = 0.3) +
  facet_wrap(~species, scales = "free") +
  theme_classic(base_size = 11) +
#  stat_smooth(se = F) +
  guides(color = F) +
  NULL

#**** Plot all data combined =======================================================
# Both together
dat$y_spec <- dat$y/dat$mass_g
dat$rate2 <- factor(dat$rate)
levels(dat$rate2) <- c("consumption [g/g/day]", "metabolism [mg O2/g/h]")

n_con <- nrow(filter(dat, rate == "consumption"))
n_met <- nrow(filter(dat, rate == "metabolism"))

dat$sample_size <- ifelse(dat$rate == "consumption",
                          paste("n=", n_con, sep=""),
                          paste("n=", n_met, sep=""))

# As functions of temperature
p8 <- ggplot(dat, aes(x = temp_arr, y = y_spec, fill = log10(mass_g))) + 
  geom_point(size = 2, alpha = 0.6, color = "white", shape = 21) +
  facet_wrap(~rate2, scales = "free", nrow = 2) +
  theme_classic(base_size = 11) +
  scale_y_log10() +
  labs(y = "mass-specific rate",
       x = "Arrhenius temperature [1/kT]") +
  scale_fill_viridis(option = "magma", name = "log10(mass)\n[g]") +
  geom_text(aes(x = Inf, y = Inf, label = sample_size, hjust = 1.05, vjust = 1.5), 
            color = "black") + 
  NULL
pWord <- p8 + theme_classic() + theme(text = element_text(size = 12))
ggsave("figures/supp/data/meta_cons_rate_temp.png", width = 6.5, height = 6.5, dpi = 600)
  

# As functions of mass
p9 <- ggplot(dat, aes(x = mass_g, y = y_spec, fill = temp_arr)) + 
  geom_point(size = 2, alpha = 0.6, color = "white", shape = 21) +
  facet_wrap(~rate2, scales = "free", nrow = 2) +
  theme_classic(base_size = 11) +
  labs(y = "mass-specific rate",
       x = "mass [g]") +
  scale_fill_viridis(option = "magma", name = "Arrhenius\ntemperature\n[1/kT]") +
  scale_y_log10() +
  scale_x_log10() +
  geom_text(aes(x = Inf, y = Inf, label = sample_size, hjust = 1.05, vjust = 1.5), 
            color = "black") + 
  NULL
pWord <- p9 + theme_classic() + theme(text = element_text(size = 12))
ggsave("figures/supp/data/meta_cons_rate_mass.png", width = 6.5, height = 6.5, dpi = 600)


# As functions of mass - color by species
dat %>% filter(rate == "consumption" 
               #& mass_g < 3
               & above_peak_temp == "N"
               ) %>% 
  ggplot(., aes(x = mass_g, y = y_spec, fill = species)) + 
  geom_point(size = 2, alpha = 1, color = "white", shape = 21) +
  facet_wrap(~rate2, scales = "free", nrow = 2) +
  theme_classic(base_size = 11) +
  labs(y = "mass-specific rate",
       x = "mass [g]") +
  scale_fill_viridis(#option = "magma", 
                     name = "Arrhenius\ntemperature\n[1/kT]", discrete = TRUE) +
  scale_y_log10() +
  scale_x_log10() +
  geom_text(aes(x = Inf, y = Inf, label = sample_size, hjust = 1.05, vjust = 1.5), 
            color = "black") + 
  NULL



#** Loop through species and save plots ==========================================================
# This is for identifying which data points are below optimum and for inspecting data by species
t <- c()

# Consumption
s_datc <- filter(dat, rate == "consumption")
for(i in unique(s_datc$common_name)) {
  
  t <- filter(s_datc, common_name == i)
  t$round_mass_g <- round(t$mass_g, digits = 0)
  t$size_cl <- cut(t$round_mass_g, 10)
  title <- t$common_name[1]
  
  p <- t %>% 
    ggplot(., aes(temp_c, y, color = factor(size_cl)))+ 
    #geom_jitter(size = 4, height = 0) +
    geom_point(size = 4) +
    geom_line() +
    theme_classic(base_size = 11) + 
    scale_color_viridis(discrete = TRUE) +
    ggtitle(title)
    
  path <- paste("figures/supp/data/species_plots/con_met/con_", i, ".pdf", sep = "")
  ggsave(path, plot = p, scale = 1, width = 20, height = 20, units = "cm")
  
}  

# Metabolism
t <- c()

s_datm <- filter(dat, rate == "metabolism")
for(i in unique(s_datm$common_name)) {
  
  t <- filter(s_datm, common_name == i)
  t$round_mass_g <- round(t$mass_g, digits = 0)
  t$size_cl <- cut(t$round_mass_g, 10)
  title <- t$common_name[1]
  
  p <- t %>% 
    ggplot(., aes(temp_c, y, color = factor(size_cl)))+ 
    #geom_jitter(size = 4, height = 0) +
    geom_point(size = 4) +
    geom_line() +
    theme_classic(base_size = 11) + 
    scale_color_viridis(discrete = TRUE) +
    ggtitle(title)
  
  path <- paste("figures/supp/data/species_plots/con_met/met_", i, ".pdf", sep = "")
  ggsave(path, plot = p, scale = 1, width = 20, height = 20, units = "cm")
  
}  


# Calculate number of unique studies
# > length(unique(s_datm$reference))
# [1] 26
# > length(unique(s_datc$reference))
# [1] 20


# C. SAVE DATA =====================================================================
glimpse(s_datm)
s_datm %>%
 select(y, mass_g, log_mass, mass_norm, log_mass_norm, temp_c, temp_arr, median_temp, 
        above_peak_temp, common_name, species, species_ab, unit, type) #%>%
#write_csv(., "data/met_analysis.csv", ";")


glimpse(s_datc) 
s_datc %>%
  select(y, mass_g, log_mass, mass_norm, log_mass_norm, temp_c, temp_arr, median_temp, 
         above_peak_temp, common_name, species, species_ab, unit, type) #%>%
#write_csv(., "data/con_analysis.csv", ";")


