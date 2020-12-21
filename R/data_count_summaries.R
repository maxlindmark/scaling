#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# Short script for counting summaries of the data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(RCurl)
library(readxl)
library(magrittr)


# B. READ IN DATA ==================================================================
# Metabolism and consumption =======================================================
# Read in your data file
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

length(unique(con$common_name))

test <- con %>% filter(above_peak_temp == "Y")
length(unique(test$common_name))

# Filter data points at below optimum temperatures
met <- met %>% filter(above_peak_temp == "N")
con <- con %>% filter(above_peak_temp == "N")

# Count data
length(unique(met$common_name))
length(unique(con$common_name))
nrow(met)
nrow(con)
mean(met$temp_c)
mean(con$temp_c)
n_stand <- nrow(filter(met, type == "Standard"))
n_rout_rest <- nrow(filter(met, type %in% c("Resting", "Routine")))
n_stand / (n_stand + n_rout_rest)
1 - (n_stand / (n_stand + n_rout_rest))

summary(met$mass_g)
summary(con$mass_g)

met %>% group_by(common_name, temp_c) %>% summarise(n = n()) %>% as.data.frame()

str(met)

# How many temperatures? 
met %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  as.data.frame()

con %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  as.data.frame()

# Average # of temperatures per species?
met %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_temp))

con %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_temp))

# How many masses? 
met %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  as.data.frame()

con %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  as.data.frame()

# Average # of masses per species?
met %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_mass))

con %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_mass))


# Growth ===========================================================================
dat <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

str(dat)

# Count data
length(unique(dat$common_name))
nrow(dat)

# How many temperatures? 
dat %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  as.data.frame()

# Average # of temperatures per species?
dat %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_temp))

# How many masses? 
dat %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  as.data.frame()

# Average # of masses per species?
dat %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_mass))


# T_opt ============================================================================
# Read in your data file
topt <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/topt_analysis.csv"))

# Count data
length(unique(topt$common_name))
nrow(topt)
