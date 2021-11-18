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
library(tidyverse)


# B. READ IN DATA ==================================================================
# Metabolism and consumption =======================================================
# Read in your data file
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

# Which species are in both data sets?
con_spec <- sort(unique(con$common_name))
met_spec <- sort(unique(met$common_name))
common_spec <- met_spec[met_spec %in% con_spec]

test <- con %>% filter(above_peak_temp == "Y")
length(unique(test$common_name))

# Filter data points at below optimum temperatures
met <- met %>% filter(above_peak_temp == "N")

# Count data first
nrow(con)
con <- con %>% filter(above_peak_temp == "N")
nrow(con)

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


con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))
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
gro <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

str(gro)

# Which species are in the growth data and any of the metabolism and consumption data sets?
gro_spec <- sort(unique(gro$common_name))
con_spec <- sort(unique(con$common_name))
met_spec <- sort(unique(met$common_name))

all_con_met_species <- c(con_spec, met_spec)

gro_spec[gro_spec %in% all_con_met_species]

# Count data
length(unique(gro$common_name))
nrow(gro)
gro2 <- gro %>% filter(above_peak_temp == "N")
nrow(gro2)

# How many temperatures? 
gro %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  as.data.frame()

# Average # of temperatures per species?
gro %>%
  group_by(species) %>% 
  summarise(n_unique_temp = length(unique(round(temp_c, digits = 0)))) %>% 
  ungroup() %>% 
  summarise(mean_n = mean(n_unique_temp))

# How many masses? 
gro %>%
  group_by(species) %>% 
  summarise(n_unique_mass = length(unique(round(mass_g, digits = 0)))) %>% 
  as.data.frame()

# Average # of masses per species?
gro %>%
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

# Total # of unique species
met_spec <- met %>% dplyr::select(species_ab) %>% distinct()
con_spec <- con %>% dplyr::select(species_ab) %>% distinct()
gro_spec <- gro %>% dplyr::select(species_ab) %>% distinct()
topt_spec <- topt %>% dplyr::select(species_ab) %>% distinct()

all_spec <- rbind(met_spec, con_spec, gro_spec, topt_spec)

length(unique(all_spec$species_ab))


# Total ============================================================================
# Total # of data points
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

gro <-
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

topt <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/topt_analysis.csv"))

length(unique(c(met$species_ab, con$species_ab, topt$species_ab)))


#nrow(met) + nrow(con) + nrow(gro)
nrow(met) + nrow(con) + nrow(gro)


# Raw data =========================================================================
con <- read_excel("data/consumption_data.xlsx")
met <- read_excel("data/metabolism_data.xlsx")
topt <- read_excel("data/growth_data_Topt.xlsx")
gro <- read_excel("data/growth_data.xlsx")


# Print unique studies
all_studies <- c(unique(c(met$reference, con$reference, topt$reference, gro$reference)))

sort(all_studies)

# Across rate
length(unique(c(met$reference, con$reference, topt$reference, gro$reference)))
length(unique(c(met$species, con$species, topt$species, gro$species)))

# Per rate
length(unique(met$species))
length(unique(con$species))
length(unique(topt$species))
length(unique(gro$species))

length(unique(met$species)) + length(unique(con$species)) + length(unique(topt$species)) + length(unique(gro$species))
length(unique(met$reference)) + length(unique(con$reference)) + length(unique(topt$reference)) + length(unique(gro$reference))

# metabolism and consumption references:
met_con_refs <- unique(c(met$reference, con$reference))

# Which of the growth references are in there?
# Use the datafile for analysis, since some are filtered from the raw data
met <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/met_analysis.csv"))

con <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/con_analysis.csv"))

gro <-
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_analysis.csv"))

topt <- 
  read.csv(text = getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/topt_analysis.csv"))

topt_spec <- unique(topt$species_ab)

topt_refs[topt_refs %in% unique(met$species_ab)]
topt_refs[topt_refs %in% unique(con$species_ab)]
topt_refs


