#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - Explore body growth rate data
# 
# A. Load libraries & read data
#
# B. Explore growth data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#======== A. LOAD LIBRARIES & READ DATA ============================================
rm(list = ls())

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

# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

dat <- read_excel("data/growth_data.xlsx")

glimpse(dat)

cols = c(1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)


#======== B. EXPLORE DATA ==========================================================
#==** General ======================================================================
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
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

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


#==** Growth rate ==================================================================
# Normalize size by creating new column with size relative to max. In most cases this will be geometric mean, but can also be size class. I will then check if there are any NA's (Walleye) that doesn't have either size and make a data.

# Mass for analysis
dat$mass <- dat$geom_mean_mass_g

# Replace NA with -9...
dat$mass[is.na(dat$mass)] <- -9

dat$mass <- ifelse(dat$mass == -9,
                   dat$mid_mass_sizecl_g,
                   dat$mass)

dat$mass

dat$mass[is.na(dat$mass)] <- -9

# This is a single size-species. Will not bother finding a temperature for this as the main exploration here is on intraspecific growth
test <- dat %>% filter(mass == -9)

dat <- dat %>% filter(mass > 0)

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass / dat$w_max_published_g

#====**** All data (also n=1 species) ==============================================
# Non-normalized mass
ggplot(dat, aes(mass, opt_temp_c)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  NULL

# Log10 mass
ggplot(dat, aes(log10(mass), opt_temp_c)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  NULL

# Normalized mass
ggplot(dat, aes(mass_norm, opt_temp_c)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  NULL

# Log10 normalized mass
ggplot(dat, aes(log10(mass_norm), opt_temp_c)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  NULL

# Log10 normalized mass color species
ggplot(dat, aes(log10(mass_norm), opt_temp_c, color = common_name)) + 
  geom_point(size = 5, alpha = 0.8) +
  geom_line(size = 3, alpha = 0.8) +
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

#====**** Intraspecific data (only n>1 species) ====================================
s_dat <- dat %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Split window by species
ggplot(s_dat, aes(log10(mass_norm), opt_temp_c, color = common_name)) + 
  geom_line(size = 4, alpha = 0.4) +
  geom_point(size = 4, alpha = 1) +
  facet_wrap(~common_name, scales = "free") +
  theme_bw(base_size = 15) +
  guides(color = FALSE) +
  scale_color_viridis(discrete = TRUE) +
  NULL

# Log10 normalized mass color species
ggplot(s_dat, aes(log10(mass_norm), opt_temp_c, color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  geom_line(size = 3, alpha = 0.2) +
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

ggplot(s_dat, aes(log10(mass_norm), opt_temp_c, color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

# I should normalize temp, because the previous plots assume that a certain size relative to max will have a certain optimum. In other words, the optimum is shared across species! Here I do it by experimental temperature, because what else? I could do environmental temperature like for consumption and metabolism, but this is a more relevant reference.

# Calculate mean optimum temperature within species
s_dat$mean_opt_temp_c <- ave(s_dat$opt_temp_c, s_dat$common_name)

# Center each size class' optimum temperature by mean optimum temperature for that species
s_dat$opt_temp_c_ct <- s_dat$opt_temp_c - s_dat$mean_opt_temp_c

# Define more readable palette
nb.cols <- length(unique(s_dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat, aes(log10(mass_norm), opt_temp_c_ct, 
                  color = common_name)) + 
  geom_point(size = 5) +
  #geom_line(size = 3, alpha = 0.4) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = mycolors) +  
  NULL

summary(lm(s_dat$opt_temp_c_ct ~ log10(s_dat$mass_norm)))


# TO THINK ABOUT:
# what is the rationale for not mixing species from different studies here? For Cmax we don't, but that's because they are so tricky to measure.. growth should be easier to measure. And by grouping species and sharing information, we do say the response variables are comparable. Note I also have much fewer species with dublicates in Cmax. 

#====**** Conclusion ===============================================================
# Do intraspecific analysis, using the same filters as above. Fit LMM to start with to account for non-independence? Here' I'm not interested in the actual slopes of each species but only the overall effect, so should be straightforward


#====**** Extra ====================================================================
# Growth rate at optimum over size 
ggplot(s_dat, aes(log10(mass_norm), `growth_rate_%day-1`, 
                  color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm", se = FALSE, size = 2, alpha = 0.2)+ 
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

ggplot(s_dat, aes(log(mass), log(`growth_rate_%day-1`), 
                  color = common_name)) + 
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm", se = FALSE, size = 2, alpha = 0.2)+ 
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL

summary(lm(log(s_dat$growth_rate) ~ log(s_dat$mass)))
# Mass specific growth size scaling exponen: is -0.34 + 1
# No one has done this for intraspecific growth and not using only optimum temperatures!
# See Barneche's new paper for references on this...

# Growth rate over temp
ggplot(s_dat, aes(opt_temp_c_ct, growth_rate, 
                  color = common_name)) + 
  geom_point(size = 7, alpha = 0.7) +
  stat_smooth(method = "lm", se = FALSE, size = 2, alpha = 0.2)+ 
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = TRUE) +
  NULL








