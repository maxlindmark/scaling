#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.27: Max Lindmark
#
# - The following code plots VBGE growth and L_inf with standard arrhenius temp dep
#
# A. Morita
#
# B. Arrhenius
#
# C. Summary / Conclusions
#      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Will clean this script later. Figures here are more conceptual anyway

# A. LOAD LIBRARIES & READ DATA ====================================================
rm(list = ls())

# Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "ggplot2",
          "viridis",
          "RCurl")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# pkg_info(pkgs)

# package loadedversion
# 1   dplyr         0.8.1
# 2 ggplot2         3.1.1
# 3   RCurl     1.95-4.12
# 4 tidylog         0.1.0
# 5   tidyr         0.8.3
# 6 viridis         0.5.1


# B. Growth ========================================================================
#** Simulate optimum ===============================================================
# In contrast to asymptotic size, we don't have an expression for how T_opt depend
# on allometric functions, so here we need to simulate.
# Create data frame (notice we now vary the activation energy for anabolism relative 
# to catabolism with parameter E)
dat <- data.frame(expand.grid(mass = seq(0, 1000, 10),
                              temp = seq(273.15 + 10, 273.15 + 30, 1),
                              E    = c(0.5, 0.6, 0.7),
                              ca   = c(-0.003, 0)))


tref <- 273.15 + 20

# Scale rates with size and temperature
dat$anab <- (0.1*dat$mass^((2/3)+(dat$ca*(dat$temp-tref)))) * exp((dat$E * (dat$temp - tref)) / ((8.617332e-05) * dat$temp * tref))

dat$cata <- (0.01*dat$mass^1) * exp((0.6 * (dat$temp - tref)) / ((8.617332e-05) * dat$temp * tref))

# Calculate difference in energy gains and losses
dat$growth <- dat$anab - dat$cata

# No c-effect
dat %>% 
  filter(ca == 0 & growth > 0) %>% 
  ggplot(., aes(temp, growth, color = mass, group = factor(mass))) + 
  geom_line(size = 3, alpha = 0.5) + 
  facet_wrap(~E, scales = "free_y") +
  theme_classic(base_size = 18) +
  scale_color_viridis(discrete = F) +
  NULL

# Negative c-effect also
dat$E_p <- as.factor(dat$E)
levels(dat$E_p) <- c("E=0.5", "E=0.6", "E=0.7")

dat$ca_p <- as.factor(dat$ca)
levels(dat$ca_p) <- c("ca=-0.003", "ca=0")

dat %>% 
  filter(growth > 0 & mass > 200 & mass < 800) %>% 
  mutate(tempc = temp-273.15) %>% 
  ggplot(., aes(tempc, growth, color = mass, group = factor(mass))) + 
  geom_line(size = 3, alpha = 0.5) + 
  facet_wrap(ca_p ~ E_p, scales = "free_y") +
  labs(x = "Temperature [C]", y = "Growth rate") +
  theme_classic(base_size = 15) +
  scale_color_viridis(discrete = F) +
  NULL

#ggsave("figs/winf_species_heatmap.pdf", plot = last_plot(), scale = 1, width = 22, height = 22, units = "cm")

# Trying to filter out optimum temperatures for each size...
# STRONG DECLINES IN OPT WITH THESE PARAMETERS
opt_dat <- dat %>% 
  filter(growth > 0) %>% # Positive growth only
  group_by(factor(mass), factor(E), factor(ca)) %>% # Find temp at max growth
  filter(growth == max(growth) & temp > tref & temp < (273.15 + 30)) # Filter only highest growth (optimum) and remove size that have highest max at ref temp (optimum is lower!)
  
ggplot(opt_dat, aes(mass, temp, color = factor(ca), shape = factor(E))) + 
  geom_point(size = 3, alpha = 0.5) +
  theme_classic(base_size = 18) +
  labs(x = "Mass [g]", y = "Optimum temperature") +
  scale_color_viridis(discrete = T) +
  NULL

ggplot(opt_dat, aes(x = mass, y = E)) + 
  geom_raster(aes(fill = temp, z = temp), interpolate = F) +
  scale_fill_viridis(begin = 0, end = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~ca, nrow = 1) +
  labs(x = "mass", 
       y = "E",
       fill = "Optimum growth temperature") +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1) +
  NULL


#** Asymptotic size ================================================================
asymmass <- data.frame(expand.grid(temp = seq(273.15 + 10, 273.15 + 30, 1),
                                   Eratio = seq(0.3, 0.9, 0.01),
                                   ca = seq(-0.005, 0, 0.001)))

# Scale allometric constant
asymmass$a <- 0.1 * exp((asymmass$Eratio * (asymmass$temp - tref)) / ((8.617332e-05) * asymmass$temp * tref))
asymmass$b <- 0.01 * exp((0.6 * (asymmass$temp - tref)) / ((8.617332e-05) * asymmass$temp * tref))

# Scale C-max exponent
asymmass$y <- (2/3)+(asymmass$ca*(asymmass$temp-tref))

# Caluclate asymptotic mass (am^y – bm^z-  -->  (a/b)^1/(z−y))
asymmass$wInf <- (asymmass$a/asymmass$b) ^ 1/(1−asymmass$y)

head(asymmass)

# Calculate ref size:
refsize <- (0.1/0.01)^1/(1−(2/3))

asymmass$mnorm <- asymmass$wInf / refsize

asymmass$Eratio_f <- as.factor(asymmass$Eratio)

asymmass %>% filter(Eratio_f %in% c("0.4", "0.6", "0.8")) %>%  
  ggplot(., aes(temp, mnorm, color = ca, group = ca)) +
  geom_line(size = 3, alpha = 0.8) + 
  facet_wrap(~Eratio, scales = "free_y") +
  theme_classic(base_size = 18) +
  scale_color_viridis(discrete = F) +
  NULL

asymmass %>% filter(Eratio_f %in% c("0.4", "0.6", "0.8")) %>%  
ggplot(., aes(x = temp, y = ca)) + 
  geom_raster(aes(fill = mnorm, z = wInf), interpolate = F) +
  scale_fill_viridis(begin = 0, end = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~Eratio, nrow = 1) +
  labs(x = bquote('Temperature'~(C^-1)), 
       y = bquote('Change in Cmax Exponent'~(C^-1)),
       fill = "Proportion change\nin asymptotic size") +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1) +
  NULL

# Heatmap with species "exponent-change" as points
dat <- data.frame(expand.grid(mass = seq(0, 1000, 10),
                              temp = c(273.15 + 10, 273.15 + 12),
                              ca   = seq(-0.015, 0.004, 0.0005),
                              cm   = seq(-0.015, 0.004, 0.0005)))

# Set exponents
bm <- 1
bi <- 2/3 
tref <- 273.15 + 10

# Scale parameters with size and temperature
# * NOTE I should use my parameter values for temperature scaling once I have good params. Now I have generic ones
dat$a_c <- 0.1 * exp(0.4 * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
dat$b_c <- bi + dat$ca * (dat$temp - (273.15 + 10))

dat$a_m <- 0.01 * exp(0.6 * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
dat$b_m <- bm + dat$cm * (dat$temp - (273.15 + 10))

# Caluclate asymptotic mass (am^y – bm^z-  -->  (a/b)^1/(z−y))
dat$w_inf <- (dat$a_c/dat$a_m)^(1/(dat$b_m−dat$b_c))


# Add species combinations of cm and ca to be plotted on top of heatmap
# Roach
roach_m <- data.frame(bm = c(0.8, 0.795, 0.684, 0.729, 0.816), 
                      temp = c(5, 10, 15, 20, 23))

summary(lm(roach_m$bm ~ roach_m$temp))

roach_c <- data.frame(bc = c(0.92, 0.83, 0.81, 0.76),
                      temp = c(5, 10, 15, 20))

summary(lm(roach_c$bc ~ roach_c$temp))

roach <- data.frame(species = "roach",
                    cm = -0.001278,
                    cm_upr = -0.001278 + (0.004374 * 1.96),
                    cm_lwr = -0.001278 - (0.004374 * 1.96),
                    ca = -0.001278, 
                    ca_upr = -0.001278 + (0.001897 * 1.96),
                    ca_lwr = -0.001278 - (0.001897 * 1.96))

# Bream
bream_m <- data.frame(bm = c(0.651, 0.717, 0.726, 0.679, 0.692), 
                      temp = c(5, 10, 15, 20, 23))

summary(lm(bream_m$bm ~ bream_m$temp))

bream_c <- data.frame(bc = c(0.96, 0.892, 0.934, 0.742),
                      temp = c(5, 10, 15, 20))

summary(lm(bream_c$bc ~ bream_c$temp))

bream <- data.frame(species = "bream",
                    cm = 0.001041,
                    cm_upr = 0.001041 + (0.002304 * 1.96),
                    cm_lwr = 0.001041 - (0.002304 * 1.96),
                    ca = -0.012240,
                    ca_upr = -0.012240 + (0.006248 * 1.96),
                    ca_lwr = -0.012240 - (0.006248 * 1.96))

# Plaice
plaice_m <- data.frame(bm = c(0.790, 0.774, 0.799, 0.782, 0.778, 0.784), 
                       temp = c(2, 6, 10, 14, 18, 22))

summary(lm(plaice_m$bm ~ plaice_m$temp))

plaice_c <- data.frame(bc = c(0.683, 0.92, 0.761, 0.711, 0.679, 0.701),
                       temp = c(2, 6, 10, 14, 18, 22))

summary(lm(plaice_c$bc ~ plaice_c$temp))

plaice <- data.frame(species = "plaice",
                     cm = -0.000250,
                     cm_upr = -0.000250 + (0.000584 * 1.96),
                     cm_lwr = -0.000250 - (0.000584 * 1.96),
                     ca = -0.004879,
                     ca_upr = -0.004879 + (0.005628 * 1.96),
                     ca_lwr = -0.004879 - (0.005628 * 1.96))

# Flounder
flounder_m <- data.frame(bm = c(0.739, 0.757, 0.799, 0.775, 0.799, 0.795), 
                         temp = c(2, 6, 10, 14, 18, 22))

summary(lm(flounder_m$bm ~ flounder_m$temp))

flounder_c <- data.frame(bc = c(0.79, 0.859, 0.776, 0.812, 0.767, 0.783),
                         temp = c(2, 6, 10, 14, 18, 22))

summary(lm(flounder_c$bc ~ flounder_c$temp))

flounder <- data.frame(species = "flounder",
                       cm = 0.0027286,
                       cm_upr = 0.0027286 + (0.0009704 * 1.96),
                       cm_lwr = 0.0027286 - (0.0009704 * 1.96),
                       ca = -0.001964,
                       ca_upr = -0.001964 + (0.002020 * 1.96),
                       ca_lwr = -0.001964 - (0.002020 * 1.96))


# Add all species together
all_spec <- rbind(bream, roach, flounder, plaice)

dat$mass_cent <- dat$w_inf/1000

# Zero-change isocline
sdat <- subset(dat, temp == 285.15 & mass_cent < 1.01 & mass_cent > 0.99)

# Plot heatmap
pal <- viridis_pal()

# Data for 80% credible interval 
m80 <- data.frame(ca = c(-0.00677, -0.00029),
                  cm = c(-0.00619, 0.00406))

ggplot(filter(dat, temp == 273.15 + 12), aes(x = cm, y = ca)) + 
  geom_raster(aes(fill = mass_cent, z = mass_cent), interpolate = F) +
  scale_fill_viridis(begin = 0, end = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = bquote('Change in Metabolism Exponent'~(C^-1)), 
       y = bquote('Change in Imax Exponent'~(C^-1)),
       fill = "Proportion change\nin asymptotic size\n(10C-12C)") +
  geom_line(data = sdat, aes(cm, ca), size = 1, col = "gray95", linetype = 2, z = NULL) +
  #geom_abline(intercept = 0, slope = 1, col = "gray95", linetype = 2) +
  geom_point(data = all_spec, aes(cm, ca), size = 5, col = "white") +
  # Add c-effects from all species-analysis
  geom_segment(data = m80, aes(y = min(ca), yend = max(ca), x = -0.00106, xend = -0.00106), 
               size = 1, col = "red") +
  geom_segment(data = m80, aes(y = -0.00345, yend = -0.00345, x = min(cm), xend = max(cm)), 
               size = 1, col = "red") +
  annotate("text", y = c(-0.0035, -0.01), x = c(-0.011, -0.011), size = 4.3, fontface = 3, 
           color = "white", label = c("max. size increasing", "max. size decreasing"),
           angle = 45) +
  theme_classic(base_size = 22) +
  theme(aspect.ratio = 1) +
  NULL

#ggsave("figs/winf_species_heatmap.pdf", plot = last_plot(), scale = 1, width = 22, height = 22, units = "cm")

#** Summary ========================================================================

# Morita et.al showed that optimum decreases with size. They don't show how an optimum emerges, and we can safely assume growth eventually reaches an optimum. However, it is important to show how and when growth has an optimum, not least in the temperature range where people use arrhenius-scaling.

# These results show that optimum declines together with asymptotic size when:
# - Exponential
# 1) when exponent-difference decline (NEW!)
# 2) the activation energy is lower for intake (can't evalute yet, lit. says both ways)

# - Unimodal
# 3) Cmax is unimodal with lower T_opt than metabolism
# 4) Same, but stronger effect if Cmax opt is size-dependent
# (both partly new, don't really know!)

# To do:

# - do i make a heatmsp of topt as well? 
# - fit activation energy, update the story above ()
# - fit garcia garcia model and predict when optimum difference is highest?

# The T_opt figure is an independent evaluation of prediction (irrespective if how it emerges!!!)
# The pred-vs-obs is trickier and more conceptual at this stage. The main point is probably just to confirm we can use those rates... getting it comparable with the cmax graph is not straightforward.

