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
          "RCurl",
          "RColorBrewer")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/packageInfo.R", ssl.verifypeer = FALSE)
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

# --- When I do a normal VBGE with arrhenius temperature, optimum is extremely dependent 
#     on mass! So for the purpose of illustrating, I will use temperature-dependence of 
#     Morita et al 2010, based on Atkinson 1994, and add temperature effects on exponents.
#     Will need to figure this out though. Perhaps I should scale constants separately?

dat <- data.frame(expand.grid(mass = seq(0, 500, 1),
                              temp = seq(273.15 + 10, 273.15 + 25, 0.01),
                              ca = c(0, -0.003)))
tref <- 273.15 + 10

dat$anab <- (3+0.7*(dat$temp - tref))*dat$mass^(2/3 + (dat$ca*(dat$temp - tref)))
dat$cata <- (0.01*(dat$temp - tref)^3)*dat$mass^1

# Calculate difference in energy gains and losses
dat$growth <- dat$anab - dat$cata

# Plot growth rates over temperature for specific sizes
pal <- viridis(n = 5)

p1 <- dat %>% 
  filter(growth > 0 & mass %in% c(100, 200, 300, 400, 500)) %>% 
  ggplot(., aes(temp-273.15, growth, color = factor(mass), linetype = factor(ca))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "Growth rate", x = "Temperature [C]") +
  scale_color_manual(values = pal[c(1:5)], 
                       name = "Mass [g]") +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.position = c(0.92, 0.65),
        legend.text=element_text(size = 11)) +
  NULL

# Trying to filter out optimum temperatures for each size...
p2 <- dat %>% 
  filter(growth > 0 & mass > 10) %>%
  group_by(factor(mass), factor(ca)) %>% # Find temp at max growth
  filter(growth == max(growth)) %>% 
  ungroup() %>% 
  ggplot(., aes(mass, temp-273.15, linetype = factor(ca))) + 
  geom_line(size = 2) + 
  theme_classic(base_size = 16) +
  scale_color_viridis(discrete = TRUE) +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = FALSE) +
  labs(y = expression(paste(T[opt], " [C]")), 
       x = "Mass [g]") +
  theme(aspect.ratio = 4/5) +
  NULL

p1 / p2
  
#ggsave("figs/t_opt_model.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm")


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
  theme_classic(base_size = 1) +
  theme(aspect.ratio = 1) +
  NULL

#**** Heatmap with species "exponent-change" as points =============================
dat <- data.frame(expand.grid(mass = seq(0, 1000, 10),
                              temp = c(273.15 + 10, 273.15 + 12),
                              ca   = seq(-0.015, 0.004, 0.0005),
                              cm   = seq(-0.015, 0.004, 0.0005)))

# Set exponents
bm <- 1
bi <- 2/3 
tref <- 273.15 + 10

# Scale parameters with size and temperature
dat$a_c <- 0.1 * exp(0.6 * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
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
m80 <- data.frame(ca = c(-0.0040, 0.0001),
                  cm = c(-0.0059, 0.0011))

cm_est <- -0.0024
ca_est <- -0.0020

ggplot(filter(dat, temp == 273.15 + 12), aes(x = cm, y = ca)) + 
  geom_raster(aes(fill = mass_cent, z = mass_cent), interpolate = F) +
  scale_fill_viridis(begin = 0, end = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = bquote('Change in Metabolism exponent'~(C^-1)), 
       y = bquote('Change in Cmax exponent'~(C^-1)),
       fill = "Proportion change\nin asymptotic size\n(10C-12C)") +
  geom_line(data = sdat, aes(cm, ca), size = 1, col = "gray95", linetype = 2, z = NULL) +
  geom_point(data = all_spec, aes(cm, ca), size = 5, col = "white") +
  # Add c-effects from all species-analysis
  # geom_segment(data = m80, aes(y = min(ca), yend = max(ca), x = cm_est, xend = cm_est), 
  #              size = 1, col = "red") +
  # geom_segment(data = m80, aes(y = ca_est, yend = ca_est, x = min(cm), xend = max(cm)), 
  #              size = 1, col = "red") +
  annotate("text", y = c(-0.012, -0.0012), x = c(-0.005, -0.010), size = 4.3, fontface = 3, 
           color = "white", label = c("max. size decreasing", "max. size increasing"),
           angle = 45) +
  theme_classic(base_size = 22) +
  theme(aspect.ratio = 1) +
  NULL

#ggsave("figs/winf_species_heatmap.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm")


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

