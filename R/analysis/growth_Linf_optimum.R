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

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)

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
                              temp = seq(273.15 + 10, 273.15 + 30, 1
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

# Caluclate asymptotic mass (am^y – bm^z-  -->  (a/b)1/(z− y))
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








