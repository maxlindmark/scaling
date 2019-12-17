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
# ==================================================================================
dat <- data.frame(expand.grid(mass = seq(0, 500, 1),
                              temp = seq(273.15 + 8, 273.15 + 16, 1)))
tref <- 273.15 + 12
a1 <- 0.1
a2 <- 0.01
b1 <- 0.6
b2 <- 0.77
Ea_con <- 0.63
Ea_met <- 0.57

dat$anab <- (a1*dat$mass^b1) * exp(Ea_con * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
dat$cata <- (a2*dat$mass^b2) * exp(Ea_met * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))

# Calculate difference in energy gains and losses
dat$growth <- dat$anab - dat$cata

# Plot growth rates over temperature for specific sizes
pal <- viridis(n = 5)

pal <- brewer.pal(name = "RdYlBu", n = length(unique(factor(dat$temp-273.15))))

p1 <- dat %>% 
  ggplot(., aes(mass, growth, color = factor(temp-273.15))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "dG/dt", 
       x = "Mass [g]",
       color =  "Temperature [C]") +
  #scale_color_brewer(palette = "RdYlBu") +
  scale_color_manual(values = rev(pal)) +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 8)) +
  NULL

p1

# Log scale of growth (allometry)
allom <- dat %>% 
  filter(temp %in% c(273.15 + 8, 273.15 + 16)) %>% 
  ggplot(., aes(log(mass), log(growth/mass), color = factor(temp-273.15))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "dG/dt", 
       x = "Mass [g]",
       color =  "Temperature [C]") +
  #scale_color_brewer(palette = "RdYlBu") +
  scale_color_manual(values = rev(pal)) +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 8)) +
  NULL

allom
#ggsave("figures/growth_rate.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm")

# ==================================================================================
## Mass instead on X
dat <- data.frame(expand.grid(mass = seq(0, 500, 50),
                              temp = seq(273.15 + 8, 273.15 + 16, 0.1)))
tref <- 273.15 + 12
a1 <- 0.1
a2 <- 0.01
b1 <- 0.6
b2 <- 0.77
Ea_con <- 0.63
Ea_met <- 0.57

dat$anab <- (a1*dat$mass^b1) * exp(Ea_con * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
dat$cata <- (a2*dat$mass^b2) * exp(Ea_met * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))

# Calculate difference in energy gains and losses
dat$growth <- dat$anab - dat$cata

p2 <- dat %>% 
  filter(mass > 0) %>% 
  ggplot(., aes(temp-273.15, growth, color = factor(mass))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "dG/dt", 
       x = "Temperature [C]",
       color =  "Mass [g]") +
  #scale_color_brewer(palette = "RdYlBu") +
  scale_color_viridis(discrete = TRUE, option = "cividis") +
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 8)) +
  NULL

p1/p2
#ggsave("figures/growth_rate.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm")

# Calculate the % increase in growth from 8C
dat$growth_cold <- rep(dat$growth[1:11], length(seq(273.15 + 8, 273.15 + 16, 0.1)))
dat$growth_rel <- dat$growth / dat$growth_cold

dat %>% 
  filter(mass > 0) %>% 
  ggplot(., aes(temp-273.15, growth_rel, color = factor(mass))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "dG/dt", 
       x = "Temperature [C]",
       color =  "Mass [g]") +
  scale_color_brewer(palette = "RdYlBu") +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 11)) +
  NULL


## Plot relative increase in growth rate

# ==================================================================================
## Different activation energies
dat <- data.frame(expand.grid(mass = seq(0, 500, 50),
                              temp = seq(273.15 + 8, 273.15 + 16, 0.1)))
tref <- 273.15 + 12
a1 <- 0.1
a2 <- 0.01
b1 <- 0.6
b2 <- 0.77
Ea_con <- 0.5
Ea_met <- 0.7

dat$anab <- (a1*dat$mass^b1) * exp(Ea_con * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
dat$cata <- (a2*dat$mass^b2) * exp(Ea_met * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))

# Calculate difference in energy gains and losses
dat$growth <- dat$anab - dat$cata

dat %>% 
  filter(mass > 0) %>% 
  ggplot(., aes(temp-273.15, growth, color = factor(mass))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "dG/dt", 
       x = "Temperature [C]",
       color =  "Mass [g]") +
  scale_color_brewer(palette = "RdYlBu") +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 11)) +
  NULL

# Calculate the % increase in growth from 8C
dat$growth_cold <- rep(dat$growth[1:11], length(seq(273.15 + 8, 273.15 + 16, 0.1)))
dat$growth_rel <- dat$growth / dat$growth_cold

dat %>% 
  filter(mass > 0) %>% 
  ggplot(., aes(temp-273.15, growth_rel, color = factor(mass))) + 
  geom_line(size = 1.5) + 
  theme_classic(base_size = 16) +
  labs(y = "dG/dt", 
       x = "Temperature [C]",
       color =  "Mass [g]") +
  scale_color_brewer(palette = "RdYlBu") +
  scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 11)) +
  NULL
d