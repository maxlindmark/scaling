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
# # ==================================================================================
# dat <- data.frame(expand.grid(mass = seq(0, 500, 1),
#                               temp = seq(273.15 + 8, 273.15 + 16, 1)))
# tref <- 273.15 + 12
# a1 <- 0.1
# a2 <- 0.01
# b1 <- 0.6
# b2 <- 0.77
# Ea_con <- 0.63
# Ea_met <- 0.57
# 
# dat$anab <- (a1*dat$mass^b1) * exp(Ea_con * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
# dat$cata <- (a2*dat$mass^b2) * exp(Ea_met * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
# 
# # Calculate difference in energy gains and losses
# dat$growth <- dat$anab - dat$cata
# 
# # Plot growth rates over temperature for specific sizes
# pal <- viridis(n = 5)
# 
# pal <- brewer.pal(name = "RdYlBu", n = length(unique(factor(dat$temp-273.15))))
# 
# p1 <- dat %>% 
#   ggplot(., aes(mass, growth, color = factor(temp-273.15))) + 
#   geom_line(size = 1.5) + 
#   theme_classic(base_size = 16) +
#   labs(y = "dG/dt", 
#        x = "Mass [g]",
#        color =  "Temperature [C]") +
#   #scale_color_brewer(palette = "RdYlBu") +
#   scale_color_manual(values = rev(pal)) +
#   scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
#   guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
#   theme(aspect.ratio = 4/5,
#         legend.text=element_text(size = 8)) +
#   NULL
# 
# p1
# 
# # Log scale of growth (allometry)
# allom <- dat %>% 
#   filter(temp %in% c(273.15 + 8, 273.15 + 16)) %>% 
#   ggplot(., aes(log(mass), log(growth/mass), color = factor(temp-273.15))) + 
#   geom_line(size = 1.5) + 
#   theme_classic(base_size = 16) +
#   labs(y = "dG/dt", 
#        x = "Mass [g]",
#        color =  "Temperature [C]") +
#   #scale_color_brewer(palette = "RdYlBu") +
#   #scale_color_manual(values = rev(pal)) +
#   scale_color_manual(values = rev(pal[c(1, 9)])) +
#   scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
#   guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
#   theme(aspect.ratio = 4/5,
#         legend.text=element_text(size = 8)) +
#   NULL
# 
# allom
# #ggsave("figures/growth_rate.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm")
# 
# # ==================================================================================
# ## Temp instead on X
# dat <- data.frame(expand.grid(mass = seq(0, 500, 50),
#                               temp = seq(273.15 + 8, 273.15 + 16, 0.1)))
# tref <- 273.15 + 12
# a1 <- 0.1
# a2 <- 0.01
# b1 <- 0.6
# b2 <- 0.77
# Ea_con <- 0.63
# Ea_met <- 0.57
# 
# dat$anab <- (a1*dat$mass^b1) * exp(Ea_con * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
# dat$cata <- (a2*dat$mass^b2) * exp(Ea_met * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
# 
# # Calculate difference in energy gains and losses
# dat$growth <- dat$anab - dat$cata
# 
# p2 <- dat %>% 
#   filter(mass > 0) %>% 
#   ggplot(., aes(temp-273.15, growth, color = factor(mass))) + 
#   geom_line(size = 1.5) + 
#   theme_classic(base_size = 16) +
#   labs(y = "dG/dt", 
#        x = "Temperature [C]",
#        color =  "Mass [g]") +
#   #scale_color_brewer(palette = "RdYlBu") +
#   scale_color_viridis(discrete = TRUE, option = "cividis") +
#   guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
#   theme(aspect.ratio = 4/5,
#         legend.text=element_text(size = 8)) +
#   NULL
# 
# p1/p2
# #ggsave("figures/growth_rate.pdf", plot = last_plot(), scale = 1, width = 23, height = 23, units = "cm")
# 
# # Calculate the % increase in growth from 8C
# dat$growth_cold <- rep(dat$growth[1:11], length(seq(273.15 + 8, 273.15 + 16, 0.1)))
# dat$growth_rel <- dat$growth / dat$growth_cold
# 
# dat %>% 
#   filter(mass > 0) %>% 
#   ggplot(., aes(temp-273.15, growth_rel, color = factor(mass))) + 
#   geom_line(size = 1.5) + 
#   theme_classic(base_size = 16) +
#   labs(y = "dG/dt", 
#        x = "Temperature [C]",
#        color =  "Mass [g]") +
#   scale_color_brewer(palette = "RdYlBu") +
#   scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
#   guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
#   theme(aspect.ratio = 4/5,
#         legend.text=element_text(size = 11)) +
#   NULL


# ==================================================================================
## Random activation energies
# Sample random activation energies
n = 100

# Metabolic rate
met_u <- 0.62
met_sd <- 0.03
met_q <- qnorm(c(0.025, 0.975), met_u, met_sd)
met <- rnorm(n = n, met_u, met_sd)

# Maximum consumption rate
int_u <- 0.69
int_sd <- 0.08
int_q <- qnorm(c(0.025, 0.975), int_u, int_sd)
int <- rnorm(n = n, int_u, int_sd)

Ea_met <- met
Ea_con <- int

# Create data.frame
# dat <- data.frame(expand.grid(mass = seq(0, 500, 1),
#                               temp = c(273.15 + 8, 273.15 + 10, 273.15 + 12)))

tref <- 273.15 + 7.5
a1 <- 0.1
a2 <- 0.01
b1 <- 0.6
b2 <- 0.77

t <- c()
data_list <- list()

for(i in 1:n) {
  
  t <- data.frame(expand.grid(mass = seq(0, 500, 1),
                              temp = c(273.15 + 5, 273.15 + 10)))
  
  t$anab <- (a1*t$mass^b1) * exp(Ea_con[i] * ((t$temp - tref) / ((8.617332e-05) * t$temp * tref)))
  t$cata <- (a2*t$mass^b2) * exp(Ea_met[i] * ((t$temp - tref) / ((8.617332e-05) * t$temp * tref)))
  t$Ea_met <- Ea_met[i]
  t$Ea_con <- Ea_con[i]
  t$sim <- paste(Ea_met[i], Ea_con[i], sep = ",")
    
  data_list[[i]] <- t
    
}

t_dat <- dplyr::bind_rows(data_list)

# Calculate difference in energy gains and losses
t_dat$growth <- t_dat$anab - t_dat$cata
t_dat$growth_s <- t_dat$growth / t_dat$mass

t_dat$group <- paste(t_dat$temp, t_dat$sim)

sum_df <- t_dat %>% 
  group_by(mass, temp) %>% 
  summarize(mean_g = mean(growth),
            mean_g_s = mean(growth_s))

head(sum_df)

#pal <- brewer.pal(name = "RdYlBu", n = length(unique(factor(dat$temp-273.15))))
pal <- brewer.pal(name = "Set1", n = 8)

p1 <- t_dat %>% filter(mass > 0) %>% 
  ggplot(., aes(mass, growth, color = factor(temp), group = group)) + 
  geom_line(size = 0.75, alpha = 0.1) + 
  theme_classic(base_size = 12) +
  labs(y = "Growth [dG/dt]", 
       x = "Mass [g]",
       color = "Temperature [C]") +
  scale_color_manual(values = pal[c(2, 1)]) +
  #guides(color = guide_legend(override.aes = list(alpha = 1))) +
  guides(color = FALSE) +
  geom_line(data = filter(sum_df, temp == 273.15 + 5), inherit.aes = FALSE, 
            aes(mass, mean_g), size = 0.8, linetype = 2, color = "white") +
  # geom_line(data = filter(sum_df, temp == 283.15), inherit.aes = FALSE, 
  #           aes(mass, mean_g), size = 1, linetype = 2, color = pal[3]) + 
  geom_line(data = filter(sum_df, temp == 273.15 + 10), inherit.aes = FALSE, 
            aes(mass, mean_g), size = 0.8, linetype = 2, color = "white") + 
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 8)) +
  NULL
p1

p2 <- t_dat %>% filter(mass > 0) %>% 
  ggplot(., aes(log(mass), log(growth_s), color = factor(temp-273.15), group = group)) + 
  geom_line(size = 0.75, alpha = 0.1) + 
  theme_classic(base_size = 12) +
  labs(y = "ln(specific growth [dG/dt])", 
       x = "ln(Mass [g])",
       #color = "Temperature [C]"
       color = expression(paste("Temperature [", degree*C, "]"))) +
  scale_color_manual(values = pal[c(2, 1)]) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  geom_line(data = filter(sum_df, temp == 273.15 + 5), inherit.aes = FALSE, 
            aes(log(mass), log(mean_g_s)), size = 0.8, linetype = 2, color = "white") +
  # geom_line(data = filter(sum_df, temp == 283.15), inherit.aes = FALSE, 
  #           aes(log(mass), log(mean_g_s)), size = 1, linetype = 2, color = pal[3]) + 
  geom_line(data = filter(sum_df, temp == 273.15 + 10), inherit.aes = FALSE, 
            aes(log(mass), log(mean_g_s)), size = 0.8, linetype = 2, color = "white") + 
  theme(aspect.ratio = 4/5,
        #legend.text=element_text(size = 8),
        legend.position = "bottom") +
  NULL
p2

p1/p2

#ggsave("figures/growth_rate.pdf", plot = last_plot(), scale = 1, width = 17, height = 17, units = "cm")



# # ==================================================================================
# ## Different activation energies
# dat <- data.frame(expand.grid(mass = seq(0, 500, 50),
#                               temp = seq(273.15 + 8, 273.15 + 16, 0.1)))
# tref <- 273.15 + 12
# a1 <- 0.1
# a2 <- 0.01
# b1 <- 0.6
# b2 <- 0.77
# Ea_con <- 0.5
# Ea_met <- 0.7
# 
# dat$anab <- (a1*dat$mass^b1) * exp(Ea_con * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
# dat$cata <- (a2*dat$mass^b2) * exp(Ea_met * ((dat$temp - tref) / ((8.617332e-05) * dat$temp * tref)))
# 
# # Calculate difference in energy gains and losses
# dat$growth <- dat$anab - dat$cata
# 
# dat %>% 
#   filter(mass > 0) %>% 
#   ggplot(., aes(temp-273.15, growth, color = factor(mass))) + 
#   geom_line(size = 1.5) + 
#   theme_classic(base_size = 16) +
#   labs(y = "dG/dt", 
#        x = "Temperature [C]",
#        color =  "Mass [g]") +
#   scale_color_brewer(palette = "RdYlBu") +
#   scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
#   guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
#   theme(aspect.ratio = 4/5,
#         legend.text=element_text(size = 11)) +
#   NULL
# 
# # Calculate the % increase in growth from 8C
# dat$growth_cold <- rep(dat$growth[1:11], length(seq(273.15 + 8, 273.15 + 16, 0.1)))
# dat$growth_rel <- dat$growth / dat$growth_cold
# 
# dat %>% 
#   filter(mass > 0) %>% 
#   ggplot(., aes(temp-273.15, growth_rel, color = factor(mass))) + 
#   geom_line(size = 1.5) + 
#   theme_classic(base_size = 16) +
#   labs(y = "dG/dt", 
#        x = "Temperature [C]",
#        color =  "Mass [g]") +
#   scale_color_brewer(palette = "RdYlBu") +
#   scale_linetype_manual(values = c("twodash", "solid"), name = "c") + 
#   guides(linetype = guide_legend(override.aes = list(size = 1.3))) +
#   theme(aspect.ratio = 4/5,
#         legend.text=element_text(size = 11)) +
#   NULL
