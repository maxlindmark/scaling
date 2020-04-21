#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to exemplify size-and temperature dependence of growth using Putter model
# 
# A. Load libraries
#
# B. Generate data and growth
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES & READ DATA ====================================================
rm(list = ls())

# Provide package names
library("dplyr")
library("tidyr")
library("tidylog")
library("ggplot2")
library("RColorBrewer")


# ==================================================================================
## Random activation energies
# Sample random activation energies
n = 50

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

pal <- brewer.pal(name = "Set1", n = 8)

p1 <- t_dat %>% filter(mass > 0) %>% 
  ggplot(., aes(mass, growth, color = factor(temp), group = group)) + 
  geom_line(size = 0.25, alpha = 0.2) + 
  theme_classic(base_size = 12) +
  labs(y = "Growth [dG/dt]", 
       x = "Mass [g]",
       color = "Temperature [C]") +
  scale_color_manual(values = pal[c(2, 1)]) +
  guides(color = FALSE) +
  geom_line(data = filter(sum_df, temp == 273.15 + 5), inherit.aes = FALSE, 
            aes(mass, mean_g), size = 0.8, linetype = 2, color = "white") +
  geom_line(data = filter(sum_df, temp == 273.15 + 10), inherit.aes = FALSE, 
            aes(mass, mean_g), size = 0.8, linetype = 2, color = "white") + 
  theme(aspect.ratio = 4/5,
        legend.text=element_text(size = 8)) +
  NULL
p1

p2 <- t_dat %>% filter(mass > 0) %>% 
  ggplot(., aes(log(mass), log(growth_s), color = factor(temp-273.15), group = group)) + 
  geom_line(size = 0.25, alpha = 0.2) + 
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

#ggsave("figures/putter_growth_rate.pdf", plot = last_plot(), scale = 1, width = 17, height = 17, units = "cm")

