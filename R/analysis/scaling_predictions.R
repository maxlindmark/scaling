#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.12.02: Max Lindmark
#
# - Code to fit hierarchical model of maximum consumption rate as a function of 
# temperature with different group-effects and compare DIC and WAIC
# 
# A. Load libraries
#
# B. Read data
#
# C. Model selection
#
# D. Model validation
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#****** IN THIS SCRIPT I WILL REFIT BOTH MODELS BUT WITH PREDICTIONS INCLUDED! T

# THAT MIGHT MEAN I WILL HAVE TO REDEFINE THEME... STATE CLEARLY WHY AND HOW THAT
# WILL RESULTS IN NOT USING THE EXACT MODEL THAT I DEFINED IN THE MODEL SELECTION



# E. PLOT MODEL ====================================================================
js = jags.samples(jm, 
                  variable.names = c("mu_b0", "mu_b1", "mu_b2", "b3"), #, "pred_warm", "pred_cold", "sigma"), 
                  n.iter = samples, 
                  thin = n.thin)

# # Warm temp:
# pred_warm <- summary(js$pred_warm, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat
# 
# # Create data frame for the predictions
# pred_warm_df <- data.frame(lwr_95 = pred_warm[1, ],
#                            lwr_80 = pred_warm[2, ],
#                            median = pred_warm[3, ],
#                            upr_80 = pred_warm[4, ],
#                            upr_95 = pred_warm[5, ],
#                            mass = mass_pred,
#                            temp = -1.5)
# 
# # Cold temp:
# pred_cold <- summary(js$pred_cold, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# # Create data frame for the predictions
# pred_cold_df <- data.frame(lwr_95 = pred_cold[1, ],
#                            lwr_80 = pred_cold[2, ],
#                            median = pred_cold[3, ],
#                            upr_80 = pred_cold[4, ],
#                            upr_95 = pred_cold[5, ],
#                            mass = mass_pred,
#                            temp = 1.5)

# Plot data and predictions with 95% credible interval (at each x, plot as ribbon)
pal <- colorRampPalette(brewer.pal(8, "Dark2"))(8)

p0 <- ggplot(dat, aes(log_mass_norm_ct, log10(y))) +
  # geom_ribbon(data = pred_warm_df, aes(x = mass, ymin = lwr_95, ymax = upr_95), 
  #             size = 2, alpha = 0.2, inherit.aes = FALSE, fill = pal[1]) +
  # geom_ribbon(data = pred_warm_df, aes(x = mass, ymin = lwr_80, ymax = upr_80), 
  #             size = 2, alpha = 0.25, inherit.aes = FALSE, fill = pal[1]) +
  # geom_line(data = pred_warm_df, aes(mass, median), 
  #           size = 1, alpha = 0.8, color = pal[1]) +
  # geom_ribbon(data = pred_cold_df, aes(x = mass, ymin = lwr_95, ymax = upr_95), 
  #             size = 2, alpha = 0.2, inherit.aes = FALSE, fill = pal[3]) +
  # geom_ribbon(data = pred_cold_df, aes(x = mass, ymin = lwr_80, ymax = upr_80), 
  #             size = 2, alpha = 0.25, inherit.aes = FALSE, fill = pal[3]) +
  # geom_line(data = pred_cold_df, aes(mass, median), 
#           size = 1, alpha = 0.8, color = pal[3]) +
geom_point(size = 2.8, shape = 21, 
           alpha = 0.8, color = "white", fill = "grey40") +
  theme_classic(base_size = 11) + 
  labs(x = "ln(standardized mass)",
       y = "ln(growth rate)") +
  annotate("text", -Inf, Inf, label = "A", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  NULL

# Posterior of parameters
color_scheme_set("gray")
sum_dat <- data.frame(summary(cs2)[1])

# Mass-coefficient
p1 <- cs2 %>% 
  mcmc_dens(pars = "mu_b1") +
  theme_classic(base_size = 11) + 
  geom_vline(xintercept = sum_dat[1, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "B", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Mass-exponent") +
  NULL

# Temperature-coefficient
p2 <- cs2 %>% 
  mcmc_dens(pars = "mu_b2") +
  theme_classic(base_size = 11) + 
  #geom_vline(xintercept = sum_dat[3, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "C", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "Activation energy") +
  NULL

# Mass-temperature interaction
p3 <- cs2 %>% 
  mcmc_dens(pars = "b3") +
  theme_classic(base_size = 11) + 
  geom_vline(xintercept = 0, color = "red", size = 0.6, linetype = 2) +
  geom_vline(xintercept = sum_dat[5, 1], color = "white", size = 0.6, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", -Inf, Inf, label = "D", size = 4, 
           fontface = "bold", hjust = -0.5, vjust = 1.3) +
  labs(x = "M*T interaction") +
  NULL

# Plot all together
p1 / (p2 + p3 + p4) + plot_layout(ncol = 1, heights = c(2.5, 1, 1))

#ggsave("figs/growth_interaction.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")

# Calculate the proportion of the posterior that is less than zero
1-ecdf(js$b3)(0) # We are 84% certain the slope is larger than 0

summary(js$b3, quantile, c(0.025, 0.1, .5, 0.9, 0.975))$stat

# [,1]
# 2.5%  -0.9469043
# 50%   -0.5674521
# 97.5% -0.1869978
