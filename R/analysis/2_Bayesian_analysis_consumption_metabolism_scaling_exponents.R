#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.23: Max Lindmark
#
# - This code uses the rstanarm package to fit Bayesian hierarchical models using 
#   lme4 syntax. 
# 
# A. Load libraries
#
# B. Fit model for consumption exponents ~ temperature
#
# C. Fit model for metabolic rate exponents ~ temperature
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is based on the following guides:
# https://mc-stan.org/rstanarm/articles/rstanarm.html
# http://www.tqmp.org/RegularArticles/vol14-2/p099/p099.pdf
# https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html#fn2
# https://cran.r-project.org/web/packages/rstanarm/rstanarm.pdf

# See also this for plotting! 
# http://mc-stan.org/bayesplot/articles/plotting-mcmc-draws.html
# https://mjskay.github.io/tidybayes/articles/tidy-rstanarm.html
# https://www.tjmahr.com/visualizing-uncertainty-rstanarm/

# --- Q is to include or not include n=1 species. In the mixed model, they are not included, but here I suppose they are treated as missing values? in which case they are by default assigned missing values? Not sure it makes sense conceptually though...

# --- Go through explore csv again, which species should I not include initially?

# --- Add pref temp env <- pref temp if no env temp (see growth scripts)

# A. LOAD LIBRARIES & READ DATA ====================================================
rm(list = ls())

# Provide package names
pkgs <- c("mlmRev",
          "tidyr", 
          "lme4",
          "tidylog",
          "rstanarm",
          "RCurl",
          "magrittr",
          "ggplot2",
          "broom",
          "readxl",
          "dplyr",
          "RColorBrewer",
          "bayesplot",
          "modelr",
          "tidybayes",
          "viridis")

# Install packages
# if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(pkgs, rownames(installed.packages())))
# }

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
#script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/packageInfo.R", ssl.verifypeer = FALSE)
#eval(parse(text = script))
#pkg_info(pkgs)

# package loadedversion
# 1     bayesplot         1.6.0
# 2         broom         0.4.4
# 3         dplyr         0.8.1
# 4       ggplot2         3.1.1
# 5          lme4        1.1-19
# 6      magrittr           1.5
# 7        mlmRev         1.0-7
# 8        modelr         0.1.4
# 9  RColorBrewer         1.1-2
# 10        RCurl     1.95-4.12
# 11       readxl         1.3.1
# 12     rstanarm        2.18.2
# 13    tidybayes         1.0.4
# 14      tidylog         0.1.0
# 15        tidyr         0.8.3
# 16      viridis         0.5.1

# Also add patchwork
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Read data
con <- read_excel("data/consumption_scaling_data.xlsx")
met <- read_excel("data/metabolism_scaling_data.xlsx")

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- rbind(con, met)

glimpse(dat)

# Which cols to make numeric?
cols = c(1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

# Normalize variables
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate)

# Relative to mid temp in environment
dat$env_temp_mid_norm <- dat$temp_c - dat$env_temp_mid

# Create abbreviated species name for plotting
# Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Filter intraspecific data consumption
datc <- dat %>% 
  filter(significant_size == "Y" & rate == "consumption") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Filter intraspecific data metabolism
datm <- dat %>% 
  filter(significant_size == "Y" & rate == "metabolism") %>% 
  group_by(common_name) %>% 
  filter(n()>1)


# B. CONSUMPTION ===================================================================
# Plot data again
nb.cols <- length(unique(datc$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(datc, aes(env_temp_mid_norm, b, color = species)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) +
NULL

#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. 

# --- Don't rely on defaults, specify them for clarity and if defaults change...
stanlmerCon <- stan_lmer(formula = b ~ env_temp_mid_norm + (1 | species_ab), 
                         data = datc,
                         iter = 3000,
                         seed = 8194)

# Check model summary
summary(stanlmerCon, digits = 4)
summary(stanlmerCon, digits = 4, probs = c(0.1, 0.9))

# Summary of priors used
prior_summary(object = stanlmerCon)


#** Plot predictions ===============================================================
color_scheme_set("viridisC")

# Posterior-predictive checks
pp_check(stanlmerCon, nreps = 100) +
  xlab("Size-scaling exponents") +
  theme_classic(base_size = 18) +
  theme(legend.text.align = 0)

# Plot species-varying intercepts
stanlmerCon %>%
  spread_draws(`(Intercept)`, b[,group]) %>%
  mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = group, x = condition_mean)) +
  geom_vline(aes(xintercept = stanlmerCon$coefficients[1]), 
             color = "red", alpha = 0.6, size = 1) +
  geom_halfeyeh() +
  theme_classic(base_size = 16)

# Plot posterior densities with credible intervals
posteriorC <- as.array(stanlmerCon)
dfsumC <- data.frame(summary(stanlmerCon, digits = 4, probs = c(0.025, 0.975, 0.1, 0.9)))

mcmc_dens(posteriorC, pars = c("env_temp_mid_norm")) + 
  geom_vline(xintercept = dfsumC$mean[2],
             color = "black", alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = dfsumC$X2.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumC$X97.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumC$X10.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = dfsumC$X90.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = 0,
             color = viridis(5)[5], size = 0.9, linetype = 1) +
  labs(x = "Change In Mass Scaling Exponent per unit C") +
  theme_classic(base_size = 14)


#** Plot overall prediction and data =============================================
# https://www.tjmahr.com/visualizing-uncertainty-rstanarm/
fits <- stanlmerCon %>% 
  as_data_frame %>% 
  #rename(intercept = `(Intercept)`) %>% 
  select(-sigma)

glimpse(fits)

reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)

n_draws <- 500
alpha_level <- .05
col_draw <- "grey30"
col_median <- reds[4]

nb.cols <- length(unique(datc$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(datc) + 
  aes(x = env_temp_mid_norm, y = b, color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = env_temp_mid_norm), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in red
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$env_temp_mid_norm), 
              size = 1.4, color = col_median) +
  geom_point(size = 4.5) +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  annotate(geom = "text", x = 15, y = 1.1, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$env_temp_mid_norm), digits = 5), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Normalzied temperature (midpoint)",
       y = "Mass-scaling exponent (consumption)") +
  NULL

#ggsave("figs/scatter_con_b.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


# C. METABOLISM ====================================================================
# Plot data again
nb.cols <- length(unique(datm$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(datm, aes(env_temp_mid_norm, b, color = species)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) +
NULL


#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. 
datm <- datm %>% drop_na(env_temp_mid_norm)

# --- Don't rely on defaults, specify them for clarity and if defaults change...
stanlmerMet <- stan_lmer(formula = b ~ env_temp_mid_norm + (1 | species_ab), 
                         data = datm,
                         seed = 8194)

# Check model summary
summary(stanlmerMet, digits = 4)
summary(stanlmerMet, digits = 4, probs = c(0.1, 0.9))

# Summary of priors used
prior_summary(object = stanlmerMet)


#** Plot predictions ===============================================================
color_scheme_set("viridis")

# Posterior-predictive checks
# Not very good...
pp_check(stanlmerMet, nreps = 100) +
  xlab("Size-scaling exponents") +
  theme_classic(base_size = 18) +
  theme(legend.text.align = 0)

# Plot species-varying intercepts
stanlmerMet %>%
  spread_draws(`(Intercept)`, b[,group]) %>%
  mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = group, x = condition_mean)) +
  geom_vline(aes(xintercept = stanlmerMet$coefficients[1]), 
             color = "red", alpha = 0.6, size = 1) +
  geom_halfeyeh() +
  theme_classic(base_size = 16)

# Plot posterior densities with credible intervals
posteriorM <- as.array(stanlmerMet)
dfsumM <- data.frame(summary(stanlmerMet, digits = 4, probs = c(0.025, 0.975, 0.1, 0.9)))

mcmc_dens(posteriorM, pars = c("env_temp_mid_norm")) + 
  geom_vline(xintercept = dfsumM$mean[2],
             color = "black", alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = dfsumM$X2.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumM$X97.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumM$X10.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = dfsumM$X90.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = 0,
             color = viridis(5)[5], size = 0.9, linetype = 1) +
  labs(x = "Change In Mass Scaling Exponent per unit C") +
  theme_classic(base_size = 14)


#** Plot overall prediction and data =============================================
# https://www.tjmahr.com/visualizing-uncertainty-rstanarm/
fits <- stanlmerMet %>% 
  as_data_frame %>% 
  #rename(intercept = `(Intercept)`) %>% 
  select(-sigma)

glimpse(fits)

reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)

n_draws <- 500
alpha_level <- .05
col_draw <- "grey30"
col_median <- reds[4]

nb.cols <- length(unique(datm$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(datm) + 
  aes(x = env_temp_mid_norm, y = b, color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = env_temp_mid_norm), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in red
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$env_temp_mid_norm), 
              size = 1.4, color = col_median) +
  geom_point(size = 4.5) +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  annotate(geom = "text", x = -5, y = 0.5, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$env_temp_mid_norm), digits = 5), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Normalzied temperature (midpoint)",
       y = "Mass-scaling exponent (metabolism)") +
  NULL

#ggsave("figs/scatter_met_b.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


# D. PLOT TOGETHER =================================================================
color_scheme_set("viridisC")

c1 <- mcmc_dens(posteriorC, pars = c("env_temp_mid_norm")) + 
  geom_vline(xintercept = dfsumC$mean[2],
             color = "black", alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = dfsumC$X2.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumC$X97.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumC$X10.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = dfsumC$X90.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = 0,
             color = viridis(5)[5], size = 0.9, linetype = 1) +
  labs(x = "") +
  xlim(-0.012, 0.012) + 
  ggtitle("Maximum Consumption Rate") +
  theme_classic(base_size = 13)

color_scheme_set("viridis")

m1 <- mcmc_dens(posteriorM, pars = c("env_temp_mid_norm")) + 
  geom_vline(xintercept = dfsumM$mean[2],
             color = "black", alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = dfsumM$X2.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumM$X97.5.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumM$X10.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = dfsumM$X90.[2],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = 0,
             color = viridis(5)[5], size = 0.9, linetype = 1) +
  labs(x = "Change In Mass Scaling Exponent per unit C") +
  xlim(-0.012, 0.012) +
  ggtitle("Metabolic Rate") +
  theme_classic(base_size = 13)

p <- c1/m1

#ggsave("figs/posterior_c.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")
