#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.29: Max Lindmark
#
# - This code uses the rstanarm package to fit Bayesian hierarchical models using 
#   lme4 syntax. 
# 
# A. Load libraries & read data
#
# B. Explore and normalize data
#
# C. Fit model for consumption ~ temperature
#
# D. Fit model for metabolic rate ~ temperature
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

# Set seed
set.seed(8194)

# Provide package names
pkgs <- c("mlmRev",
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
# 13    tidybayes         1.1.0
# 14      tidylog         0.1.0
# 15      viridis         0.5.1

# Also add patchwork
# devtools::install_github("thomasp85/patchwork")
library(patchwork)

# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

con <- read_excel("data/consumption_data.xlsx")
met <- read_excel("data/metabolism_data.xlsx")

# Which cols to make numeric?
cols = c(1,2,3,12,13,14,15,16,17)

con[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
met[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

colnames(con)
head(con)

con <- con %>% rename(y = consumption)
met <- met %>% rename(y = metabolic_rate)

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- bind_rows(con, met)

head(dat)


# B. EXPLORE DATA ==================================================================
#** Normalize variables ============================================================
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate)

# Relative to temp in environment
dat$env_temp_mid_norm <- dat$temp_c - dat$env_temp_mid

# Create inverse temp
dat$inv_temp <- 1/((dat$temp_c + 273.15) * 8.617332e-05)

# Mean-center inverse temp
dat$inv_temp_ct <- dat$inv_temp - mean(dat$inv_temp)

# Create inverse temp with based on centered temperature within species
dat$inv_temp_sp <- 1/((dat$env_temp_mid_norm + 273.15) * 8.617332e-05)

# Mean-center inverse temp with based on centered temperature within species
dat$inv_temp_sp_ct <- dat$inv_temp_sp - mean(dat$inv_temp_sp, na.rm = T)

# Center log(mass) relative to 100 (size-coefficient - size scaling exponent - is then for 100g individuals). 
# Could also mean across species (130g)...
dat$ln_mass_g <- log(dat$mass_g)
dat$ln_mass_g_ct <- dat$ln_mass_g - mean(dat$ln_mass_g) # log(100)# log(mean(dat$mass_g))

# Calculate normalized mass & mean-centre
dat$mass_norm <- dat$mass_g / dat$w_max_published_g
dat$ln_mass_norm_ct <- log(dat$mass_norm) - log(mean(dat$mass_norm))

# Create data with with species that have also temperature repliates within species (intra-specific analysis)
datc <- data.frame(
  dat %>% 
    dplyr::filter(rate == "consumption") %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(unique_t = n_distinct(env_temp_mid_norm)) %>% 
    dplyr::filter(unique_t > 1)
)

datm <-  data.frame(
  dat %>% 
    dplyr::filter(rate == "metabolism") %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(unique_t = n_distinct(env_temp_mid_norm)) %>% 
    dplyr::filter(unique_t > 1)
)


# Plot data again (log rate vs inverse Boltzmann temp:
pal <- viridis(n = 5)
reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)

p1 <- datc %>% 
  filter(env_temp_mid_norm < 13) %>% 
  ggplot(., aes(inv_temp, log(y), color = species)) +
  scale_color_viridis(discrete = TRUE) + 
  guides(color = FALSE) +
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.4) +
  labs(x = "Inverse temperature [1/kT]", y = "log(consumption [g/day])") +
  stat_smooth(method = "lm", size = 3, alpha = 0.8, geom = "line", color = reds[4]) +
  NULL

# Metabolism
# --- Filter also the same measurement? Routine vs Standard? Or, because I have species as a random factor, I account for that?
# --- Or add "type of metabolism" as random intercept?

p2 <- datm %>% 
  filter(unit == "mg O2/h" & env_temp_mid_norm < 13) %>% 
  ggplot(., aes(inv_temp, log(y), color = species)) + 
  scale_color_viridis(discrete = TRUE) + 
  guides(color = FALSE) +
  theme_classic(base_size = 15) +
  geom_point(size = 4, alpha = 0.4) +
  labs(x = "Inverse temperature [1/kT]", y = "log(metabolic rate [mg O2/h])") +
  stat_smooth(method = "lm", size = 3, alpha = 0.8, geom = "line", color = reds[4]) +
  NULL

p1 + p2

# Subset data according to filter above (below optimum [see data exploration script], specific unit)
unique(datc$unit)
datc <- datc %>% filter(env_temp_mid_norm < 13)
datm <- datm %>% filter(unit == "mg O2/h" & env_temp_mid_norm < 12)


# C. Consumption ===================================================================
#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model

# Non-normalized mass 
stanlmerCon <- stan_lmer(formula = log(y) ~ inv_temp_sp_ct + ln_mass_g_ct + (1 |species), 
                         data = datc,
                         iter = 3000,
                         seed = 8194) 

# Check model summary
summary(stanlmerCon, digits = 4)

# Summary of priors used
prior_summary(object = stanlmerCon)


#** Plot predictions ================================================================
color_scheme_set("viridisC")

# Plot species-varying intercepts as densities and intervals
stanlmerCon %>%
  spread_draws(`(Intercept)`, b[,group]) %>%
  mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = group, x = condition_mean)) +
  geom_vline(aes(xintercept = stanlmerCon$coefficients[1]), 
             color = "red", alpha = 0.6, size = 1) +
  geom_halfeyeh() +
  ggtitle("Maximum Consumption Rate") +
  theme_classic(base_size = 16)

# Plot posterior predictive check
pp_check(stanlmerCon, nreps = 100) +
  xlab("Maximum consumption rate") +
  theme_classic(base_size = 18) +
  theme(legend.text.align = 0)

# Plot posterior densities with credible intervals
posteriorC <- as.array(stanlmerCon)
dfsumC <- data.frame(summary(stanlmerCon, digits = 4, probs = c(0.025, 0.975, 0.1, 0.9)))

head(dfsumC)

# Activation energy
c1 <- mcmc_dens(posteriorC, pars = c("inv_temp_sp_ct")) + 
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
  # geom_segment(x = dfsumC$X2.5.[2], xend = dfsumC$X2.5.[2],
  #              y = 0, yend = 2, 
  #              color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  # geom_segment(x = dfsumC$X97.5.[2], xend = dfsumC$X97.5.[2],
  #              y = 0, yend = 2, 
  #              color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  # geom_segment(x = dfsumC$X10.[2], xend = dfsumC$X10.[2],
  #              y = 0, yend = 6, 
  #              color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  # geom_segment(x = dfsumC$X90.[2], xend = dfsumC$X90.[2],
  #              y = 0, yend = 6, 
  #              color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 20),
        plot.title = element_text(size = 20)) +
  xlim(c(-0.63, -0.43)) +
  ggtitle("Maximum Consumption Rate") +
  labs(x = "")  

# Size-exponent
c2 <- mcmc_dens(posteriorC, pars = c("ln_mass_g_ct")) + 
  geom_vline(xintercept = dfsumC$mean[3],
             color = "black", alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = dfsumC$X2.5.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumC$X97.5.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumC$X10.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = dfsumC$X90.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 20)) +
  xlim(c(0.5, 0.8)) +
  labs(x = "")  


# D. Metabolism ====================================================================
#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model

# Non-normalized mass 
stanlmerMet <- stan_lmer(formula = log(y) ~ inv_temp_sp_ct + ln_mass_g_ct + (1 |species), 
                         data = datm,
                         iter = 3000,
                         seed = 8194) 

# Check model summary
summary(stanlmerMet, digits = 4)

# Summary of priors used
prior_summary(object = stanlmerMet)


#** Plot predictions ================================================================
color_scheme_set("viridis")

# Plot species-varying intercepts as densities and intervals
stanlmerMet %>%
  spread_draws(`(Intercept)`, b[,group]) %>%
  mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = group, x = condition_mean)) +
  geom_vline(aes(xintercept = stanlmerMet$coefficients[1]), 
             color = "red", alpha = 0.6, size = 1) +
  geom_halfeyeh() +
  ggtitle("Metabolic Rate") +
  theme_classic(base_size = 16)

# Plot posterior predictive check
pp_check(stanlmerMet, nreps = 100) +
  xlab("Metabolic rate") +
  theme_classic(base_size = 18) +
  theme(legend.text.align = 0)

# Plot posterior densities with credible intervals
posteriorM <- as.array(stanlmerMet)
dfsumM <- data.frame(summary(stanlmerMet, digits = 4, probs = c(0.025, 0.975, 0.1, 0.9)))

# Activation energy
m1 <- mcmc_dens(posteriorM, pars = c("inv_temp_sp_ct")) + 
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
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 20),
        plot.title = element_text(size = 20)) +
  xlim(c(-0.63, -0.43)) +
  ggtitle("Metabolic Rate") +
  labs(x = "Activation Energy")  

# Size-exponent
m2 <- mcmc_dens(posteriorM, pars = c("ln_mass_g_ct")) + 
  geom_vline(xintercept = dfsumM$mean[3],
             color = "black", alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = dfsumM$X2.5.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumM$X97.5.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 2) +
  geom_vline(xintercept = dfsumM$X10.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  geom_vline(xintercept = dfsumM$X90.[3],
             color = "black", alpha = 0.6, size = 0.9, linetype = 3) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 20)) +
  xlim(c(0.5, 0.8)) +
  labs(x = "Mass-Scaling Exponent")  


#** All together ===================================================================
p <- (c1+c2)/(m1+m2)

ggsave("figs/posterior_E_b.pdf", plot = p, scale = 1, width = 20, height = 20, units = "cm")

