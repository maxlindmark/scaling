#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.29: Max Lindmark
#
# - This code uses the rstanarm package to fit Bayesian hierarchical models using 
#   lme4 syntax. 
# 
# A. Load libraries & read data
#
# B. Fit model for consumption ~ temperature
#
# C. Fit model for metabolic rate ~ temperature
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
          "tidybayes")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

library(viridis)

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

# Relative to "preferred" temp (not all species have this info)
dat$pref_temp_mid_norm <- dat$temp_c - dat$pref_temp_mid

# Create inverse temp
dat$inv_temp <- 1/((dat$temp_c + 273.15) * 8.617332e-05)

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass_g / dat$w_max_published_g
# ---- Stechlin cisco har larger size than max*

# Create data with with species that have also temperature repliates within species (intra-specific analysis)
s_datc <- data.frame(
  dat %>% 
    filter(rate == "consumption") %>% 
    group_by(species) %>% 
    mutate(unique_t = n_distinct(env_temp_mid_norm)) %>% 
    filter(unique_t > 1)
)

s_datm <-  data.frame(
  dat %>% 
    filter(rate == "metabolism") %>% 
    group_by(species) %>% 
    mutate(unique_t = n_distinct(env_temp_mid_norm)) %>% 
    filter(unique_t > 1)
)

# Plot data again (log rate vs inverse Boltzmann temp:
pal <- viridis(n = 5)
reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)


p1 <- s_datc %>% 
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
p2 <- s_datm %>% 
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
s_datc <- s_datc %>% filter(env_temp_mid_norm < 13)
s_datm <- s_datm %>% filter(unit == "mg O2/h" & env_temp_mid_norm < 12)

summary(lm(log(y) ~ inv_temp, data = s_datm))
summary(lm(log(y) ~ inv_temp, data = s_datc))


# B. Consumption ===================================================================
#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model: lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

# Test classic mixed model
me1 <- lmer(log(y) ~ inv_temp + (inv_temp | species), s_datc)
summary(me1)
coef(me1)

# --- dont rely on defaults, specify them for clarity and if defaults change...
c1_stanlmer <- stan_lmer(formula = log(y) ~ inv_temp + (inv_temp | species), 
                         data = s_datc,
                         seed = 8194)

# Summary of priors used
prior_summary(object = c1_stanlmer)

# Print model summary
print(c1_stanlmer, digits = 3)


#** Check sampling quality & model convergence =====================================
summary(c1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)


#** Plot overall prediction and data ===============================================
fits <- c1_stanlmer %>% 
  as_data_frame %>% 
  #rename(intercept = `(Intercept)`) %>% 
  select(-sigma)

glimpse(fits)

reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)

n_draws <- 500
alpha_level <- .1
col_draw <- "grey40"
col_median <- reds[4]

nb.cols <- length(unique(s_datc$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_datc) + 
  aes(x = inv_temp, y = log(y), color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = inv_temp), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in blue
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$inv_temp), 
              size = 1.4, color = col_median) +
  geom_jitter(size = 4.5, width = 0, alpha = 0.5, height = 0) +
  annotate(geom = "text", x = 38.4, y = -3.4, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$inv_temp), digits = 2), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Inverse temperature [1/kT]", 
       y = "log(maximum consumption rate [g/d])") +
  NULL

#ggsave("figs/activation_energy_con.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


# B. Metabolism ====================================================================
#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model: lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

# Test classic mixed model
me1 <- lmer(log(y) ~ inv_temp + (inv_temp | species), s_datm)
summary(me1)
coef(me1)

# --- dont rely on defaults, specify them for clarity and if defaults change...
m1_stanlmer <- stan_lmer(formula = log(y) ~ inv_temp + (inv_temp | species), 
                         data = s_datm,
                         seed = 8194,
                         adapt_delta = 0.99)

# Summary of priors used
prior_summary(object = m1_stanlmer)

# Print model summary
print(m1_stanlmer, digits = 3)


#** Check sampling quality & model convergence =====================================
summary(m1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# Compare the observed distribution of b scores (dark blue line) to 100 simulated datasets from the posterior predictive distribution (light blue lines)
pp_check(m1_stanlmer, nreps = 100)
# pp_check(c1_stanlmer, plotfun = "stat_grouped", stat = "median", group = "species")

# Plot Rhat (1 is good)
plot(m1_stanlmer, "rhat")

# Plot ess
plot(m1_stanlmer, "ess")


#** Plot overall prediction and data ===============================================
fits <- m1_stanlmer %>% 
  as_data_frame %>% 
  #rename(intercept = `(Intercept)`) %>% 
  select(-sigma)

glimpse(fits)

reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)

n_draws <- 500
alpha_level <- .1
col_draw <- "grey40"
col_median <- reds[4]

nb.cols <- length(unique(s_datm$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_datm) + 
  aes(x = inv_temp, y = log(y), color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = inv_temp), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in blue
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$inv_temp), 
              size = 1.4, color = col_median) +
  geom_jitter(size = 4.5, width = 0.1, alpha = 0.5, height = 0) +
  annotate(geom = "text", x = 38.8, y = -1.8, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$inv_temp), digits = 2), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Inverse temperature [1/kT]", 
       y = expression(paste("log(metabolic rate [mg  ", O[2], "/h])"))) +
  NULL

#ggsave("figs/activation_energy_meta.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")



