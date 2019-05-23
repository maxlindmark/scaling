#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.23: Max Lindmark
#
# - This code uses the rstanarm package to fit Bayesian hierarchical models using 
#   lme4 syntax. 
# 
# A. Load libraries
#
# B. Fit model
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

# *** Q is to include or not include n=1 species. In the mixed model, they are not included, but here I suppose they are treated as missing values? in which case they are by default assigned missing values? Not sure it makes sense conceptually though...

# *** Go through explore csv again, which species should I not include initially?

# *** Add pref temp env <- pref temp if no env temp (see growth scripts)

#======== A. LOAD LIBRARIES & READ DATA ============================================
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
          "plyr",
          "RColorBrewer",
          "bayesplot",
          "modelr",
          "tidybayes")

library(viridis)

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
#script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
#eval(parse(text = script))
#pkg_info(pkgs)

# package loadedversion
# 1     bayesplot         1.6.0
# 2         broom         0.4.4
# 3         dplyr       0.8.0.1
# 4       ggplot2         3.1.1
# 5          lme4        1.1-19
# 6      magrittr           1.5
# 7        mlmRev         1.0-7
# 8        modelr         0.1.4
# 9          plyr         1.8.4
# 10 RColorBrewer         1.1-2
# 11        RCurl     1.95-4.12
# 12       readxl         1.3.1
# 13     rstanarm        2.18.2
# 14    tidybayes         1.0.4
# 15      tidylog         0.1.0

# Read data
con <- read_excel("data/consumption_scaling_data.xlsx")
met <- read_excel("data/metabolism_scaling_data.xlsx")

con$rate <- "consumption"
met$rate <- "metabolism"

dat <- rbind(con, met)

glimpse(dat)

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
s_con <- dat %>% 
  filter(significant_size == "Y" & rate == "consumption") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Filter intraspecific data metabolism
s_met <- dat %>% 
  filter(significant_size == "Y" & rate == "metabolism") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# If env temp == NA, use mid point of pref.


#======== B. CONSUMPTION ===========================================================
# Plot data again
nb.cols <- length(unique(s_con$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_con, aes(env_temp_mid_norm, b, color = species)) + 
  geom_point(size = 5) +
  # geom_line(size = 2) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors)
NULL


#====**** Set up stan model ========================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model, with the following syntax: lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

# *** dont rely on defaults, specify them for clarity and if defaults change...
s_con <- s_con %>% drop_na(env_temp_mid_norm)

c1_stanlmer <- stan_lmer(formula = b ~ env_temp_mid_norm + (env_temp_mid_norm | species_ab), 
                         data = s_con,
                         seed = 8194)

# Summary of priors used
prior_summary(object = c1_stanlmer)

# Print model summary
print(c1_stanlmer, digits = 3)


#====**** Check sampling quality & model convergence ===============================
summary(c1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# Compare the observed distribution of b scores (dark blue line) to 100 simulated datasets from the posterior predictive distribution (light blue lines)
pp_check(c1_stanlmer, nreps = 100)
# pp_check(c1_stanlmer, plotfun = "stat_grouped", stat = "median", group = "species")

# Plot Rhat (1 is good)
plot(c1_stanlmer, "rhat")

# Plot ess
plot(c1_stanlmer, "ess")

# *** See metabolism code for proper diagnostics


#====**** Extract the posterior draws for all parameters ===========================
# *** testing alternative way of plotting posterior including densities!
# posterior <- as.array(c1_stanlmer)
# 
# test <- data.frame(species = sort(unique(s_con$species_ab)))
# test <- test[-which(test$species == "E.coioides"), ]
# 
# 
# par <- paste("b[env_temp_mid_norm species_ab:", 
#              test, 
#              "]", 
#              sep = "")
# 
# mcmc_areas(
#   posterior,
#   pars = c(par),
#   prob = 0.8, # 80% intervals
#   prob_outer = 0.99, # 99%
#   point_est = "mean"
# ) + theme_classic(base_size = 12)


sims <- as.matrix(c1_stanlmer)

para_name <- colnames(sims)

para_name # *** need to know all parameters here (check sigma)

# Create df of parameters I want to plot:
species_b_c <- data.frame(species = sort(unique(s_con$species_ab)))

# Paste parameter name + species to get all species-specific parameters (slope in this case)
species_b_c$param <- paste("b[env_temp_mid_norm species:", 
                           sort(unique(species_b_c$species)), 
                           "]", 
                           sep = "")

# Extract each species slope and 90% CI
summaryc1_90 <- tidy(c1_stanlmer, intervals = TRUE, prob =.9, 
                     parameters = "varying")

summaryc1_90 <- summaryc1_90 %>% 
  filter(term == "env_temp_mid_norm")

# Estimate
species_b_c$pred_slope <- summaryc1_90$estimate

# 90% CI
species_b_c$pred_slope_lwr90 <- summaryc1_90$lower
species_b_c$pred_slope_upr90 <- summaryc1_90$upper

# 50% CI
summaryc1_50 <- tidy(c1_stanlmer, intervals = TRUE, prob =.5, 
                     parameters = "varying")

summaryc1_50 <- summaryc1_50 %>% 
  filter(term == "env_temp_mid_norm")

species_b_c$pred_slope_lwr50 <- summaryc1_50$lower
species_b_c$pred_slope_upr50 <- summaryc1_50$upper


#====**** Plot species-slopes ======================================================
blues <- colorRampPalette(brewer.pal(5, "Blues"))(5)
# pal2 <- viridis(n = 4)
pal <- colorRampPalette(brewer.pal(8, "Paired"))(12)
pal2 <- c("#fdb863", "#b2abd2", "#5e3c99", "#e66101") # From colorbrewer website

# Add overall mean  and 90% & 50% CI
g_mean_c <- summary(c1_stanlmer, probs = c(0.05, 0.95), digits = 5)[2] # Get the slope
g_mean_ci90 <- posterior_interval(c1_stanlmer, prob = 0.90, pars = "env_temp_mid_norm")
g_mean_ci50 <- posterior_interval(c1_stanlmer, prob = 0.50, pars = "env_temp_mid_norm")

# Scale parameters to make them slope, not diff from mean!
pdat_c <- species_b_c
pdat_c[, 3:7] <- sweep(species_b_c[, 3:7], 2, g_mean_c, FUN = "+")

# Group species with "stronger" evidence
pdat_c$sig <- ifelse(pdat_c$pred_slope < 0 & pdat_c$pred_slope_upr50 < 0, 
                     2, -9)

pdat_c$sig <- ifelse(pdat_c$pred_slope > 0 & pdat_c$pred_slope_lwr50 > 0, 
                     1, pdat_c$sig)

pdat_c$sig <- ifelse(pdat_c$pred_slope < 0 & pdat_c$pred_slope_upr90 < 0, 
                     3, pdat_c$sig)

pdat_c$sig <- ifelse(pdat_c$pred_slope > 0 & pdat_c$pred_slope_lwr90 > 0, 
                     4, pdat_c$sig)

pdat_c$sig_f <- as.factor(pdat_c$sig)

# Group positive and negative slopes
pdat_c$sign <- ifelse(pdat_c$pred_slope < 0, "neg", "pos")

# Reorder based on predicted slope and plot!
ggplot(pdat_c, aes(reorder(species, pred_slope), pred_slope, 
                   shape = sign, color = sig_f)) +
  geom_rect(data = pdat_c, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci90[1], ymax = g_mean_ci90[2],
            color = "grey93", fill = "grey93") + # overplotting many rectangles here..
  geom_rect(data = pdat_c, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci50[1], ymax = g_mean_ci50[2],
            color = "grey83", fill = "grey83") + # overplotting many rectangles here..
  geom_hline(yintercept = 0, col = "grey1", linetype = 1, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = g_mean_c, col = pal[8], linetype = "twodash", size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr90, ymax = pred_slope_upr90), 
                color = "gray70", size = 0.5, width = 0) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr50, ymax = pred_slope_upr50), 
                color = "gray50", size = 0.8, width = 0) +
  scale_color_manual(values = pal[c(4,1,2)]) + # Need to watch out for this, depends on # levels
  #scale_color_manual(values = c("#b2abd2", "#fdb863", "#e66101", "#5e3c99")) + 
  #scale_color_manual(values = pal2) + # Need to watch out for this, depends on levels
  geom_point(size = 4, col = "gray50") +
  geom_point(data = filter(pdat_c, sig > 0), 
             aes(reorder(species, pred_slope), pred_slope, shape = sign, color = sig_f), size = 4) +
  xlab("Species") + 
  ylab("Change in mass-scaling exponent per T (C)") +
  coord_flip() +
  guides(shape = FALSE, color = FALSE) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  theme(aspect.ratio = 1) +
  #ylim(min(pdat$pred_slope_lwr95), -1*min(pdat$pred_slope_lwr95)) +
  #theme(axis.text.y = element_blank()) + # If I end up with too many species
  NULL

ggsave("figs/intraspec_con.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#====**** Plot intercept (exponent at mid temp) ====================================
# *** Do this differently... perhaps the same way as above? But still, I want the overall mean slope as well...
inter <- tidy(c1_stanlmer, intervals = TRUE, prob =.95, 
              parameters = "varying")

inter <- inter %>% 
  filter(term == "(Intercept)")

inter$a <- inter$estimate + summary(c1_stanlmer, probs = c(0.025, 0.975), digits = 5)[1]

ggplot(inter, aes(a)) +
  geom_density(fill = pal[2], color = NA, alpha = 0.8) + 
  coord_cartesian(xlim = c(0.39, 0.91), expand = c(0,0)) +
  theme_classic(base_size = 15) +
  geom_vline(xintercept = summary(c1_stanlmer, probs = c(0.025, 0.975), digits = 5)[1], 
             size = 1.5, color = "gray20") +
  xlab("Intercept") +
  theme(aspect.ratio = 1) +
  NULL
  
# *** Need to check the arguments of this function...
pp_check(c1_stanlmer) + 
  theme(aspect.ratio = 1) +
  theme_classic(base_size = 15)


#====**** Plot overall prediction and data =========================================
# https://www.tjmahr.com/visualizing-uncertainty-rstanarm/
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

nb.cols <- length(unique(s_con$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_con) + 
  aes(x = env_temp_mid_norm, y = b, color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = env_temp_mid_norm), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in blue
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$env_temp_mid_norm), 
              size = 1.4, color = col_median) +
  geom_point(size = 4.5) +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Normalzied temperature midpoint",
       y = "Mass-scaling exponent (consumption)") +
  NULL

ggsave("figs/scatter_con_b.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#======== B. METABOLISM ============================================================
# Plot data again
nb.cols <- length(unique(s_met$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_met, aes(env_temp_mid_norm, b, color = species)) + 
  geom_point(size = 5) +
  geom_line(size = 2) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) +
NULL

# Remove species with na temperature!
s_met2 <- s_met
s_met2$env_temp_mid_norm
s_met2 <- s_met2 %>% drop_na(env_temp_mid_norm)
s_met2$env_temp_mid_norm
s_met2$b


#====**** Set up stan model ========================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model, with the following syntax: lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

# *** dont rely on defaults, specify them for clarity and if defaults change...
m1_stanlmer <- stan_lmer(formula = b ~ env_temp_mid_norm + (env_temp_mid_norm | species_ab), 
                         data = s_met2,
                         seed = 8194,
                         adapt_delta = 0.999999999999) # This is essentially reducing the step size. Could be the reason for the autocorrelated samples... We did not need to do this for consumption. This is likely not optimal. Should wait for more data though before digging into changing error structure..

# Summary of priors used
prior_summary(object = m1_stanlmer)

# Print model summary
print(m1_stanlmer, digits = 3)


#====**** Check sampling quality & model convergence ===============================
summary(m1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# Compare the observed distribution of b scores (dark blue line) to 100 simulated datasets from the posterior predictive distribution (light blue lines)
pp_check(m1_stanlmer, nreps = 100)
# pp_check(c1_stanlmer, plotfun = "stat_grouped", stat = "median", group = "species")

# Plot Rhat (1 is good)
plot(m1_stanlmer, "rhat")

# Plot ess
# *** This suggest low number of effective simulations, high autocorrelation...
plot(m1_stanlmer, "ess")

posterior <- as.array(m1_stanlmer)

# *** Add this kind of diagnostics to Cmax as well...
mcmc_trace(posterior, pars = c("env_temp_mid_norm"),
           facet_args = list(ncol = 1, strip.position = "left"))


#====**** Extract the posterior draws for all parameters ===========================
sims <- as.matrix(m1_stanlmer)

para_name <- colnames(sims)

para_name # *** need to know all parameters here (check sigma)

# Create df of parameters I want to plot:
species_b_m <- data.frame(species = sort(unique(s_met2$species_ab))) # s_met2 because of filtering!

# Paste parameter name + species to get all species-specific parameters (slope in this case)
species_b_m$param <- paste("b[env_temp_mid_norm species:", 
                           sort(unique(species_b_m$species_ab)), 
                           "]", 
                           sep = "")

# Extract each species slope and 90% CI
summarym1_90 <- tidy(m1_stanlmer, intervals = TRUE, prob =.9, 
                     parameters = "varying")

summarym1_90 <- summarym1_90 %>% 
  filter(term == "env_temp_mid_norm")

# Estimate
species_b_m$pred_slope <- summarym1_90$estimate

# 90% CI
species_b_m$pred_slope_lwr90 <- summarym1_90$lower
species_b_m$pred_slope_upr90 <- summarym1_90$upper

# 50% CI
summarym1_50 <- tidy(m1_stanlmer, intervals = TRUE, prob =.5, 
                     parameters = "varying")

summarym1_50 <- summarym1_50 %>% 
  filter(term == "env_temp_mid_norm")

species_b_m$pred_slope_lwr50 <- summarym1_50$lower
species_b_m$pred_slope_upr50 <- summarym1_50$upper


#====**** Plot species-slopes ======================================================
#blues <- colorRampPalette(brewer.pal(8, "Blues"))(9)
# pal2 <- viridis(n = 4)
pal <- colorRampPalette(brewer.pal(8, "Paired"))(12)

# Add overall mean  and 90% & 50% CI
g_mean_m <- summary(m1_stanlmer, probs = c(0.05, 0.95), digits = 5)[2] # Get the slope
g_mean_ci90 <- posterior_interval(m1_stanlmer, prob = 0.90, pars = "env_temp_mid_norm")
g_mean_ci50 <- posterior_interval(m1_stanlmer, prob = 0.50, pars = "env_temp_mid_norm")

# Scale parameters to make them slope, not diff from mean!
pdat_m <- species_b_m
pdat_m[, 3:7] <- sweep(species_b_m[, 3:7], 2, g_mean_m, FUN = "+")

# Group species with "stronger" evidence
pdat_m$sig <- ifelse(pdat_m$pred_slope < 0 & pdat_m$pred_slope_upr50 < 0, 
                     2, -9)

pdat_m$sig <- ifelse(pdat_m$pred_slope > 0 & pdat_m$pred_slope_lwr50 > 0, 
                     1, pdat_m$sig)

pdat_m$sig <- ifelse(pdat_m$pred_slope < 0 & pdat_m$pred_slope_upr90 < 0, 
                     3, pdat_m$sig)

pdat_m$sig <- ifelse(pdat_m$pred_slope > 0 & pdat_m$pred_slope_lwr90 > 0, 
                     4, pdat_m$sig)

pdat_m$sig_f <- as.factor(pdat_m$sig)

# Group positive and negative slopes
pdat_m$sign <- ifelse(pdat_m$pred_slope < 0, "neg", "pos")

# Reorder based on predicted slope and plot!
ggplot(pdat_m, aes(reorder(species, pred_slope), pred_slope, 
                   shape = sign, color = sig_f)) +
  geom_rect(data = pdat_m, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci90[1], ymax = g_mean_ci90[2],
            color = "grey93", fill = "grey93") + # overplotting many rectangles here..
  geom_rect(data = pdat_m, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci50[1], ymax = g_mean_ci50[2],
            color = "grey83", fill = "grey83") + # overplotting many rectangles here..
  geom_hline(yintercept = 0, col = "grey1", linetype = 1, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = g_mean_m, col = pal[8], linetype = "twodash", size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr90, ymax = pred_slope_upr90), 
                color = "gray70", size = 0.5, width = 0) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr50, ymax = pred_slope_upr50), 
                color = "gray50", size = 0.8, width = 0) +
  scale_color_manual(values = pal[c(1,1,4)]) + # Need to watch out for this, depends on # levels
  #scale_color_manual(values = pal2[c(2,3,4)]) + # Need to watch out for this, depends on # levels "pink"
  geom_point(size = 4, col = "gray50") +
  geom_point(data = filter(pdat_m, sig > 0), 
             aes(reorder(species, pred_slope), pred_slope, shape = sign, color = sig_f), size = 4) +
  xlab("Species") + 
  ylab("Change in mass-scaling exponent per T (C)") +
  coord_flip() +
  guides(shape = FALSE, color = FALSE) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  theme(aspect.ratio = 1) +
  ylim(-0.025, 0.02) + # This was not need for Cmax!
  #ylim(min(pdat$pred_slope_lwr95), -1*min(pdat$pred_slope_lwr95)) +
  #theme(axis.text.y = element_blank()) + # If I end up with too many species
  NULL

ggsave("figs/intraspec_met.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#====**** Plot intercept (exponent at mid temp) ====================================
inter <- tidy(m1_stanlmer, intervals = TRUE, prob =.95, 
              parameters = "varying")

inter <- inter %>% 
  filter(term == "(Intercept)")

inter$a <- inter$estimate + summary(m1_stanlmer, probs = c(0.025, 0.975), digits = 5)[1]

ggplot(inter, aes(a)) +
  geom_density(fill = pal[2], color = NA, alpha = 0.8) + 
  coord_cartesian(xlim = c(0.815, 0.865), expand = c(0,0)) +
  theme_classic(base_size = 15) +
  geom_vline(xintercept = summary(m1_stanlmer, probs = c(0.025, 0.975), digits = 5)[1], 
             size = 1.5, color = "gray20") +
  xlab("Intercept") +
  theme(aspect.ratio = 1) +
  NULL

# *** Need to check the arguments of this function...
pp_check(m1_stanlmer) + 
  theme(aspect.ratio = 1) +
  theme_classic(base_size = 15)


#====**** Plot overall prediction and data =========================================
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

nb.cols <- length(unique(s_met$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_met) + 
  aes(x = env_temp_mid_norm, y = b, color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = env_temp_mid_norm), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in blue
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$env_temp_mid_norm), 
              size = 1.4, color = col_median) +
  geom_point(size = 4.5) +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  ylim(0.4, 1.25) + # removing outlier
  labs(x = "Normalzied temperature midpoint",
       y = "Mass-scaling exponent (metabolism)") +
  NULL

ggsave("figs/scatter_met_b.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")






#------------------------------------ Testing raw stan 
# The idea is to exactly reproduce the rstanarm model!

# Check this for how to extract posterior distributions: http://www.maths.bath.ac.uk/~jjf23/stan/rbd.html

# This is the main guide (paper): http://jakewestfall.org/misc/SorensenEtAl.pdf
# and blog: https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html

# More examples here: http://www.maths.bath.ac.uk/~jjf23/stan/

# Mathy paper here: https://arxiv.org/pdf/1506.06201.pdf and here https://pdfs.semanticscholar.org/d828/2c0cd55431b4fc7eaa0bceaa880ba1bafc34.pdf

# Write a Stan Program
stanmodel1 = "
data {
int<lower=0> J;          // number of schools 
real y[J];               // estimated treatment effects
real<lower=0> sigma[J];  // s.e. of effect estimates 
}
parameters {
real mu; 
real<lower=0> tau;
vector[J] eta;
}
transformed parameters {
vector[J] theta;
theta = mu + tau * eta;
}
model {
target += normal_lpdf(eta | 0, 1);
target += normal_lpdf(y | theta, sigma);
}
"

# Format data for Stan:
stanDat <- list(subj = as.integer(s_con$subj),
                item = as.integer(s_con$item),
                rt = s_con$rt,
                so = s_con$so,
                N = nrow(s_con),
                J = nlevels(s_con$subj),
                K = nlevels(s_con$item))


# Sample from the Posterior Distribution
library(rstan)

fit1 <- stan(
  file = "schools.stan",  # Stan program
  data = schools_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
)





#------------------------------------ Continue here -----------------------------

#------------------------------------# Freq. Mixed model
me1 <- lmer(b ~ env_temp_mid_norm + (env_temp_mid_norm | species), s_con)
summary(me1)
coef(me1)

me2 <- lmer(b ~ env_temp_mid_norm + (env_temp_mid_norm | species), s_met2)
summary(me1)
coef(me1)

## Singular fit!

#------------------------------------



#------------------------------------
# Multiple LM's
# Test same data with multiple lm's instead:
p <- c()
r <- c()
s <- c()
t <- data.frame()

for (i in unique(species_b$species)){
  p  <- data.frame(subset(s_con, species == i))  
  
  m1 <- summary(lm(b ~ env_temp_mid_norm, data = p))
  s  <- data.frame(slope = m1$coefficients[2],
                   species = i,
                   p = m1$coefficients[,4][2],
                   se = m1$coefficients[,2][2])
  t  <- rbind(t, s)
}

t

# Now plot coefficients from species-specific lm's
t$upper <- t$slope + t$se * 1.96 
t$lower <- t$slope - t$se * 1.96

ggplot(t, aes(reorder(species, slope), slope)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, col = "red", linetype = 2, size = 1.3) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), 
                width = 0, color = "gray55", size = 1.5) + 
  guides(color = FALSE) + 
  xlab("Species") + 
  ylab("Slope") +
  coord_flip() +
  theme_classic(base_size = 15) +
  NULL
#------------------------------------





