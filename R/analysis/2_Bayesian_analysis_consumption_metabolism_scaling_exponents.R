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
s_con <- dat %>% 
  filter(significant_size == "Y" & rate == "consumption") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Filter intraspecific data metabolism
s_met <- dat %>% 
  filter(significant_size == "Y" & rate == "metabolism") %>% 
  group_by(common_name) %>% 
  filter(n()>1)


# B. CONSUMPTION ===================================================================
# Plot data again
nb.cols <- length(unique(s_con$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_con, aes(env_temp_mid_norm, b, color = species)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors) +
NULL


#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model: lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

s_con <- s_con %>% drop_na(env_temp_mid_norm)

# --- dont rely on defaults, specify them for clarity and if defaults change...
c1_stanlmer <- stan_lmer(formula = b ~ env_temp_mid_norm + (env_temp_mid_norm | species_ab), 
                         data = s_con,
                         seed = 8194)

# Summary of priors used
prior_summary(object = c1_stanlmer)

# Print model summary
print(c1_stanlmer, digits = 3)


#** Check sampling quality & model convergence =====================================
summary(c1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# Work on this...


#** Extract the posterior draws for all parameters =================================
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


#** Plot species-slopes ============================================================
blues <- colorRampPalette(brewer.pal(5, "Blues"))(5)
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

#ggsave("figs/intraspec_con.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#** Plot overall prediction and data =============================================
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
  annotate(geom = "text", x = 15, y = 1.1, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$env_temp_mid_norm), digits = 4), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Normalzied temperature (midpoint)",
       y = "Mass-scaling exponent (consumption)") +
  NULL

#ggsave("figs/scatter_con_b.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


# C. METABOLISM ====================================================================
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


#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model, with the following syntax: lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

s_met <- s_met %>% drop_na(env_temp_mid_norm)

# --- dont rely on defaults, specify them for clarity and if defaults change...
m1_stanlmer <- stan_lmer(formula = b ~ env_temp_mid_norm + (env_temp_mid_norm | species_ab), 
                         data = s_met,
                         seed = 8194,
                         adapt_delta = 0.999999999999) # This is essentially reducing the step size. Could be the reason for the autocorrelated samples... We did not need to do this for consumption. This is likely not optimal. Should wait for more data though before digging into changing error structure..

# Summary of priors used
prior_summary(object = m1_stanlmer)

# Print model summary
print(m1_stanlmer, digits = 3)


#** Check sampling quality & model convergence =====================================
summary(m1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# Work on this...


#** Extract the posterior draws for all parameters =================================
# Create df of parameters I want to plot:
species_b_m <- data.frame(species = sort(unique(s_met$species_ab)))

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


#** Plot species-slopes ============================================================
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

#ggsave("figs/intraspec_met.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


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
  annotate(geom = "text", x = 10, y = 0.5, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$env_temp_mid_norm), digits = 4), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  ylim(0.4, 1.25) + # removing outlier
  labs(x = "Normalzied temperature (midpoint)",
       y = "Mass-scaling exponent (metabolism)") +
  NULL

#ggsave("figs/scatter_met_b.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#** Plot intercepts (exponent at mid temp) ==========================================
color_scheme_set("teal")
con_post <- as.array(c1_stanlmer)
dimnames(con_post)

pc <- mcmc_areas(con_post,
           pars = c("(Intercept)"),
           prob = 0.8, # 80% intervals
           prob_outer = 0.99, # 99%
           point_est = "mean") +
  theme_classic(base_size = 13) +
  coord_cartesian(xlim = c(0.6, 0.95)) +
  theme(axis.text.y = element_blank(),
        aspect.ratio = 4/5,
        plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle("Consumption") +
  theme() +
  NULL

color_scheme_set("red")
met_post <- as.array(m1_stanlmer)
dimnames(met_post)
pm <- mcmc_areas(met_post,
                 pars = c("(Intercept)"),
                 prob = 0.8, # 80% intervals
                 prob_outer = 0.99, # 99%
                 point_est = "mean") +
  theme_classic(base_size = 13) +
  coord_cartesian(xlim = c(0.6, 0.95)) +
  labs(x = "Scaling exponent", y = "") +
  theme(axis.text.y = element_blank(),
        aspect.ratio = 4/5,
        plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle("Metabolism") +
  # scale_color_manual(values = mycolors[4]) + # remove index if I skip viridis
  # scale_fill_manual(values = mycolors[4]) + 
  NULL

pc/pm

ggsave("figs/posterior_intercept_exponent.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")

