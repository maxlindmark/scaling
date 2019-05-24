#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.23: Max Lindmark
#
# - Analyze growth rate data over temperature and body size
# 
# A. Load libraries & read data
#
# B. Analyze growth data
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
          "plyr",
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

# Print package version
# script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# pkg_info(pkgs)

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

# Will crate a csv that one can read directly once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/growth_data.xlsx"))

dat <- read_excel("data/growth_data.xlsx")

glimpse(dat)

# Which cols to make numeric?
cols = c(1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))


#** Prepare data ===================================================================
# Normalize size by creating new column with size relative to max. In most cases this will be geometric mean, but can also be size class. I will then check if there are any NA's (Walleye) that doesn't have either size and make a data.

# Create abbreviated species name for plotting
# Get first part of name
sp1 <- substring(dat$species, 1, 1)

# Get species name
sp2 <- gsub( ".*\\s", "", dat$species )

dat$species_ab <- paste(sp1, sp2, sep = ".")

# Mass for analysis
dat$mass <- dat$geom_mean_mass_g

# Replace NA with -9...
dat$mass[is.na(dat$mass)] <- -9

dat$mass <- ifelse(dat$mass == -9,
                   dat$mid_mass_sizecl_g,
                   dat$mass)

dat$mass

dat$mass[is.na(dat$mass)] <- -9

dat <- dat %>% filter(mass > 0)

# Now normalize mass with respect to max mass
dat$mass_norm <- dat$mass / dat$w_max_published_g

# Intraspecific data (n>1 species)
s_dat <- dat %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Log mass
s_dat$log10_mass_n <- log10(s_dat$mass_norm)

# Calculate mean optimum temperature within species
s_dat$mean_opt_temp_c <- ave(s_dat$opt_temp_c, s_dat$common_name)

# Center each size class' optimum temperature by mean optimum temperature for that species
s_dat$opt_temp_c_ct <- s_dat$opt_temp_c - s_dat$mean_opt_temp_c


# B. Growth optimum model ==========================================================
# Plot data
nb.cols <- length(unique(s_dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat, aes(log10(mass_norm), opt_temp_c_ct, 
                  color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = mycolors) +  
  NULL


#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model: lmer(opt_temp_c_ct ~ log10_mass_n + (log10_mass_n | species), s_dat)

# *** dont rely on defaults, specify them for clarity and if defaults change...
m1stanlmer <- stan_lmer(formula = opt_temp_c_ct ~ log10_mass_n + (log10_mass_n | species), 
                        data = s_dat,
                        seed = 8194,
                        adapt_delta = 0.999)

# Summary of priors used
prior_summary(object = m1stanlmer)

# Print model summary
print(m1stanlmer, digits = 3)


#** Check sampling quality & model convergence =====================================
summary(m1stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# Compare the observed distribution of b scores (dark blue line) to 100 simulated datasets from the posterior predictive distribution (light blue lines)
pp_check(m1stanlmer, nreps = 100)
# pp_check(c1_stanlmer, plotfun = "stat_grouped", stat = "median", group = "species")

# Plot Rhat (1 is good)
plot(m1stanlmer, "rhat")

# Plot ess
plot(m1stanlmer, "ess")


#** Extract the posterior draws for all parameters =================================
sims <- as.matrix(m1stanlmer)

para_name <- colnames(sims)

para_name # *** need to know all parameters here (check sigma)

# Create df of parameters I want to plot:
df <- data.frame(species = sort(unique(s_dat$species_ab)))

# Paste parameter name + species to get all species-specific parameters (slope in this case)
df$param <- paste("b[log10_mass_n species:", 
                  sort(unique(df$species_ab)), 
                  "]", 
                  sep = "")

# Extract each species slope and 90% CI
summarym1_90 <- tidy(m1stanlmer, intervals = TRUE, prob =.9, 
                     parameters = "varying")

summarym1_90 <- summarym1_90 %>% 
  filter(term == "log10_mass_n")

# Estimate
df$pred_slope <- summarym1_90$estimate

# 90% CI
df$pred_slope_lwr90 <- summarym1_90$lower
df$pred_slope_upr90 <- summarym1_90$upper

# 50% CI
summarym1_50 <- tidy(m1stanlmer, intervals = TRUE, prob =.5, 
                     parameters = "varying")

summarym1_50 <- summarym1_50 %>% 
  filter(term == "log10_mass_n")

df$pred_slope_lwr50 <- summarym1_50$lower
df$pred_slope_upr50 <- summarym1_50$upper


#** Plot species-slopes ============================================================
#blues <- colorRampPalette(brewer.pal(8, "Blues"))(9)
# pal2 <- viridis(n = 4)
pal <- colorRampPalette(brewer.pal(8, "Paired"))(12)
pal2 <- c("#fdb863", "#b2abd2", "#5e3c99", "#e66101") # From colorbrewer website

# Add overall mean  and 90% & 50% CI
g_mean <- summary(m1stanlmer, probs = c(0.05, 0.95), digits = 5)[2] # Get the slope
g_mean_ci90 <- posterior_interval(m1stanlmer, prob = 0.90, pars = "log10_mass_n")
g_mean_ci50 <- posterior_interval(m1stanlmer, prob = 0.50, pars = "log10_mass_n")

# Scale parameters to make them slope, not diff from mean!
pdat <- df
pdat[, 3:7] <- sweep(df[, 3:7], 2, g_mean, FUN = "+")

# Group species with "stronger" evidence
pdat$sig <- ifelse(pdat$pred_slope < 0 & pdat$pred_slope_upr50 < 0, 
                   2, -9)

pdat$sig <- ifelse(pdat$pred_slope > 0 & pdat$pred_slope_lwr50 > 0, 
                   1, pdat$sig)

pdat$sig <- ifelse(pdat$pred_slope < 0 & pdat$pred_slope_upr90 < 0, 
                   3, pdat$sig)

pdat$sig <- ifelse(pdat$pred_slope > 0 & pdat$pred_slope_lwr90 > 0, 
                   4, pdat$sig)

pdat$sig_f <- as.factor(pdat$sig)

# Group positive and negative slopes
pdat$sign <- ifelse(pdat$pred_slope < 0, "neg", "pos")

# Reorder based on predicted slope and plot!
ggplot(pdat, aes(reorder(species, pred_slope), pred_slope, 
                 shape = sign, color = sig_f)) +
  geom_rect(data = pdat, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci90[1], ymax = g_mean_ci90[2],
            color = "grey93", fill = "grey93") + # overplotting many rectangles here..
  geom_rect(data = pdat, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci50[1], ymax = g_mean_ci50[2],
            color = "grey83", fill = "grey83") + # overplotting many rectangles here..
  #geom_hline(yintercept = 0, col = "grey1", linetype = 1, size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = g_mean, col = pal[8], linetype = "twodash", size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr90, ymax = pred_slope_upr90), 
                color = "gray70", size = 0.5, width = 0) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr50, ymax = pred_slope_upr50), 
                color = "gray50", size = 0.8, width = 0) +
  #scale_color_manual(values = pal[c(4,1,2)]) + # Need to watch out for this, depends on # levels
  #scale_color_manual(values = c("#b2abd2", "#fdb863", "#e66101", "#5e3c99")) + # Need to watch out for this, depends on # levels "pink"
  scale_color_manual(values = pal2) + # Need to watch out for this, depends on # levels "pink"
  geom_point(size = 4, col = "gray50") +
  geom_point(data = filter(pdat, sig > 0), 
             aes(reorder(species, pred_slope), pred_slope, shape = sign, color = sig_f), size = 4) +
  xlab("Species") + 
  ylab("Change in T_opt (C) for growth / order of magnitude change in mass") +
  coord_flip() +
  guides(shape = FALSE, color = FALSE) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  theme(aspect.ratio = 1) +
  #ylim(min(pdat$pred_slope_lwr95), -1*min(pdat$pred_slope_lwr95)) +
  #theme(axis.text.y = element_blank()) + # If I end up with too many species
  NULL

ggsave("figs/growth_species.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#** Plot overall prediction and data ===============================================
fits <- m1stanlmer %>% 
  as_data_frame %>% 
  #rename(intercept = `(Intercept)`) %>% 
  select(-sigma)

glimpse(fits)

reds <- colorRampPalette(brewer.pal(5, "Reds"))(5)

n_draws <- 500
alpha_level <- .1
col_draw <- "grey40"
col_median <- reds[4]

nb.cols <- length(unique(s_dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat) + 
  aes(x = log10_mass_n, y = opt_temp_c_ct, color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = log10_mass_n), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in blue
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$log10_mass_n), 
              size = 1.4, color = col_median) +
  geom_point(size = 4.5) +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "Normalzied log10(mass)",
       y = "Normalized optimum growth temperature") +
  NULL

ggsave("figs/growth_scatter.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")



