#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.11: Max Lindmark
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
          "RColorBrewer")

library(viridis)

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/package_info.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
pkg_info(pkgs)

# package loadedversion
# 1         broom         0.4.4
# 2         dplyr       0.8.0.1
# 3       ggplot2         3.1.1
# 4          lme4        1.1-19
# 5      magrittr           1.5
# 6        mlmRev         1.0-7
# 7          plyr         1.8.4
# 8  RColorBrewer         1.1-2
# 9         RCurl     1.95-4.12
# 10       readxl         1.3.1
# 11     rstanarm        2.18.2
# 12      tidylog         0.1.0

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

unique(dat$species)

# Normalize variables
# Inspect temperatures
ggplot(dat, aes(env_temp_mid)) +
  geom_histogram() + 
  facet_wrap(~rate)

# Relative to mid temp in environment
dat$env_temp_mid_norm <- dat$temp_c - dat$env_temp_mid

# Filter intraspecific data consumption
s_con <- dat %>% 
  filter(significant_size == "Y" & rate == "consumption") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# Filter intraspecific data metabolism
s_met <- met %>% 
  filter(significant_size == "Y"& rate == "metabolism") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

# If env temp == NA, use mid point of pref.

#======== B. CONSUMPTION ================================================
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
c1_stanlmer <- stan_lmer(formula = b ~ env_temp_mid_norm + (env_temp_mid_norm | species), 
                         data = s_con,
                         seed = 349)

# Summary of priors used
prior_summary(object = c1_stanlmer)

# Print model summary
print(c1_stanlmer, digits = 3)


#====**** Check sampling quality & model convergence ===============================
summary(c1_stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)

# *** What is difference from the previous mixed model? Which species have I removed? Fit mixed model and look at overall and species slope (not their CI, which is hard)
# *** The exploring-data set should end with the selection I will use! Did I forget to list potentially wrong species?

# Compare the observed distribution of b scores (dark blue line) to 100 simulated datasets from the posterior predictive distribution (light blue lines)
pp_check(c1_stanlmer, nreps = 100)
# pp_check(c1_stanlmer, plotfun = "stat_grouped", stat = "median", group = "species")

# Plot Rhat (1 is good)
plot(c1_stanlmer, "rhat")

# Plot ess
plot(c1_stanlmer, "ess")


#====**** Extract the posterior draws for all parameters ===========================
sims <- as.matrix(c1_stanlmer)

para_name <- colnames(sims)

para_name # *** need to know all parameters here (check sigma)

# Create df of parameters I want to plot:
species_b <- data.frame(species = sort(unique(s_con$species)))

species_b$param <- paste("b[env_temp_mid_norm species:", 
                         sort(unique(species_b$species)), 
                         "]", 
                         sep = "")

# For now need to remove this species because I haven't put in a temperature and NAs are removed in Rstan
species_b <- species_b[-which(species_b$species == "Epinephelus coioides"), ]

# Extract each species slope and 95% CI
summaryc1_95 <- tidy(c1_stanlmer, intervals = TRUE, prob =.95, 
                     parameters = "varying")

summaryc1_95 <- summaryc1_95 %>% 
  filter(term == "env_temp_mid_norm")

species_b$pred_slope <- summaryc1_95$estimate
species_b$pred_slope_lwr95 <- summaryc1_95$lower
species_b$pred_slope_upr95 <- summaryc1_95$upper

# Extract each species slope and 90% CI
summaryc1_90 <- tidy(c1_stanlmer, intervals = TRUE, prob =.9, 
                     parameters = "varying")

summaryc1_90 <- summaryc1_90 %>% 
  filter(term == "env_temp_mid_norm")

species_b$pred_slope <- summaryc1_95$estimate
species_b$pred_slope_lwr90 <- summaryc1_90$lower
species_b$pred_slope_upr90 <- summaryc1_90$upper

# Extract each species slope and 50% CI
summaryc1_50 <- tidy(c1_stanlmer, intervals = TRUE, prob =.5, 
                     parameters = "varying")

summaryc1_50 <- summaryc1_50 %>% 
  filter(term == "env_temp_mid_norm")

species_b$pred_slope <- summaryc1_95$estimate
species_b$pred_slope_lwr50 <- summaryc1_50$lower
species_b$pred_slope_upr50 <- summaryc1_50$upper

#====**** Plot species-slopes ==================================================
#blues <- colorRampPalette(brewer.pal(8, "Blues"))(9)
pal <- colorRampPalette(brewer.pal(8, "Set1"))(10)
pal2 <- viridis(n = 4)

# Add overall mean 
g_mean <- summary(c1_stanlmer, probs = c(0.025, 0.975), digits = 5)[2] # Get the slope
g_mean_ci95 <- posterior_interval(c1_stanlmer, prob = 0.95, pars = "env_temp_mid_norm")

# Scale parameters to make them slope, not diff from mean!
pdat <- species_b
pdat[, 3:9] <- sweep(species_b[, 3:9], 2, g_mean, FUN = "+")

# Group species with "stronger" evidence
pdat$sig <- ifelse(pdat$pred_slope < 0 & pdat$pred_slope_upr50 < 0, 
                        2, -9)
pdat$sig <- ifelse(pdat$pred_slope > 0 & pdat$pred_slope_lwr50 > 0, 
                        1, pdat$sig )

# Group positive and negative slopes
pdat$sign <- ifelse(pdat$pred_slope < 0, "neg", "pos")

# Reorder based on predicted slope
ggplot(pdat, aes(reorder(species, pred_slope), pred_slope, shape = sign)) +
  geom_rect(data = species_b, aes(reorder(species, pred_slope), pred_slope),
            xmin = 0, xmax = 100, ymin = g_mean_ci95[1], ymax = g_mean_ci95[2],
            color = "grey93", fill = "grey93") + # overplotting many rectangles here..
  geom_hline(yintercept = 0, col = "grey30", linetype = 1, size = 0.8, alpha = 0.5) +
  geom_hline(yintercept = g_mean, col = pal[5], linetype = "twodash", size = 0.8, alpha = 0.7) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr95, ymax = pred_slope_upr95), 
                color = "gray20", size = 2, width = 0, alpha = 0.2) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr90, ymax = pred_slope_upr90), 
                color = "gray20", size = 2, width = 0, alpha = 0.3) +
  geom_errorbar(aes(reorder(species, pred_slope), 
                    ymin = pred_slope_lwr50, ymax = pred_slope_upr50), 
                color = "gray20", size = 2, width = 0, alpha = 0.4) +
  scale_fill_manual(values = pal[c(2, 1)]) +
  scale_shape_manual(values = c(21, 23)) + 
  geom_point(size = 2, fill = "gray30", alpha = 0.8, color = "white") +
  geom_point(data = filter(pdat, sig > 0), 
             aes(reorder(species, pred_slope), pred_slope, shape = sign, fill = factor(sig)), 
             size = 2, alpha = 1, color = "white") +
  xlab("Species") + 
  ylab("Slope") +
  coord_flip() +
  guides(color = FALSE, shape = FALSE, fill = FALSE) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 8, face = "italic")) +
  theme(aspect.ratio = 1) +
  ylim(min(pdat$pred_slope_lwr95), -1*min(pdat$pred_slope_lwr95)) +
  # theme(axis.text.y = element_blank()) + # If I end up with too many species
  NULL

ggsave("test.pdf", plot = last_plot(), scale = 1, width = 16, height = 16, units = "cm")


#====**** Plot intercept (exponent at mid temp) ====================================
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
             size = 1.5, color = "gray20")
  xlab("Intercept") +
  theme(aspect.ratio = 1) +
  NULL
  
# *** Need to check the arguments of this function...
pp_check(c1_stanlmer) + 
  theme(aspect.ratio = 1) +
  theme_classic(base_size = 15)


#====**** Plot overall prediction and data ====================================
# Plot data again
nb.cols <- length(unique(s_con$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

a <- summary(c1_stanlmer, probs = c(0.025, 0.975), digits = 5)[1]
b <- summary(c1_stanlmer, probs = c(0.025, 0.975), digits = 5)[2]

ggplot(s_con, aes(env_temp_mid_norm, b, color = species)) + 
  geom_point(size = 5) +
  geom_abline(intercept = a, slope = b, size = 2) +
  theme_classic(base_size = 15) +
  guides(color = FALSE) +
  scale_color_manual(values = mycolors)
NULL





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



















#------------------------------------
# Mixed model


me1 <- lmer(b ~ env_temp_mid_norm + (env_temp_mid_norm | species), s_con)
summary(me1)
coef(me1)
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






#------------------------------------

#======== B. METABOLISM =================================================


