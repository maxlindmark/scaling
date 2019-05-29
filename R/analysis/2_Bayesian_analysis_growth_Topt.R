#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.23: Max Lindmark
#
# - Analyze growth rate data over temperature and body size
# 
# A. Load libraries & read data
#
# B. Fit model for T_opt ~ temperature
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
# We use geometric mean for size, and if this is not possible we'll go with mid point 
# of size range
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
#dat$mass_norm <- dat$mass / dat$w_max_published_g

# Intraspecific data (n>1 species)
s_dat <- data.frame(
  dat %>% 
  group_by(common_name) %>% 
  filter(n()>1)
)

# Log mass
s_dat$log10_mass <- log10(s_dat$mass)

# Log mass normalized
#s_dat$log10_mass_n <- log10(s_dat$mass_norm)

# Calculate mean optimum temperature within species
s_dat$mean_opt_temp_c <- ave(s_dat$opt_temp_c, s_dat$common_name)

# Center each size class' optimum temperature by mean optimum temperature for that species
s_dat$opt_temp_c_ct <- s_dat$opt_temp_c - s_dat$mean_opt_temp_c


# B. Growth optimum model ==========================================================
# Plot data
nb.cols <- length(unique(s_dat$species))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(s_dat, aes(log10(mass), opt_temp_c_ct, 
                  color = common_name)) + 
  geom_point(size = 5) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = mycolors) +  
  NULL


#** Set up stan model ==============================================================
# The rstanarm package uses lme4 syntax. In a preliminary analysis I explored a random intercept and slope model: lmer(opt_temp_c_ct ~ log10_mass_n + (log10_mass_n | species), s_dat)

# *** dont rely on defaults, specify them for clarity and if defaults change...
m1stanlmer <- stan_lmer(formula = opt_temp_c_ct ~ log10_mass + (log10_mass | species), 
                        data = s_dat,
                        seed = 8194,
                        adapt_delta = 0.99)

# Summary of priors used
prior_summary(object = m1stanlmer)

# Print model summary
print(m1stanlmer, digits = 3)


#** Check sampling quality & model convergence =====================================
summary(m1stanlmer, 
        probs = c(0.025, 0.975),
        digits = 5)


#** Extract the posterior draws for all parameters =================================
# Create df of parameters I want to plot:
df <- data.frame(species = sort(unique(s_dat$species_ab)))

# Paste parameter name + species to get all species-specific parameters (slope in this case)
df$param <- paste("b[log10_mass species:", 
                  sort(unique(df$species_ab)), 
                  "]", 
                  sep = "")

# Extract each species slope and 90% CI
summarym1_90 <- tidy(m1stanlmer, intervals = TRUE, prob =.9, 
                     parameters = "varying")

summarym1_90 <- summarym1_90 %>% 
  filter(term == "log10_mass")

# Estimate
df$pred_slope <- summarym1_90$estimate

# 90% CI
df$pred_slope_lwr90 <- summarym1_90$lower
df$pred_slope_upr90 <- summarym1_90$upper

# 50% CI
summarym1_50 <- tidy(m1stanlmer, intervals = TRUE, prob =.5, 
                     parameters = "varying")

summarym1_50 <- summarym1_50 %>% 
  filter(term == "log10_mass")

df$pred_slope_lwr50 <- summarym1_50$lower
df$pred_slope_upr50 <- summarym1_50$upper


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
  aes(x = log10_mass, y = opt_temp_c_ct, color = species) + 
  # Plot a random sample of rows as gray semi-transparent lines
  geom_abline(aes(intercept = `(Intercept)`, 
                  slope = log10_mass), 
              data = sample_n(fits, n_draws), color = col_draw, 
              alpha = alpha_level, size = 1.1) + 
  # Plot the median values in blue
  geom_abline(intercept = median(fits$`(Intercept)`), 
              slope = median(fits$log10_mass), 
              size = 1.4, color = col_median) +
  geom_point(size = 4.5) +
  theme_classic(base_size = 15) +
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = TRUE) +
  guides(color = FALSE) +
  annotate(geom = "text", x = 0, y = -2.5, 
           label = paste("y=", round(median(fits$`(Intercept)`), digits = 2),
                         round(median(fits$log10_mass), digits = 2), "x", sep=""),
           size = 5.5, fontface = "italic") +
  theme_classic(base_size = 17) +
  theme(aspect.ratio = 4/5) +
  labs(x = "log10(mass)",
       y = "Normalized optimum growth temperature") +
  NULL

#ggsave("figs/growth_scatter.pdf", plot = last_plot(), scale = 1, width = 18, height = 18, units = "cm")


#** Plot posterior for slope ===============================================
# --- not finished
color_scheme_set("teal")
g_slope_post <- as.array(m1stanlmer)
dimnames(g_slope_post)

mcmc_dens(g_slope_post,
          pars = c("log10_mass"),
          probs = 0.8) +
  theme_classic(base_size = 13) +
  theme(aspect.ratio = 4/5,
        plot.title = element_text(hjust = 0.5, size = 13)) +
  theme() +
  NULL
