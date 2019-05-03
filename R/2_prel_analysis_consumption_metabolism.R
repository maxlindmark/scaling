#=====================================
# 2019.01.28 First analysis of Cmax 
#=====================================
rm(list=ls(all=TRUE))

##-- install packages
# devtools::install_github("elbersb/tidylog")

.libPaths()
.libPaths("C:/Program Files/R/R-3.5.0/library")
.libPaths()

##-- load libraries
library(lattice)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(lme4)
library(merTools)
library(arm)
library(lmerTest)
library(tidylog)
#display.brewer.all()

##-- set wd
setwd("//storage-og.slu.se/home$/mxli0002/My Documents/Max SLU/Papers/Scaling_exponents_temp/Data_exploration/Consumption")

dat <- read.csv("Consumption_scaling_data.csv", sep = ";")

head(dat)
str(dat)

dat$lm_all <- 1 # Grouping variable to visualize overall slope

##-- Center variables
# Mean temp in experiment
dat <- dat %>% group_by(common_name) %>% mutate(mean_temp_c = mean(temp_c))
dat$temp_c_ct <- dat$temp_c - dat$mean_temp_c

# Relative to temp in nature..
dat$temp_mid_ct <- dat$temp_c - dat$species_temp_mid


#=====================================
# Statistics
#=====================================
# I'm going to go ahead and use all data in the mixed model. A few single-level observations within the random effect wouldn't be too bad, but they contribute to the overall mean effect. https://stats.stackexchange.com/questions/24280/can-i-fit-a-mixed-model-with-subjects-that-only-have-1-observation

df <- dat %>% 
  filter(significant_size == "Y") %>% 
  group_by(common_name) %>% 
  filter(n()>1)

df_full <- dat %>% filter(significant_size == "Y") # do both full and intraspecific only data

ggplot(df, aes(temp_mid_ct, b)) + 
  ylim(0, 1.2) +
  geom_point(size = 5, alpha = 0.7) +
  stat_smooth(method = "lm") +
  theme_classic(base_size = 15) +
  NULL

# basic linear model
m1 <- lm(b ~ temp_mid_ct, data = df)

summary(m1)

par(mfrow = c(2,2))
plot(m1)

par(mfrow = c(1,1))

# check residuals
rdat <- data.frame(residuals = residuals(m1),
                   fitted = predict(m1),
                   temp_mid_ct = m1$model$temp_mid_ct,
                   species = df$species)

# residuals vs fitted
ggplot(rdat, aes(fitted, residuals(m1))) + 
  geom_hline(yintercept = 0) +
  geom_point(size = 4) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# (also by group)
colourCount <- length(unique(rdat$species)) # number of levels
getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))

ggplot(rdat, aes(reorder(species, residuals), residuals, color = species)) + 
  geom_hline(yintercept = 0) +
  theme_classic(base_size = 15) +
  geom_boxplot(size = 1.2) +
  guides(color = F) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# clearly non-independent residuals and they are grouped within species. Try fixed effect independent intercept for each species
m2 <- lm(b ~ temp_mid_ct + species, data = df)

summary(m2)

par(mfrow = c(2,2))
plot(m2)

par(mfrow = c(1,1))
rdat2 <- data.frame(residuals = residuals(m2),
                    temp_mid_ct = m2$model$temp_mid_ct,
                    species = df$species)

ggplot(rdat2, aes(reorder(species, residuals), residuals, color = species)) + 
  geom_hline(yintercept = 0) +
  theme_classic(base_size = 15) +
  geom_boxplot(size = 1.2) +
  guides(color = F) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# N_data / N_parameter
length(df$temp_c)/length(summary(m2)$coefficients[, 1])

# here the between-species variation was accounted for to improve independence of residuals. But we are not really interested in the estimate of each species intercept...

# trying an interaction to also get different slopes for each group
m3 <- lm(b ~ temp_mid_ct * species, data = df)

summary(m3)

par(mfrow = c(2,2))
plot(m3)

par(mfrow = c(1,1))
rdat3 <- data.frame(residuals = residuals(m3),
                    temp_mid_ct = m3$model$temp_mid_ct,
                    species = df$species)

ggplot(rdat3, aes(reorder(species, residuals), residuals, color = species)) + 
  geom_hline(yintercept = 0) +
  theme_classic(base_size = 15) +
  geom_boxplot(size = 1.2) +
  guides(color = F) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  NULL

# N_data / N_parameter
length(df$temp_c)/length(summary(m3)$coefficients[, 1]) # or 23...

# not economic use of data to improve independence of residuals with fixed effects...




#=====================================
# Random effects
#=====================================
# General sources for mixed models:
# lme4 "vignette" https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
# lme4 package: https://cran.r-project.org/web/packages/lme4/lme4.pdf
# GLMM FAQ: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#confidence-intervals-on-conditional-meansblupsrandom-effects
# Nice paper to cite for R2, p-values etc.. https://peerj.com/articles/4794.pdf
# Good simple intro to mixed models including algebraic equations... http://talklab.psy.gla.ac.uk/KeepItMaximalR2.pdf


# Here I could actually use all data because I can if it's not a large minority of levels within groups that have 1 obs... But residuals look worse then so I will wait with that and for now use only intra-specific data.
# https://stats.stackexchange.com/questions/24280/can-i-fit-a-mixed-model-with-subjects-that-only-have-1-observation

df_me <- df

# look at the distribution of fixed intercepts. If somehwat normal, we can model it as a random effect. It's hard to tell but not that far from normal (bit skewed)
m2_coef <- coef(lm(b ~ temp_mid_ct + species - 1, data = df))
hist(m2_coef[-1], xlab = "species intercept estimate")

# I will fit the full model (random intercept and slope) for the fishbase cent temp, but also maybe the normal temp. But the normal temp doesn't make sense since we have now idea how far from the normal temp the experiments were at, and the whole idea is to test effect of that distance!

# following Zuur will start with the full model, and then following Ben Bolker we'll skip model selection and just use what we think is valid. https://gkhajduk.github.io/2017-03-09-mixed-models/

# I use lme4 with lmerTest to get approximate DF and p-values, based on the Satterthwaite approximation. https://www.jstatsoft.org/article/view/v082i13. An alternative would be to bootstrap (see bookmark), or use a likelihood ratio test with and without the main effect *NOTE REML vs ML ESTIMATION*: https://gkhajduk.github.io/2017-03-09-mixed-models/) and also here in the main FAQ: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#is-the-likelihood-ratio-test-reliable-for-mixed-models


##-- Plot data again
# plotting temp as difference between fishbase temp
ggplot(df_me, aes(temp_mid_ct, b)) + 
  geom_point(size = 6, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  NULL

# add in species colour
ggplot(df_me, aes(temp_mid_ct, b, color = common_name)) + 
  geom_point(size = 6, alpha = 0.7) +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", se = F) +
  stat_smooth(aes(group = lm_all), method = "lm", colour = "black", size = 3) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  guides(color = F) +
  NULL

# it seems that the species-effect is on both intercept and slope...

#=====================================
# Fit random intercept and slope model
#=====================================
me1 <- lmer(b ~ temp_mid_ct + (temp_mid_ct | species), df_me)

summary(me1)
coef(me1)

# seems I can also write like this:
# test <- lmer(b ~ temp_mid_ct + (1 + temp_mid_ct | species), df_me)
# summary(test)

#=====================================
# Check importance of random effects
#=====================================

# We can use  the ranova function (note, unlike other Likelihood Ratio Tests, we use the model fitted with REML, not ML): https://stats.stackexchange.com/questions/347368/how-to-interpret-the-output-of-lmertestranova-in-r, See also https://www.rensvandeschoot.com/tutorials/lme4/

ranova(me1, reduce.terms = F) # If TRUE, the model is reduced to random intercept, but we want to contrast it to the only Fixed Effects model. This shows that the model becomes significantlyy worse when dropping the random effect. This is basically a likelihood ration test testing how much more likely the data are in either of the two models tested


#=====================================
# Check assumptions
#=====================================
# https://www.rensvandeschoot.com/tutorials/lme4/

plot(fitted(me1), resid(me1, type = "pearson")) # Residuals look slightly better with the smaller data set
abline(0, 0, col="red")

df_me$fit <- as.vector(fitted(me1)) # Can also use predict(me1) here... This seems to be with all random effects...
df_me$res <- as.vector(resid(me1))

# plot residuals by species
ggplot(df_me, aes(x = reorder(species, res), y = res, color = species)) + 
  geom_hline(yintercept = 0) +
  theme_classic(base_size = 15) +
  geom_boxplot(size = 1.2) +
  guides(color = F) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  NULL

ggplot(df_me, aes(x = fit, y = res, color = species)) + 
  geom_hline(yintercept = 0) +
  theme_classic(base_size = 15) +
  geom_point(size = 5) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  NULL

# QQ-plots (check here for functions and what the plots mean: https://stats.stackexchange.com/questions/362840/formal-definition-of-the-qqline-used-in-a-q-q-plot)
qqnorm(resid(me1)); qqline(resid(me1), col= "red") 
hist(resid(me1))

# so I have very tailed (non-normal) distributions, and there isn't much I can do about that... It will affect p-values and CI's... https://stats.stackexchange.com/questions/125856/r-lmer-model-diagnosis-qqnorm


#=====================================
# Plot model predictions
#=====================================
# first predict population-level effects without any random effects. So it's for a generic species... Note they were used for fitting, so we get a shrinkage and account for that some species have more data etc.https://stats.stackexchange.com/questions/250277/are-mixed-models-useful-as-predictive-models and https://stats.stackexchange.com/questions/262277/why-would-you-predict-from-a-mixed-effect-model-without-including-random-effects

df_me$pred_norand <- predict(me1, re.form = NA) 

confint(me1, method = "Wald") # Here I can get CI for regression coefficients, not the actualt model!! confint() is the same as as 1.96 * SE from the output. This accounts for the random factors: https://stats.stackexchange.com/questions/233800/how-can-i-get-confidence-intervals-for-fixed-effects-using-the-rlmer-function-r/233822#233822 and https://stats.stackexchange.com/questions/117641/how-trustworthy-are-the-confidence-intervals-for-lmer-objects-through-effects-pa.

# Here I'm calculating CI's based on the method in the FAQ: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#lme4, i.e. via the VCV matrix... I am assuming there are no random effects here, because BB uses this line of code together predict above, using the same argument (re.form) to supress random effects.

# The approach is here:
# figure out the model matrix X corresponding to the new data;
# matrix-multiply X by the parameter vector beta to get the predictions (or linear predictor in the case of GLM(M)s);
# extract the variance-covariance matrix of the parameters V
# compute XVX' to get the variance-covariance matrix of the predictions;
# extract the diagonal of this matrix to get variances of predictions;
# if computing prediction rather than confidence intervals, add the residual variance;
# take the square-root of the variances to get the standard deviations (errors) of the predictions;
# compute confidence intervals based on a Normal approximation;
# for GL(M)Ms, run the confidence interval boundaries (not the standard errors) through the inverse-link function.

newdat <- subset(df_me, select = c("species", "temp_mid_ct")) # I should simulate a new data instead...

newdat$b <- predict(me1, newdat, re.form = NA) # Random effect = NA, predicting for a typical species

mm <- model.matrix(terms(me1), newdat)

pvar1 <- diag(mm %*% tcrossprod(vcov(me1), mm))

tvar1 <- pvar1 + VarCorr(me1)$species[1] ## must be adapted for more complex models

cmult <- 1.96

newdat <- data.frame(
  newdat
  , plo = newdat$b - cmult*sqrt(pvar1)
  , phi = newdat$b + cmult*sqrt(pvar1)
  , tlo = newdat$b - cmult*sqrt(tvar1)
  , thi = newdat$b + cmult*sqrt(tvar1)
) # Here plo and phi are the CI's: https://biologyforfun.wordpress.com/2015/06/17/confidence-intervals-for-prediction-in-glmms/

ggplot(newdat, aes(x = temp_mid_ct, y = b)) + 
  geom_point() +
  geom_pointrange(aes(ymin = plo, ymax = phi)) +
  labs(title = "CI based on fixed-effects uncertainty ONLY")

# Now I want to add CI's
df_me$pred_norand_lower <- newdat$plo
df_me$pred_norand_upper <- newdat$phi

ggplot(df_me, aes(temp_mid_ct, b)) + 
  geom_line(data = df_me, aes(temp_mid_ct, pred_norand), size = 3) +
  #geom_ribbon(data = df_me, aes(ymin = pred_norand_lower, ymax = pred_norand_upper, 
  #                              x = temp_mid_ct), alpha = 0.3)+
  geom_point(size = 8, alpha = 0.3) +
  geom_line(data = df_me, aes(temp_mid_ct, pred_norand_upper), size = 3, color = "red", linetype = 2) +
  geom_line(data = df_me, aes(temp_mid_ct, pred_norand_lower), size = 3, color = "red", linetype = 2) +
  theme_classic(base_size = 36) +
  theme(aspect.ratio = 1) +
  labs(x = "Temperature_ct", y = "Size scaling exponent")
NULL

summary(me1)

# These are the results of a "typical" group... i.e. the inter-specific model. i.e. I focus on the population prediction, not accounting for variaiton in the random effect (Bayesian) and assuming there's now skew... https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme

# (I could try using the predictInterval()function just to see they are the same...)

#=====================================
# Intra-specific results
#=====================================
# first check the (point estimates) intercepts and slopes for each random level
# on coef and ranef: 
# https://stats.stackexchange.com/questions/214129/whats-the-interpretation-of-ranef-fixef-and-coef-in-mixed-effects-model-using
# here is why it's Bayesian what we do: https://www.r-bloggers.com/random-regression-coefficients-using-lme4/
# here is a post on the topic: https://stackoverflow.com/questions/31535747/se-for-coef-mermod-coefficients?noredirect=1&lq=1
# which is a dublicate of this, in which ben lists the two options for getting these SE's: https://stackoverflow.com/questions/26198958/extracting-coefficients-and-their-standard-error-from-lme

coef(me1)
se.fixef(me1) # this gives almost exactl the same as the summary output from lmerTest


# create new data frame to store the
pdat <- data.frame(slope   = coef(me1)$species[, 2],
                   se      = as.vector(se.coef(me1)$species[, 2]),
                   species = sort(unique(me1@frame$species)))

# plot distribution (note it's density on y)
ggplot(pdat, aes(slope)) +
  geom_density()


ggplot(pdat, aes(x = slope)) + 
  geom_histogram(aes(y = ..density..), binwidth = 0.002, colour = "black", fill = "grey") +
  geom_density(alpha = 0.2, fill = "darkorange") +
  theme_classic(base_size = 18) + 
  geom_vline(aes(xintercept = mean(slope, na.rm = T)),
             color = "red", linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = 0, color = "steelblue", linetype = "dashed", size = 1.5) +
  NULL

pdat$upr <- pdat$slope + pdat$se * 1.96
pdat$lwr <- pdat$slope - pdat$se * 1.96

ggplot(pdat, aes(reorder(species, slope), slope)) +
  geom_errorbar(aes(reorder(species, slope), ymin = lwr, ymax = upr), color = "grey70", size = 2, width = 0) +
  geom_point(size = 6) +
  guides(color = FALSE) + 
  geom_hline(yintercept = 0, col = "red", linetype = 2, size = 1.7) +
  xlab("Species") + 
  ylab("Slope") +
  coord_flip() +
  theme_classic(base_size = 26) +
  NULL


#=====================================
# Multiple LM's instead
#=====================================
# an alternative is to fit alot of regressions, but I want shrinkage to to mean, borrow information etc. In mixed models, the CI's of these group-level predictions are Bayesian...
# https://stats.stackexchange.com/questions/377687/build-confidence-intervals-for-random-effects-in-intercept-only-models

# compare this plot with the multiple regression approach:

p <- c()
r <- c()
s <- c()
t <- data.frame()

for (i in unique(df_me$species)){
  p <- data.frame(subset(df_me, species == i))  
  
  m1 <- summary(lm(b ~ temp_mid_ct, data = p))
  s <- data.frame(slope = m1$coefficients[2],
                  species = i,
                  p = m1$coefficients[,4][2],
                  se = m1$coefficients[,2][2])
  t <- rbind(t, s)
}

t


# now plot coefficients from species-specific lm's
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

# the estimates are similar, but much more stable CI's, due to borrowing information etc: https://stats.stackexchange.com/questions/377687/build-confidence-intervals-for-random-effects-in-intercept-only-models



#=====================================
# To do
#=====================================
# read this thread: https://stats.stackexchange.com/questions/147836/prediction-interval-for-lmer-mixed-effects-model-in-r

##-- EXTRA NOTES: 

# Here Ben recommends using nlme for simple LMM: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#why-doesnt-lme4-display-denominator-degrees-of-freedomp-values-what-other-options-do-i-have.
# https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode

# Maybe there aren't as many groups (species) as one would hope... read more about that here, and how bayesian models better take care of the uncertainty around those mean and variances estimates
# https://stats.stackexchange.com/questions/37647/what-is-the-minimum-recommended-number-of-groups-for-a-random-effects-factor

# See also here for how Bayesian models can account for different error structures etc
# https://stats.stackexchange.com/questions/64226/lme-and-lmer-comparison

# on R^2: https://peerj.com/articles/4794.pdf

# only random slope (this makes the random effect not significant in ranova())
# me2 <- lmer(b ~ temp_mid_ct + (0 + temp_mid_ct | species), df_me)


