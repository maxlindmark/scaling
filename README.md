# Scaling of growth, metabolism and maximum consumption with temperature and body size within species of fish
Project investigating the intraspecific temperature- and size dependence of these rates using data from a systematic literature search.

**Authors:** Max Lindmark, [Jan Ohlberger](http://janohlberger.com/Homepage/), Anna Gårdmark

This repository contains all data and code used for analyses and figures from Lindmark et al (20XX). 

## How to replicate our analyses

### Data
`data`: Data on growth, maximum consumption rate, metabolic rate (oxygen consumption) and optimum temperature for growth as acquired from authors or extracted from figures and tables, including additional biological information of the species, in .xlsx format. .csv files are used for the analysis and contain only key columns. For info on how these data were compiled, see Appendix S1.

`R/exploration`: Code for exploring growth, consumption and metabolism data, filtering data and exporting for model fitting scripts ( as .csv).

### Analysis
`R/analysis`: Code for all analysis (defining models, model selection & validation, prediciton and fit) and corresponding figures.

`R/analysis/JAGS_models`: .R files for generating .txt files of JAGS models for model selection. Split into three folders, one for each model type: log_linear, non_linear and T_opt.

`R/analysis/log_linear`: Scripts for model selection and analysis (validation, fit and prediction) for log_linear models (growth, consumption and metabolism at sub-peak temperature).

`R/analysis/non_linear_cons`: Scripts for polynomial model of consumption data including temperatures beyond peak (model selection, validation, fit and prediction).

`R/analysis/T_opt`: Scripts for linear model of optimum growth temperature as a funciton of body mass (model selection, validation, fit and prediction).

`R/analysis/concept_figure`: Conceptual figure on bioenergetics of growth using parameters estimated in this study.

`R/analysis/growth`: Hierarchical Multiple Regression model in JAGS for mass and temperature dependence of growth at sub-optimum temperatures. Code for specifying JAGS models, perform model selection, validation, predictions and evalute model fit.

`R/analysis/meta_cons`: Hierarchical Multiple Regression model in JAGS for mass and temperature dependence of metabolism and growth at sub-optimum temperatures. Code for specifying JAGS models, perform model selection, validation, predictions and evalute model fit.

`R/analysis/non_linear_model`: Quadratic model in JAGS for mass and temperature dependence of consumption using temperatures that include optimum temperatures. Code for specifying JAGS models, perform model selection, validation, predictions and evalute model fit.

### Figures
`figures`: Figures for main text.

`figures/supp`: Supplementary figures (some are not in paper)


### Key packages
```{r}
# print(sessionInfo())
# other attached packages:
# scales_1.1.0
# MCMCvis_0.14.0
# bayesplot_1.7.1
# patchwork_0.0.1
# viridis_0.5.1
# viridisLite_0.3.0
# magrittr_1.5
# readxl_1.3.1      
# RCurl_1.95-4.12
# bitops_1.0-6
# ggmcmc_1.3
# ggplot2_3.2.1
# tidyr_1.0.0
# dplyr_0.8.3
# RColorBrewer_1.1-2
# rjags_4-10        
# coda_0.19-3    
```
