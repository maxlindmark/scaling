# Scaling of growth, metabolism and maximum consumption with temperature and body size within species of fish
Project investigating the intraspecific temperature- and size dependence of these rates using data from a systematic literature search.

**Authors:** [Max Lindmark](https://maxlindmark.netlify.app/), [Jan Ohlberger](http://janohlberger.com/Homepage/), [Anna Gårdmark](https://internt.slu.se/en/cv-originals/anna-gardmark/)

This repository contains all data and code used for analyses and figures from Lindmark et al (20XX). 

## How to replicate our analyses

### Data
`data`: Data on growth, maximum consumption rate, metabolic rate (oxygen consumption) and optimum temperature for growth as acquired from authors or extracted from figures and tables, including additional biological information of the species, in .xlsx format. .csv files are used for the analysis and contain only key columns. For info on how these data were compiled, see below and Appendix S1.

`R/exploration`: Code for exploring growth, consumption and metabolism data, filtering data and exporting for model fitting scripts ( as .csv).

`data/lists_of_papers_from_lit_search`: .txt files containing all articles found through the literature search, before applying selection criteria.

### Analysis
`R/analysis`: Code for all analysis (defining models, model selection & validation, prediciton and fit) and corresponding figures.

`JAGS_models`: .R files for generating .txt files of JAGS models for model selection. Split into three folders, one for each model type: log_linear, non_linear and T_opt.

`R/analysis/log_linear`: Scripts for model selection and analysis (validation, fit and prediction) for log_linear models (growth, consumption and metabolism at sub-peak temperature).

`R/analysis/unimodal_consumption`: Scripts for fitting the Sharpe-Schoolfield equation to consumption data including temperatures beyond peak (model validation, fit and prediction).

`R/analysis/T_opt`: Scripts for linear model of optimum growth temperature as a funciton of body mass (model selection, validation, fit and prediction).

`R/analysis/concept_model`: Conceptual figure on bioenergetics of growth using parameters estimated in this study for maximum consumption rate and metabolic rate, converted to common unit (g/g/d or g/d), rescaled and standardized for plotting purposes.

`R/analysis/log_linear/growth`: Hierarchical Multiple Regression model in JAGS for mass and temperature dependence of growth at sub-optimum temperatures. Code for model selection, validation, predictions and evalute model fit.

`R/analysis/log_linear/meta_cons`: Hierarchical Multiple Regression model in JAGS for mass and temperature dependence of metabolism and growth at sub-optimum temperatures. Code for model selection, validation, predictions and evalute model fit.

### Figures
`figures`: Figures for main text.

`figures/supp`: Supplementary figures (some are not in paper nor supplement)


### Key packages
```{r}
# print(sessionInfo())
# other attached packages:
# scales_1.1.1
# MCMCvis_0.14.0
# bayesplot_1.7.2
# patchwork_1.0.1
# viridis_0.5.1
# viridisLite_0.3.0
# magrittr_2.0.1
# readxl_1.3.1      
# RCurl_1.98-1.2
# bitops_1.0-6
# ggmcmc_1.4.1
# ggplot2_3.2.2
# tidyr_1.1.2
# dplyr_1.0.2
# RColorBrewer_1.1-2
# rjags_4-10        
# coda_0.19-4    
# tidylog_1.0.2   
```


