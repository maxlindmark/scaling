# Optimum growth temperature declines with body size within fish species

**Authors:** [Max Lindmark](https://maxlindmark.netlify.app/), [Jan Ohlberger](http://janohlberger.com/Homepage/), [Anna Gårdmark](https://internt.slu.se/en/cv-originals/anna-gardmark/)

This repository contains all data and code used for analyses and figures from [Lindmark et al (2021)](https://www.biorxiv.org/content/10.1101/2021.01.21.427580v2).

Citation: Lindmark, M., Ohlberger, J. & Gårdmark, A. (2021). [Optimum growth temperature declines with body size within fish species](https://www.biorxiv.org/content/10.1101/2021.01.21.427580v2). bioRxiv, 2021.01.21.427580.

## How to replicate our analyses

### Data
`data`: Data on growth, maximum consumption rate, metabolic rate (oxygen consumption) and optimum temperature for growth as acquired from authors or extracted from figures and tables, including additional biological information of the species, in .xlsx format. .csv files are used for the analysis and contain only key columns. For info on how these data were compiled, see below and Appendix S1.

`R/exploration`: Code for exploring growth, consumption and metabolism data, filtering data and exporting data for model fitting scripts (as .csv files). 

`data/lists_of_papers_from_lit_search`: .txt files containing all articles found through the literature search, before applying selection criteria.

### Analysis
`JAGS_models`: .R files for generating .txt files of JAGS models for model selection. Split into three folders, one for each model type: log_linear, unimodal_consumption and T_opt.

`R/analysis`: Code for all analysis (defining models, model selection & validation, prediction and fit) and corresponding figures.

`R/analysis/log_linear`: Scripts for model selection and analysis (validation, fit and prediction) for log_linear models (consumption and metabolism at sub-peak temperature).

`R/analysis/T_opt`: Scripts for linear model of optimum growth temperature as a function of log body mass (model selection, validation, fit and prediction).

`R/analysis/unimodal_consumption`: Scripts for fitting the Sharpe-Schoolfield equation to consumption data including temperatures beyond peak (model validation, fit and prediction).

`R/analysis/concept_model.R`: Conceptual figure on bioenergetics of growth using parameters estimated in this study for maximum consumption rate and metabolic rate, converted to common unit (g/g/d or g/d), rescaled and standardized for plotting purposes.

### Figures
`figures`: Figures for main text.

`figures/supp`: Supplementary figures (some are not in paper nor supplement)


### Key packages
```{r}
# print(sessionInfo())
# other attached packages:
# broom_0.7.10
# nls.multstart_1.2.0
# rTPC_1.0.0
# tidylog_1.0.2
# scales_1.1.1
# MCMCvis_0.14.0
# bayesplot_1.7.2
# patchwork_1.1.1
# viridis_0.5.1
# viridisLite_0.4.0
# magrittr_2.0.1
# readxl_1.3.1       
# RCurl_1.98-1.5
# ggmcmc_1.4.1
# ggplot2_3.3.5
# tidyr_1.1.4
# dplyr_1.0.7
# RColorBrewer_1.1-2
# rjags_4-10
# coda_0.19-4     
```


