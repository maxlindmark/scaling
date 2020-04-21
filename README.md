# Scaling of growth, metabolism and maximum consumption with temperature and body size within species of fish
Project investigating the intraspecific temperature- and size dependence of these rates using data from a systematic literature search.

**Authors:** Max Lindmark, Jan Ohlberger, Anna Gårdmark

This repository contains all of the clean data and code used to generate the analyses and figures from Lindmark et al. 

## How to replicate our analyses

### Data
`data`: Contains data (.xlsx) on metabolic rate, maximum consumption rate, growth rate and optimum temperature for growth. For info on how these data were compiled, see Appendix SX.

`R/exploration`: Contains code for exploring data, filtering data and exporting for model fitting scripts ( as .csv)

### Analysis
`R/analysis`: Contains code for all analysis

`R/analysis/concept_figure`: Conceptual figure on bioenergetics of growth

`R/analysis/growth`: Hierarchical Multiple Regression model in JAGS for growth~temperature and mass using sub-optimum temperatures only. Contains code for specifying JAGS models, perform model selection, perform model validation, make predictions and evalute model fit

`R/analysis/meta_cons`: Hierarchical Multiple Regression model in JAGS for metabolism or consumprtion ~ temperature and mass using sub-optimum temperatures only. Contains code for specifying JAGS models, perform model selection, perform model validation, make predictions and evalute model fit

`R/analysis/non_linear_model`: Quadratic model in JAGS for consumprtion ~ temperature and mass using temperatures that include optimum temperatures. Contains code for specifying JAGS models, perform model selection, perform model validation, make predictions and evalute model fit

### Figures
`figures`: Contains figures for main text

`figures/supp`: Contains supplementary figures


### Key packages
```{r}
# other attached packages:
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
