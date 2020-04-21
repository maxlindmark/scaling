# Scaling of growth, metabolism and maximum consumption with temperature and body size within species of fish
Project investigating the intraspecific temperature- and size dependence of these rates using data from a systematic literature search.

**Authors:** Max Lindmark, Jan Ohlberger, Anna Gårdmark

This repository contains all of the clean data and code used to generate the analyses and figures from Lindmark et al. 

## How to replicate our analyses

### Data
'data': Contains data (.xlsx) on metabolic rate, maximum consumption rate, growth rate and optimum temperature for growth. For info on how these data were compiled, see Appendix SX.

**R/exploration**: Contains code for exploring data and filtering data for modelling

### Analysis
**R/analysis**: Contains code for all analysis

**R/analysis/concept_figure**: Conceptual figure on bioenergetics of growth

**R/analysis/growth**: Hierarchical Multiple Regression model in JAGS for growth~temperature and mass using sub-optimum temperatures only. Contains code for specifying JAGS models, perform model selection, perform model validation, make predictions and evalute model fit

**R/analysis/meta_cons**: Hierarchical Multiple Regression model in JAGS for metabolism or consumprtion ~ temperature and mass using sub-optimum temperatures only. Contains code for specifying JAGS models, perform model selection, perform model validation, make predictions and evalute model fit

**R/analysis/non_linear_model**: Quadratic model in JAGS for consumprtion ~ temperature and mass using temperatures that include optimum temperatures. Contains code for specifying JAGS models, perform model selection, perform model validation, make predictions and evalute model fit

### Figures
**figures**: Contains figures for main text

**figures/supp**: Contains supplementary figures
