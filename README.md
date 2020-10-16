# SurvBin Package

## Description

The **SurvBin** package contains four key functions: 

  - *lstats* to compute the standardized L-statistic  (with pooled/unpooled variance estimators);
  - *lstats_boots* to compute the standardized L-statistic (with bootstrap variance estimator);
  - *survtest* to compute the Weighted Kaplan-Meier based-test;
  - *bintest* to compute the proportions' test. 

The **SurvBin** package also provides the function *survbinCov* that can be used to calculate the covariance between the survival and binary statistics; and the function *simsurvbin* for simulating bivariate binary and survival data.

More information can be found at: https://martabofillroig.github.io/SurvBin/

## Installation

``` r
# install.packages("devtools")
library(devtools)
install_github("MartaBofillRoig/SurvBin")
```

## References

This repository also contains the source files of the preprint:

- "A class of two-sample nonparametric statistics for binary and time-to-event outcomes". Marta Bofill Roig, Guadalupe Gómez Melis. (2020). 
https://arxiv.org/abs/2002.01369

Specifically, in the folder CODE_paper/CaseStudy, there is the source code for reproduce the illustration; and in the folder CODE_paper/Extension_Simulation, there is the code to reproduce the simulation study (see Extension_Simulation subfolder).
