# spatstat.model

## Parametric statistical modelling of spatial data for the spatstat family

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spatstat.model)](http://CRAN.R-project.org/package=spatstat.model) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.model)](https://github.com/spatstat/spatstat.model)

The original `spatstat` package has been split into
several sub-packages (See [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.model` is one of the
sub-packages. It contains all the main user-level functions that perform
**parametric statistical modelling** of spatial data,
with the exception of data on linear networks.

Most of the functionality is for spatial point patterns in two dimensions.
There is a very modest amount of functionality for 3D and higher dimensional patterns
and space-time patterns.

### Overview 

`spatstat.model` supports

- parametric modelling (fitting models to point pattern data, model selection, model prediction)
- formal inference (hypothesis tests, confidence intervals)
- informal validation (model diagnostics)

### Detailed contents

For a full list of functions, see the help file for `spatstat.model-package`.

#### Parametric modelling 
- fitting Poisson point process models to point pattern data (`ppm`)
- fitting spatial logistic regression models to point pattern data (`slrm`)
- fitting Cox point process models to point pattern data (`kppm`)
- fitting Neyman-Scott cluster process models to point pattern data (`kppm`)
- fitting Gibbs point process models to point pattern data (`ppm`)
- fitting recursively partitioned models to point patterns (`rppm`)
- class support for fitted models (`update`, `summary`, `predict`, `plot`, `simulate`, `coef`, `confint`, `vcov`, `anova`)
- minimum contrast estimation
- simulation of fitted point process models

#### Formal inference

- hypothesis tests (quadrat test, Clark-Evans test, Berman test, Diggle-Cressie-Loosmore-Ford test, scan test, studentised permutation test, segregation test, ANOVA tests of fitted models, adjusted composite
likelihood ratio test, envelope tests, Dao-Genton test, balanced independent two-stage test)
- confidence intervals for parameters of a model
- prediction intervals for point counts

#### Informal validation

- leverage
- influence
- partial residuals
- added variable plot
- diagnostic plots
- pseudoscore residual plots
- model compensators
- Q-Q plots

