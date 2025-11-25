# spatstat.model

## Parametric statistical modelling of spatial data for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.model)](http://CRAN.R-project.org/package=spatstat.model) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.model)](https://github.com/spatstat/spatstat.model)

The original `spatstat` package has been split into
several sub-packages
(see [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.model` is one of the
sub-packages. It contains all the main user-level functions that perform
**parametric statistical modelling** of spatial data,
with the exception of data on linear networks.

Most of the functionality is for spatial point patterns in two dimensions.
There is a very modest amount of functionality for 3D and higher dimensional patterns
and space-time patterns.

You are viewing the GitHub repository which holds
the latest **development version** of `spatstat.model`.
For the latest public release on CRAN, click the green badge above.

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
- fitting determinantal point process models to point pattern data (`dppm`)
- fitting recursively partitioned models to point patterns (`rppm`)
- class support for fitted models (`update`, `print`, `summary`, `predict`, `plot`, `simulate`, `coef`, `confint`, `vcov`, `anova`, `residuals`, `fitted`, `deviance`, `AIC`, `logLik`, `terms`, `formula`, `model.matrix`)
- minimum contrast estimation (generic algorithm)
- simulation of fitted point process models

#### Formal inference

- hypothesis tests (quadrat test, Clark-Evans test, Berman test, Diggle-Cressie-Loosmore-Ford test, scan test, studentised permutation test, segregation test, ANOVA tests of fitted models, adjusted composite
likelihood ratio test, envelope tests, Dao-Genton test, balanced independent two-stage test)
- confidence intervals for parameters of a model
- prediction intervals for point counts

#### Informal validation

- residuals
- leverage
- influence
- partial residual plot
- added variable plot
- diagnostic plots
- pseudoscore residual plots
- model compensators of summary functions
- Q-Q plots

### Installing the package

This repository contains the _development version_ of
`spatstat.model`. The easiest way to install the development version
is to start R and type

```R
repo <- c('https://spatstat.r-universe.dev', 'https://cloud.r-project.org')
install.packages("spatstat.model", dependencies=TRUE, repos=repo)
```

To install the latest _public release_ of `spatstat.model`,
type

```R
install.packages("spatstat.model")
```

