# spatstat.model

## Parametric statistical modelling of spatial data for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.model)](http://CRAN.R-project.org/package=spatstat.model) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.model)](https://github.com/spatstat/spatstat.model)

You are viewing the GitHub repository which holds
the latest **development version** of `spatstat.model`.
For the latest public release on CRAN, click the green badge above.

 - [Overview of `spatstat.model`](#overview)
 - [Detailed contents of package](#detailed)
 - [Installing the package](#installing)
 - [Bug reports](#bugreports)
 - [Questions](#questions)
 - [Proposing changes to code](#proposing)
 - [Future development](#future)

___

### <a name="overview"></a> Overview of `spatstat.model`

The original `spatstat` package has been split into
several sub-packages
(see [spatstat/spatstat](https://github.com/spatstat/spatstat))

This package `spatstat.model` is one of the
sub-packages. It contains all the main user-level functions that perform
**parametric statistical modelling** of spatial data,
with the exception of data on linear networks.

Most of the functionality is for spatial point patterns in two dimensions.
There is a very modest amount of functionality for 3D
and higher dimensional patterns and space-time patterns.

`spatstat.model` supports

- parametric modelling (fitting models to point pattern data, model selection, model prediction)
- formal inference (hypothesis tests, confidence intervals)
- informal validation (model diagnostics)

___

### <a name="detailed"></a> Detailed contents of `spatstat.model`

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

___

### <a name="installing"></a> Installing the package

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

___

## <a name="bugreports"></a> Bug reports 

Users are encouraged to report bugs.
If you find a bug in a `spatstat` function,
please identify the sub-package containing that function.
Visit the GitHub repository for the sub-package, 
click the `Issues` tab at the top of the page, 
and press *new issue* to start a new bug report, documentation correction
or feature request.

**Please do not post questions** on the Issues pages,
because they are too clunky for correspondence.

## <a name="questions"></a> Questions about spatstat

For questions about the `spatstat` package family, first check 
the question-and-answer website
[stackoverflow](http://stackoverflow.com/questions/tagged/spatstat)
to see whether your question has already been asked and answered.
If not, you can either post your question at stackoverflow, or
email the authors.

## <a name="proposing"></a> Proposing changes to the code

Feel free to fork `spatstat.model`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

## <a name="future"></a> Future development

The `spatstat` package family is the result of 30 years of software development
and contains over 200,000 lines of code.
It is still under development,
motivated by the needs of researchers in many fields,
and driven by innovations in statistical science.
We welcome contributions of code, and suggestions
for improvements.

