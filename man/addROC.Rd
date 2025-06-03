\name{addROC}
\alias{addROC}
\title{
  ROC Curves for Single Term Additions to a Model
}
\description{
  Given a fitted point process model,
  consider adding new explanatory variables,
  and compute the ROC curve for each new variable.
}
\usage{
addROC(object, scope, high=TRUE, \dots)
}
\arguments{
  \item{object}{
    A fitted point process model (object of class \code{"ppm"},
    \code{"kppm"}, \code{"dppm"}, \code{"slrm"} or \code{"lppm"})
    specifying the model to be extended. 
  }
  \item{scope}{
    A formula or a character vector specifying the variable or variables to be
    considered for addition, or a fitted point process model containing
    all of these variables.
  }
  \item{high}{
    Argument passed to \code{\link[spatstat.explore]{roc}} to specify whether the
    ROC curves should be based on high or low values of the covariates.
    Either a logical value, or a logical vector of the same length
    as \code{scope}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} to control the
    pixel resolution used for calculation.
  }
}
\details{
  This function is like \code{\link[stats]{add1}}
  in that it considers each possible term that could be
  added to the model \code{object}
  (or only the terms listed in the \code{scope} argument),
  adds each such term to the model, and measures the change in the model.
  In this case the change is measured by computing the ROC curve
  for the added covariate, using the original model \code{object} as a
  baseline.

  Either \code{object} or \code{scope} should be a fitted point process
  model, and the other argument may be a fitted point process model or a
  formula. If \code{object} is a fitted model then \code{scope} may be a
  character vector of the names of variables to be added.
}
\value{
  A named list containing the ROC curves for each new explanatory variable.
  The individual entries belong to class \code{"fv"},
  so they can be plotted.
  The list belongs to the class \code{"anylist"}
  so it can be plotted in its entirety.
}
\author{
  \rocketAuthors.
}
\references{
  \rocpaper
}
\seealso{
  \code{\link[spatstat.model]{addapply}},
  \code{\link[spatstat.model]{roc.ppm}}.
}
\examples{
  dimyx <- if(interactive()) NULL else 32
  fit0 <- ppm(bei ~ 1, data=bei.extra)

  z <- addROC(fit0, . ~ grad + elev, dimyx=dimyx)

  plot(z)

  ## how to compute AUC for each curve
  sapply(z, auc)
}
\keyword{spatial}
\keyword{nonparametric}
