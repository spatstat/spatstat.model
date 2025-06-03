\name{dropROC}
\alias{dropROC}
\title{
  ROC Curves for all Single Term Deletions from a Model
}
\description{
  Given a fitted point process model,
  consider dropping each possible term in the model,
  and compute the ROC curve for the dropped covariate.
}
\usage{
dropROC(object, scope = NULL, high=TRUE, \dots)
}
\arguments{
  \item{object}{
    A fitted point process model (object of class \code{"ppm"},
    \code{"kppm"}, \code{"dppm"}, \code{"slrm"} or \code{"lppm"}).
  }
  \item{scope}{
    A formula or a character vector specifying the terms to be
    considered for deletion. The default is all possible terms.
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
  This function is like \code{\link[stats]{drop1}}
  in that it considers each possible term in the model \code{object}
  (or only the terms listed in the \code{scope} argument),
  deletes each such term from the model, and measures the change in the model.
  In this case the change is measured by computing the ROC curve
  for the deleted covariate, using the updated model as a baseline.
}
\value{
  A named list containing the ROC curves for each possible deletion.
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
  \code{\link[spatstat.model]{dropply}},
  \code{\link[spatstat.model]{addROC}},
  \code{\link[spatstat.model]{roc.ppm}}.
}
\examples{
  dimyx <- if(interactive()) NULL else 32
  fut <- ppm(bei ~ grad + elev, data=bei.extra)
  z <- dropROC(fut, dimyx=dimyx)
  plot(z)

  ## how to compute AUC for each 
  sapply(z, auc)
}
\keyword{spatial}
\keyword{nonparametric}
