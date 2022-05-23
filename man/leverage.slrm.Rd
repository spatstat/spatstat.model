\name{leverage.slrm}
\alias{leverage.slrm}
\alias{influence.slrm}
\alias{dfbetas.slrm}
\alias{dffit.slrm}
\title{
  Leverage and Influence Diagnostics for Spatial Logistic Regression
}
\description{
  For a fitted spatial logistic regression model, these functions
  compute diagnostics of leverage and influence.
}
\usage{
\method{leverage}{slrm}(model, \dots)
\method{influence}{slrm}(model, \dots)
\method{dfbetas}{slrm}(model, \dots)
\method{dffit}{slrm}(object, \dots)
}
\arguments{
  \item{model,object}{
    A fitted spatial logistic regression model
    (object of class \code{"slrm"}).
  }
  \item{\dots}{
    Arguments passed to methods.
  }
}
\details{
  These functions are methods for the generics
  \code{\link{leverage}}, \code{\link[stats]{influence}},
  \code{\link[stats]{dfbetas}} and \code{\link{dffit}}
  for the class \code{"slrm"}.
  
  These functions adapt the standard diagnostics for logistic regression
  (see \code{\link[stats]{influence.measures}}) to a
  fitted spatial logistic regression model (object of class
  \code{"slrm"}). This adaptation was described by Baddeley, Chang
  and Song (2013).

  \code{leverage.slrm} computes the leverage value (diagonal of the hat
  matrix) for the covariate data in each pixel.
  The result is a pixel image.

  \code{influence.slrm} computes the likelihood influence 
  for the data (covariates and presence/absence of points) in each pixel.
  The result is a pixel image.

  \code{dfbetas.slrm} computes the parameter influence 
  for the data (covariates and presence/absence of points) in each pixel.
  The result is a list of pixel images, one image for each of the
  model coefficients in \code{coef(model)}. The list can be plotted
  immediately.

  \code{dffit.slrm} computes the total influence 
  for the data (covariates and presence/absence of points) in each pixel.
  The result is a pixel image.
}
\value{
  A pixel image, or a list of pixel images.
}
\references{
  Baddeley, A., Chang, Y.M. and Song, Y. (2013)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 86--104.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats]{influence.measures}}.
  
  \code{\link{leverage.ppm}},
  \code{\link{influence.ppm}},
  \code{\link{dfbetas.ppm}},
  \code{\link{dffit.ppm}}
}
\examples{
   H <- unmark(humberside)
   fit <- slrm(H ~ x+y, dimyx=32)
   plot(leverage(fit))
   plot(influence(fit))
   plot(dfbetas(fit))
   plot(dffit(fit))
}
\keyword{spatial}
\keyword{models}
\concept{diagnostics}
