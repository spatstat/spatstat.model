\name{methods.dppm}
\alias{methods.dppm} %DoNotExport
\alias{coef.dppm}
\alias{formula.dppm}
\alias{print.dppm}
\alias{terms.dppm}
\alias{labels.dppm}
\title{
  Methods for Determinantal Point Process Models
}
\description{
  These are methods for the class \code{"dppm"}.
}
\usage{
\method{coef}{dppm}(object, \dots)
\method{formula}{dppm}(x, \dots)
\method{print}{dppm}(x, ...)
\method{terms}{dppm}(x, \dots)
\method{labels}{dppm}(object, \dots)
}
\arguments{
  \item{x,object}{
    An object of class \code{"dppm"},
    representing a fitted determinantal point process model.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link[stats]{coef}},
  \code{\link[stats]{formula}},
  \code{\link[base]{print}},
  \code{\link[stats]{terms}} and 
  \code{\link[base]{labels}}
  for the class \code{"dppm"}.

  An object of class \code{"dppm"} represents a fitted
  determinantal point process model. It is obtained from \code{\link{dppm}}.

  The method \code{coef.dppm} returns the vector of
  \emph{regression coefficients} of the fitted model.
  It does not return the interaction parameters.
}
\value{
  See the help files for the corresponding generic functions.
}
\author{
  \adrian
}
\seealso{
  \code{\link{dppm}},  \code{\link{plot.dppm}},
  \code{\link{predict.dppm}}, \code{\link{simulate.dppm}},
  \code{\link{as.ppm.dppm}}.
}
\examples{
  fit <- dppm(swedishpines ~ x + y, dppGauss, method="c")
  coef(fit)
  formula(fit)
  tf <- terms(fit)
  labels(fit)
}
\keyword{spatial}
\keyword{methods}
