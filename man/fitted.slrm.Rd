\name{fitted.slrm}
\Rdversion{1.1}
\alias{fitted.slrm}
\title{
  Fitted Probabilities for Spatial Logistic Regression
}
\description{
  Given a fitted Spatial Logistic Regression model,
  this function computes the fitted probabilities for each pixel,
  or the fitted probabilities at each original data point.
}
\usage{
  \method{fitted}{slrm}(object, \dots, type="probabilities",
                 dataonly=FALSE, leaveoneout=FALSE)
}
\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Ignored.
  }
  \item{type}{
    Character string (partially) matching one of 
    \code{"probabilities"}, \code{"intensity"} or \code{"link"}
    determining the quantity that should be predicted.
  }
  \item{dataonly}{
    Logical. If \code{TRUE}, then values will only be computed
    at the points of the data point pattern. If \code{FALSE}, then
    values will be computed at the pixels used to fit the model.
  }
  \item{leaveoneout}{
    Logical value
    specifying whether to perform a leave-one-out calculation
    when \code{dataonly=TRUE}. If \code{leaveoneout=TRUE},
    the fitted value at each data point \code{X[i]} is calculated
    by re-fitting the model to the data with \code{X[i]} removed.
  }
}
\details{
  This is a method for the generic function
  \code{\link[stats:fitted.values]{fitted}}
  for spatial logistic regression models
  (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  By default, the algorithm computes
  the fitted probabilities of the presence of a random point in each
  pixel,
  and returns them as an image.

  If \code{dataonly=TRUE}, the algorithm computes the fitted
  presence probabilities only at the locations of the original data points.
}
\value{
  A pixel image (object of class \code{"im"}) containing the
  fitted probability for each pixel,
  or a numeric vector containing the fitted probability at each data point.
}
\seealso{
  \code{\link{slrm}},
  \code{\link[stats:fitted.values]{fitted}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  plot(fitted(fit))
  fitted(fit, dataonly=TRUE)
}
\author{\adrian
  and \rolf.
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
