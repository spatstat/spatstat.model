\name{intensity.mppm}  
\alias{intensity.mppm}
\title{
  Intensity of Fitted Point Process Model for Multiple Patterns
}
\description{
  Computes the intensity of a point process model
  fitted to multiple point patterns.
}
\usage{
 \method{intensity}{mppm}(X, \dots)
}
\arguments{
  \item{X}{
    A point process model fitted to multiple point pattern
    datasets (object of class \code{"mppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{intensity.ppm}}.
  }
}
\details{
  This is a method for the generic function
  \code{\link[spatstat.geom]{intensity}} for objects of class \code{"mppm"}.
  Such an object is created
  by the model-fitting function \code{\link{mppm}} and represents
  a point process model fitted to multiple point patterns.

  The intensity of a point process model is the expected
  number of random points per unit area. 

  For each of the original point pattern datasets to which \code{X} was
  fitted, the intensity of the fitted model will be computed
  by \code{\link[spatstat.model]{intensity.ppm}}. These results will be coerced to
  an appropriate format.

  For stationary point process models,
  the result of \code{intensity.mppm(X)}
  will be a numeric vector giving the fitted intensity
  values for each dataset, or (if the data were multitype point
  patterns) a numeric matrix giving the fitted
  intensity for each type of point for each dataset.
  
  For non-stationary models, the result of \code{intensity.mppm(X)} will
  be a list of pixel images giving the fitted intensity for each dataset,
  or a list of lists of pixel images giving the fitted intensity of each
  type of point for each dataset.

  When additional arguments \code{\dots}
  are given, the result will be a list containing the results of
  \code{intensity.ppm} for each dataset.
}
\value{
  A numeric vector or matrix, a list of pixel images,
  or a list of lists of pixel images, or a list containing data of
  another type.
}
\references{
  \baddrubaturnbook
}
\seealso{
  \code{\link[spatstat.geom]{intensity}},
  \code{\link[spatstat.model]{intensity.ppm}}.
  
  Type \code{methods(intensity)} to see methods for other classes.
}
\examples{
  h <- hyperframe(Bugs=waterstriders)
  m1 <- mppm(Bugs ~ 1, data=h)
  intensity(m1)
  mx <- mppm(Bugs ~ x, data=h)
  intensity(mx)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
