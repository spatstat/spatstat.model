\name{methods.zgibbsmodel}
\alias{methods.zgibbsmodel} % DoNotExport
\alias{as.interact.zgibbsmodel}
\alias{as.isf.zgibbsmodel}
\alias{interactionorder.zgibbsmodel}
\alias{intensity.zgibbsmodel}
\alias{is.poisson.zgibbsmodel}
\alias{is.stationary.zgibbsmodel}
\alias{print.zgibbsmodel}
\title{
  Methods for Gibbs Models
}
\description{
  Methods for the experimental class of Gibbs models
}
\usage{
\method{as.interact}{zgibbsmodel}(object)
\method{as.isf}{zgibbsmodel}(object)
\method{interactionorder}{zgibbsmodel}(object)
\method{is.poisson}{zgibbsmodel}(x)
\method{is.stationary}{zgibbsmodel}(x)
\method{print}{zgibbsmodel}(x, \dots)
\method{intensity}{zgibbsmodel}(X, \dots, approx=c("Poisson", "DPP"))
}
\arguments{
  \item{object,x,X}{
    Object of class \code{"zgibbsmodel"}.
  }
  \item{\dots}{Additional arguments.}
  \item{approx}{
    Character string (partially matched) specifying the type
    of approximation.
  }
}
\details{
  Experimental.
}
\value{
  Same as for other methods.
}
\author{
  \adrian
}
\seealso{
  \code{\link{zgibbsmodel}}
}
\examples{
  m <- zgibbsmodel(10, Strauss(0.1), -0.5)
  m
  is.poisson(m)
  is.stationary(m)
  interactionorder(m)
  as.interact(m)
  as.isf(m)
  intensity(m)
  intensity(m, approx="D")
}
\keyword{spatial}
\keyword{models}
