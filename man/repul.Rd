\name{repul}
\alias{repul}
\alias{repul.dppm}
\alias{repul.detpointprocfamily}
\alias{repul.clusterprocess}
\alias{repul.kppm}
\alias{repul.ppm}
\alias{repul.slrm}
\title{
  Repulsiveness Index of a Point Process Model
}
\description{
  Computes a measure of the degree of repulsion between points
  in a point process model.
}
\usage{
repul(model, \dots)

\method{repul}{dppm}(model, \dots)

\method{repul}{detpointprocfamily}(model, \dots)

\method{repul}{clusterprocess}(model, \dots)

\method{repul}{kppm}(model, \dots)

\method{repul}{ppm}(model, \dots)

\method{repul}{slrm}(model, \dots)
}
\arguments{
  \item{model}{
    A fitted point process model
    (object of class \code{"dppm"},
    \code{"kppm"}, \code{"ppm"} or \code{"slrm"})
    or a point process specified by its
    parameters (object of class \code{"detpointprocfamily"}
    or \code{"clusterprocess"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The repulsiveness index \eqn{\mu}{mu} of a
  point process model was defined by Lavancier, \Moller and Rubak (2015) as
  \deqn{
    \mu = \lambda \int (1- g(x)) \, dx
  }{
    mu = lambda * integral[(1-g(x))]
  }
  where \eqn{\lambda}{\lambda} is the intensity of the model and
  \eqn{g(x)} is the pair correlation function, and
  the integral is taken over all two-dimensional vectors \eqn{x}.

  Values of \eqn{\mu} are dimensionless.
  Larger positive values of \eqn{\mu}{mu} indicate stronger repulsion
  between points. 

  The repulsiveness index was originally defined for
  determinantal point processes, but the same definition can be applied
  to other kinds of point process models. The repulsiveness index is
  positive for a determinantal point process model, zero for a Poisson
  process or spatial logistic regression model, negative for
  a cluster process or Cox process, and typically positive for a
  Gibbs process.
  
  If the model is stationary, the result is a single number.

  If the model is not stationary, 
  the result is a pixel image (obtained by multiplying
  the spatially-varying intensity by the integral defined above).

  For Gibbs models of class \code{"ppm"} the calculation uses the
  Poisson-saddlepoint approximation to the pair correlation function,
  and is only implemented for stationary unmarked processes.
}
\value{
  A numeric value or a pixel image.
}
\references{
  Lavancier, F., \Moller, J.  and Rubak, E. (2015),
  Determinantal point process models and statistical inference.
  \emph{Journal of the Royal Statistical Society:
    Series B (Statistical Methodology)},
  \bold{77}, 853--877.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{dppm}},
  \code{\link{psib}}.
}
\examples{
  jpines <- residualspaper$Fig1
  \testonly{
     # smaller dataset for testing
    jpines <- jpines[c(TRUE,FALSE)]
  }
  fit <- dppm(jpines ~ 1, dppGauss)
  repul(fit)

  m <- dppGauss(lambda=10, alpha=0.1, d=2)
  repul(m)
}
\keyword{spatial}
\keyword{models}
