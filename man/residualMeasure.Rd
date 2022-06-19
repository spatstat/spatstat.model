\name{residualMeasure}
\alias{residualMeasure}
\title{
  Residual Measure for an Observed Point Pattern and a Fitted Intensity
}
\description{
  Given a point pattern and an estimate of its intensity function
  obtained in any fashion, compute the residual measure.
}
\usage{
residualMeasure(Q, lambda,
               type = c("raw", "inverse", "Pearson", "pearson"),
               ...)
}
\arguments{
  \item{Q}{
    A point pattern (object of class \code{"ppp"})
    or quadrature scheme (object of class \code{"quad"}).
  }
  \item{lambda}{
    Predicted intensity. An image (object of class \code{"im"})
    or a list of images.
  }
  \item{type}{
    Character string (partially matched) specifying the type of residuals.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{quadscheme}}
    if \code{Q} is a point pattern.
  }
}
\details{
  This command constructs the residual measure
  for the model in which \code{Q} is the observed point pattern
  or quadrature scheme, and \code{lambda} is the estimated intensity
  function obtained in some fashion.
}
\value{
  A measure (object of class \code{"msr"}).
}
\references{
  Baddeley, A., Turner, R., \Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Baddeley, A., \Moller, J. and Pakes, A.G. (2008) 
  Properties of residuals for spatial point processes.
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{60}, 627--649.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{residuals.ppm}}
}
\examples{
  ## nonparametric regression estimate of intensity
  ## as a function of terrain elevation
  f <- rhohat(bei, bei.extra$elev)
  ## predicted intensity as a function of location
  lam <- predict(f)
  ## residuals
  res <- residualMeasure(bei, lam)
  res
  plot(res)
}
\keyword{spatial}
\keyword{models}
