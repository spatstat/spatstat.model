\name{residuals.rppm}
\alias{residuals.rppm}
\title{
  Residuals for Recursively Partitioned Point Process Model
}
\description{
  Given a point process model
  that was fitted to a point pattern
  by recursive partitioning (regression tree) methods,
  compute the residual measure.
}
\usage{
  \method{residuals}{rppm}(object,
                           type=c("raw", "inverse", "Pearson"),
                           \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
    for which residuals should be calculated.
  }
  \item{type}{
    String (partially matched)
    indicating the type of residuals to be calculated.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A measure (object of class \code{"msr"}).
}
\details{
  This function computes the residual measure for a
  point process model that was fitted to point pattern data
  by recursive partitioning of the covariates.

  The argument \code{object} must be a fitted model object
  of class \code{"rppm"}). Such objects are
  created by the fitting algorithm \code{\link{rppm}}.

  The type of residual is chosen by the argument \code{type}.
}
\seealso{
  \code{\link{residuals.ppm}}
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
\examples{
  fit <- rppm(bei ~ elev + grad, data=bei.extra)
  res <- residuals(fit)
  plot(res)
}
\author{
  \adrian.
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
