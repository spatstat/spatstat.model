\name{persist}
\alias{persist}
\title{
  Spatial Persistence Index
}
\description{
  Given a cluster process model fitted to point pattern data,
  calculate the spatial persistence index for the model
  over the observation window.
}
\usage{
persist(object, W=Window(object))
}
\arguments{
  \item{object}{
    Fitted cluster process model (object of class \code{"kppm"}).
  }
  \item{W}{
    Spatial region (object of class \code{"owin"}) for which the
    persistence index should be calculated.
    Default is the window in which the original point pattern dataset
    was observed.
  }
      
}
\details{
  The spatial persistence index of a cluster process model
  is a numerical index which expresses the spatial scale of the model
  relative to the size of the window in which the data were observed.
  It is defined as (Baddeley et al., 2022, section 10.2)
  \deqn{
    v = \frac{g(d) - 1}{g(0)-1}
  }{
    v = (g(d) - 1)/(g(0) - 1)
  }
  where \eqn{g} is the pair correlation function of the
  cluster process, and \eqn{d} is the diameter of the observation window
  of the original point pattern dataset to which the model was fitted.

  The index \eqn{v} is dimensionless and takes values between 0 and 1.
  It depends on both the fitted cluster process, and on the
  window in which the original data were observed. 
  (The user can specify a different observation window \code{W},
  for which the persistence index should be calculated.)

  The spatial persistence index effectively measures the
  size of a typical cluster in the cluster process
  (observed within the observation window)
  as a fraction of the size of the observation window.
  Values of \eqn{v} close to 1 indicate that the
  clusters are so large that the model (observed within the observation window)
  is effectively a mixed Poisson process.
}
\value{
  A single numerical value between 0 and 1.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{psib}},
  \code{\link{panysib}},
  \code{\link{clusterstrength}},
  \code{\link{repul}}
}
\references{
  Baddeley, A., Davies, T.M., Hazelton, M.L., Rakshit, S. and Turner, R. 
  (2022)
  Fundamental problems in fitting spatial cluster process models.
  \emph{Spatial Statistics} \bold{52}, 100709.
  DOI: \code{10.1016/j.spasta.2022.100709}
}
\examples{
  #' non-clustered region of full redwood data
  X <- redwoodfull[redwoodfull.extra$regionI]
  fit <- kppm(X)

  persist(fit)

  psib(fit)
  panysib(fit)

  ## conclusion: fitted model is approximately a mixed Poisson process
}
\keyword{spatial}
\keyword{models}
