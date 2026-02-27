\name{clusterstrength}
\alias{clusterstrength}
\title{
  Cluster Strength Index
}
\description{
  Calculate the cluster strength index \eqn{\varphi}{\phi}
  for a cluster point process or Cox point process model.
}
\usage{
clusterstrength(object)
}
\arguments{
  \item{object}{
    Clustered point process model.
    Either an object of class \code{"kppm"}
    representing a cluster point process or Cox point process model
    fitted to point pattern data, or
    an object of class \code{"clusterprocess"}
    representing a cluster process model with specified parameters.
  }
}
\details{
  The cluster strength index of a cluster process model 
  is a numerical index which expresses the strength of clustering.
  It is defined as (Baddeley et al., 2022, section 10.2)
  \deqn{
    \varphi = g(0)-1
  }{
    \phi = g(0)-1
  }
  where \eqn{g} is the pair correlation function of the
  cluster process. 

  The index \eqn{\varphi}{\phi} is dimensionless and takes
  non-negative values. Values close to zero indicate that the
  process is close to a Poisson process.

  For a cluster process, \eqn{\varphi}{\phi} is related to the
  sibling probability \eqn{p} by
  \eqn{p = \varphi/(1+\varphi)}{p = \phi/(1+\phi)}.
  
  For a Cox process with driving random intensity
  \eqn{\Lambda(x)}{Lambda(x)},
  \deqn{
    \varphi = \frac{\mbox{var}(\Lambda(0))}{E[\Lambda(0)]^2}
  }{
    \phi = var(Lambda(0))/(E(Lambda(0))^2)
  }
  is a measure of the variability of the random intensity.
}
\value{
  A single numerical value greater than or equal to zero.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{psib}},
  \code{\link{panysib}},
  \code{\link{persist}},
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
  #' Fit a model to the clustered region of full redwood data
  X <- redwoodfull[redwoodfull.extra$regionII]
  fit <- kppm(X)
  clusterstrength(fit)

  #' Create a Thomas model
  m <- clusterprocess("Thomas", kappa=10, mu=5, scale=0.1)
  clusterstrength(m)
}
\keyword{spatial}
\keyword{models}
