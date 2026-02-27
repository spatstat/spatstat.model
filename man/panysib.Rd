\name{panysib}
\alias{panysib}
\title{
  Probability that a Point Has Any Siblings
}
\description{
  Given a cluster process model, calculate the probability
  that a point of the process has any siblings.
}
\usage{
panysib(object)
}
\arguments{
  \item{object}{
    Fitted cluster process model (object of class \code{"kppm"})
    or theoretical cluster process model
    (object of class \code{"clusterprocess"}).
  }
}
\details{
  In a Poisson cluster process, two points are called \emph{siblings}
  if they belong to the same cluster, that is, if they had the same
  parent point. This function computes the probability that a 
  given random point has \code{any} siblings.

  If \code{object} is a stationary point process, the result is a single
  number, which is the probability that a typical point of the process
  has any siblings.
  If this number is small,
  then the process is approximately a homogeneous Poisson process
  (complete spatial randomness). The converse is not true (Baddeley et al, 2022).

  Otherwise, the result is a pixel image, in which the value at any
  location \code{u} is the conditional probability,
  given there is a point of the
  process at \code{u}, that this point has any siblings.
  If the pixel values are all small,
  then the process is approximately an inhomogeneous Poisson process.

  This concept was proposed by Baddeley et al (2022).
}
\value{
  A single number (if \code{object} is a stationary point process)
  or a pixel image (otherwise).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{psib}}, \code{\link{persist}}, \code{\link{repul}}.
}
\references{
  Baddeley, A., Davies, T.M., Hazelton, M.L., Rakshit, S. and Turner, R. 
  (2022)
  Fundamental problems in fitting spatial cluster process models.
  \emph{Spatial Statistics} \bold{52}, 100709.
  DOI: \code{10.1016/j.spasta.2022.100709}
}
\examples{
  #' clustered region of full redwood data
  Xc <- redwoodfull[redwoodfull.extra$regionII]
  fitc <- kppm(Xc)
  panysib(fitc)

  if(interactive()) {
     psib(fitc)
     ##' conclusion: strong clustering in region II

     Xr <- redwoodfull[redwoodfull.extra$regionI]
     fitr <- kppm(Xr)
     panysib(fitr)
     psib(fitr)
     ##' panysib is high, psib is low
     unlist(parameters(fitr))
     ##' 
     ##' conclusion: not clustered in region I
  } 

  fitxy <- kppm(redwoodfull ~ polynom(x,y,2))
  plot(panysib(fitxy))
}
\keyword{spatial}
\keyword{models}
