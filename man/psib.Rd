\name{psib}
\alias{psib}
\alias{psib.kppm}
\title{
  Sibling Probability of Cluster Point Process
}
\description{
  Computes the sibling probability of a cluster point process model.
}
\usage{
  psib(object)

  \method{psib}{kppm}(object)
}
\arguments{
  \item{object}{
    Fitted cluster point process model
    (object of class \code{"kppm"}).
  }
}
\details{
  In a Poisson cluster process, two points are called \emph{siblings}
  if they belong to the same cluster, that is, if they had the same
  parent point. If two points of the process are
  separated by a distance \eqn{r}, the probability that
  they are siblings is \eqn{p(r) = 1 - 1/g(r)} where \eqn{g} is the
  pair correlation function of the process.
  
  The value \eqn{p(0) = 1 - 1/g(0)} is the probability that,
  if two points of the process are situated very close to each other,
  they came from the same cluster. This probability
  is an index of the strength of clustering, with high values
  suggesting strong clustering.

  This concept was proposed in Baddeley, Rubak and Turner (2015, page 479)
  and Baddeley (2017). It was shown in Baddeley et al (2022) that
  the sibling probability is directly related to the strength of clustering.
}
\value{
  A single number.
}
\references{
  Baddeley, A. (2017)
  Local composite likelihood for spatial point processes.
  \emph{Spatial Statistics} \bold{22}, 261--295.

  \baddrubaturnbook

  Baddeley, A., Davies, T.M., Hazelton, M.L., Rakshit, S. and Turner, R. 
  (2022)
  Fundamental problems in fitting spatial cluster process models.
  \emph{Spatial Statistics} \bold{52}, 100709.
  DOI: \code{10.1016/j.spasta.2022.100709}
}
\author{
  \adrian.
}
\seealso{
  \code{\link{kppm}}, \code{\link{panysib}}, \code{\link{persist}},
  \code{\link{repul}}.
}
\examples{
  fit <- kppm(redwood3 ~1, "Thomas")
  psib(fit)
}
\keyword{spatial}
\keyword{models}
