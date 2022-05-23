\name{interactionorder}
\alias{interactionorder}
\alias{interactionorder.isf}
\alias{interactionorder.interact}
\alias{interactionorder.fii}
\alias{interactionorder.ppm}
\title{
  Determine the Order of Interpoint Interaction in a Model
}
\description{
  Given a point process model, report
  the order of interpoint interaction.
}
\usage{
interactionorder(object)

\method{interactionorder}{ppm}(object)

\method{interactionorder}{interact}(object)

\method{interactionorder}{isf}(object)

\method{interactionorder}{fii}(object)
}
\arguments{
  \item{object}{
    A point process model (class \code{"ppm"})
    or similar information.
  }
}
\details{
  This function determines the order of interpoint interaction
  in a Gibbs point process model (or a related object).

  The interaction order is defined as the largest number \code{k} such
  that the probability density of the model contains terms involving \code{k} points at a time.
  For example, in a pairwise interaction
  process such as the Strauss process, the probability density contains
  interaction terms between each pair of points, but does not contain
  any terms that involve three points at a time, so the interaction order is 2.
  
  Poisson point processes have interaction order 1.
  Pairwise-interaction processes have interaction order 2.
  Point processes with the triplet interaction \code{\link{Triplets}}
  have interaction order 3. The Geyer saturation model
  \code{\link{Geyer}} and the area-interaction model
  \code{\link{AreaInter}} have infinite order of interaction.
}
\value{
  A positive integer, or \code{Inf}.
}
\author{
  \spatstatAuthors.
}
\examples{
  interactionorder(ppm(cells ~ 1))
  interactionorder(Strauss(0.1))
  interactionorder(Triplets(0.1))
  interactionorder(Geyer(0.1, 2))
  interactionorder(Hybrid(Strauss(0.1), Triplets(0.2)))
}
\keyword{spatial}
\keyword{models}
