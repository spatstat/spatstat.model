\name{hardcoredist}
\alias{hardcoredist}
\alias{hardcoredist.fii}
\alias{hardcoredist.ppm}
\title{
  Extract the Hard Core Distance of a Point Process Model
}
\description{
  Extract or compute the hard core distance of a point process model.
}
\usage{
hardcoredist(x, \dots)

\method{hardcoredist}{fii}(x, \dots, epsilon = 0)

\method{hardcoredist}{ppm}(x, \dots, epsilon = 0)
}
\arguments{
  \item{x}{
    An object representing a point process model
    (class \code{"ppm"})
    or the interaction structure of a point process
    (class \code{"fii"}) or similar.
  }
  \item{\dots}{
    Additional arguments passed to methods.
  }
  \item{epsilon}{
    Tolerance for defining the hard core.
  }
}
\details{
  A point process model has a hard core distance \code{h}
  if it is impossible for two random points to lie closer than the
  distance \code{h} apart.

  The function \code{hardcoredist} is generic, with methods for
  objects of class \code{"ppm"} (point process models)
  and \code{"fii"} (fitted point process interactions).
  It extracts or computes the hard core distance.

  If \code{epsilon} is specified, then the code calculates the
  largest distance at which the interaction factor is smaller than
  \code{epsilon}, implying that points are unlikely to occur closer
  than this distance.

  The result is zero if the model does not have a hard core distance.
}
\value{
  A single numeric value, or for multitype point processes, a numeric matrix giving the
  hard core distances for each pair of types of points.
}
\author{
  \spatstatAuthors.
}
\examples{
  m <- ppm(cells~1, Hardcore())
  hardcoredist(m)
}
\keyword{spatial}
\keyword{models}
