\name{clusterfield.kppm}
\alias{clusterfield.kppm}
\title{Field of clusters}
\description{
  Calculate the superposition of cluster kernels at the location of a
  point pattern.
}
\usage{
  \method{clusterfield}{kppm}(model, locations = NULL, \dots)
}
\arguments{
  \item{model}{
    Cluster model. Either a fitted cluster model (object of class
    \code{"kppm"}), a character string specifying the type of cluster
    model, or a function defining the cluster kernel. See Details.
  }
  \item{locations}{
    A point pattern giving the locations of the kernels. Defaults to the
    centroid of the observation window for the \code{"kppm"} method and
    to the center of a unit square otherwise.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[spatstat.explore]{density.ppp}} or the
    cluster kernel. See Details.
  }
}
\details{
  The function \code{clusterfield} is generic,
  with a method for \code{"kppm"} (described here)
  and methods for \code{"character"} and \code{"function"}.
  
  The method \code{clusterfield.kppm} extracts the relevant information
  from the fitted model and calls \code{\link[spatstat.random]{clusterfield.function}}.

  The calculations are performed by \code{\link[spatstat.explore]{density.ppp}} and
  \code{\dots} arguments are passed thereto for control over the pixel
  resolution etc.
  (These arguments are then passed on to \code{\link[spatstat.geom]{pixellate.ppp}}
  and \code{\link[spatstat.geom]{as.mask}}.)

}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{kppm}},

  \code{\link[spatstat.random]{clusterfield}}.
}
\examples{
  fit <- kppm(redwood~1, "Thomas")
  Z <- clusterfield(fit, eps = 0.01)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
