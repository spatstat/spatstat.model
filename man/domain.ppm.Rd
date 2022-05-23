\name{domain.ppm}
\alias{domain.ppm}
\alias{domain.kppm}
\alias{domain.dppm}
\alias{domain.slrm}
\alias{domain.msr}
\title{
  Extract the Domain of any Spatial Object
}
\description{
  Given a spatial object such as a point pattern, in any number of dimensions,
  this function extracts the spatial domain in which the object is defined.
}
\usage{

 \method{domain}{ppm}(X, \dots, from=c("points", "covariates"))

 \method{domain}{kppm}(X, \dots, from=c("points", "covariates"))

 \method{domain}{dppm}(X, \dots, from=c("points", "covariates"))

 \method{domain}{slrm}(X, \dots, from=c("points", "covariates"))

 \method{domain}{msr}(X, \dots)




}
\arguments{
  \item{X}{
    A spatial object such as a point pattern (in any number
    of dimensions), line segment pattern or pixel image.
  }
  \item{\dots}{
    Extra arguments. They are ignored by all the methods listed here.
  }
  \item{from}{Character string. See Details.}
}
\details{
  The function \code{\link[spatstat.geom]{domain}} is generic.

  For a spatial object \code{X} in any number of dimensions, 
  \code{domain(X)} extracts the spatial domain in which \code{X} is
  defined.

  For a two-dimensional object \code{X}, typically \code{domain(X)}
  is the same as \code{Window(X)}.

  Exceptions occur for methods related to linear networks.

  The argument \code{from} applies when \code{X} is a fitted 
  point process model
  (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
  If \code{from="data"} (the default),
  \code{domain} extracts the window of the original point
  pattern data to which the model was fitted.
  If \code{from="covariates"} then \code{domain} returns the
  window in which the spatial covariates of the model were provided.
}
\value{
  A spatial object representing the domain of \code{X}.
  Typically a window (object of class \code{"owin"}),
  a three-dimensional box (\code{"box3"}), a multidimensional
  box (\code{"boxx"}) or a linear network (\code{"linnet"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.geom]{domain}},
  \code{\link[spatstat.geom]{domain.quadratcount}},
  \code{\link[spatstat.explore]{domain.quadrattest}},
  \code{\link[spatstat.random]{domain.rmhmodel}},
  \code{\link[spatstat.linnet]{domain.lpp}}.
  \code{\link{Window}},
  \code{\link{Frame}}.
}
\examples{
  domain(ppm(redwood ~ 1))
}
\keyword{spatial}
\keyword{manip}
