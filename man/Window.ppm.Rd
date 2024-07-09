\name{Window.ppm}
\alias{Window.ppm}
\alias{Window.kppm}
\alias{Window.dppm}
\alias{Window.slrm}
\alias{Window.msr}
\title{Extract Window of Spatial Object}
\description{
  Given a spatial object (such as a point pattern or pixel image)
  in two dimensions, these functions extract the window
  in which the object is defined.
}
\usage{
 \method{Window}{ppm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{kppm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{dppm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{slrm}(X, \dots, from=c("points", "covariates"))

 \method{Window}{msr}(X, \dots)




}
\arguments{
  \item{X}{A spatial object.}
  \item{\dots}{Ignored.}
  \item{from}{Character string. See Details.}
}
\value{
  An object of class \code{"owin"} (see \code{\link[spatstat.geom]{owin.object}})
  specifying an observation window.
}
\details{
  These are methods for the generic function \code{\link[spatstat.geom]{Window}}
  which extract the spatial window in which the object \code{X}
  is defined. 
  The argument \code{from} applies when \code{X} is a fitted 
  two-dimensional point process model
  (object of class \code{"ppm"}, \code{"kppm"}, \code{"slrm"} or \code{"dppm"}).
  If \code{from="data"} (the default),
  \code{Window} extracts the window of the original point
  pattern data to which the model was fitted.
  If \code{from="covariates"} then \code{Window} returns the
  window in which the spatial covariates of the model were provided.
}
\seealso{
  \code{\link[spatstat.geom]{Window}},
  \code{\link[spatstat.geom]{Window.ppp}},
  \code{\link[spatstat.geom]{Window.psp}}.
  
  \code{\link[spatstat.geom]{owin.object}}
}
\examples{
   A <- ppm(cells ~ 1)
   Window(A)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
 