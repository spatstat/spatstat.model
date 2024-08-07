\name{as.owin.ppm}
\alias{as.owin.ppm}
\alias{as.owin.kppm}
\alias{as.owin.dppm}
\alias{as.owin.slrm}
\alias{as.owin.msr}
\title{Convert Data To Class owin}
\description{
Converts data specifying an observation window
in any of several formats, into an object of class \code{"owin"}.
}
\usage{
 \method{as.owin}{ppm}(W, \dots, from=c("points", "covariates"), fatal=TRUE)

 \method{as.owin}{kppm}(W, \dots, from=c("points", "covariates"), fatal=TRUE)

 \method{as.owin}{dppm}(W, \dots, from=c("points", "covariates"), fatal=TRUE)

 \method{as.owin}{slrm}(W, \dots, from=c("points", "covariates"))

 \method{as.owin}{msr}(W, \dots, fatal=TRUE)
}
\arguments{
  \item{W}{
    Data specifying an observation window, in any of several formats
    described under \emph{Details} below.
  }
  \item{fatal}{
    Logical value determining what to do
    if the data cannot be converted to an observation window.
    See Details.
  }
  \item{\dots}{Ignored.}
  \item{from}{Character string. See Details.}
}
\value{
  An object of class \code{"owin"} (see \code{\link[spatstat.geom]{owin.object}})
  specifying an observation window.
}
\details{
  The class \code{"owin"} is a way of specifying the observation window
  for a point pattern. See \code{\link[spatstat.geom]{owin.object}} for an overview.
 
  The generic function \code{\link[spatstat.geom]{as.owin}}
  converts data in any of several formats 
  into an object of class \code{"owin"} for use by the \pkg{spatstat}
  package. The function \code{\link[spatstat.geom]{as.owin}}
  is generic, with methods
  for different classes of objects, and a default method.

  The argument \code{W} may be
  \itemize{
    \item
    an object of class \code{"owin"}
    \item
    a structure with entries \code{xrange}, \code{yrange} specifying the 
    \eqn{x} and \eqn{y} dimensions of a rectangle
    \item
    a structure with entries named \code{xmin}, \code{xmax}, \code{ymin},
    \code{ymax} (in any order) 
    specifying the \eqn{x} and \eqn{y} dimensions of a rectangle.
    This will accept objects of class \code{bbox} in the \code{sf} package.
    \item
    a numeric vector of length 4
    (interpreted as \code{(xmin, xmax, ymin, ymax)} in that order)
    specifying the \eqn{x} and \eqn{y} dimensions of a rectangle
    \item
    a structure with entries named \code{xl}, \code{xu}, \code{yl}, \code{yu}
    (in any order)
    specifying the \eqn{x} and \eqn{y} dimensions of a rectangle
    as \code{(xmin, xmax) = (xl, xu)} and 
    \code{(ymin, ymax) = (yl, yu)}. This will accept objects of
    class \code{spp} used in the Venables and Ripley \pkg{spatial}
    package.
    \item
    an object of class \code{"ppp"} representing a point pattern.
    In this case, the object's \code{window} structure will be
    extracted.
    \item
    an object of class \code{"psp"} representing a line segment pattern.
    In this case, the object's \code{window} structure will be
    extracted.
    \item
    an object of class \code{"tess"} representing a tessellation.
    In this case, the object's \code{window} structure will be
    extracted.
    \item
    an object of class \code{"quad"} representing a quadrature scheme.
    In this case, the window of the \code{data} component will be
    extracted.
    \item
    an object of class \code{"im"} representing a pixel image.
    In this case, a window of type \code{"mask"} will be returned,
    with the same pixel raster coordinates as the image.
    An image pixel value of \code{NA}, signifying that the pixel
    lies outside the window, is transformed into the logical value
    \code{FALSE}, which is the corresponding convention for window masks.
    \item
    an object of class \code{"ppm"}, \code{"kppm"}, \code{"slrm"}
    or \code{"dppm"}
    representing a fitted point process
    model. In this case, if \code{from="data"} (the default),
    \code{as.owin} extracts the  original point
    pattern data to which the model was fitted, and returns the
    observation window of this point pattern. If
    \code{from="covariates"} then \code{as.owin} extracts the
    covariate images to which the model was fitted,
    and returns a binary mask window that specifies the pixel locations.
    \item
    an object of class \code{"lpp"}
    representing a point pattern on a linear network.
    In this case, \code{as.owin} extracts the linear network
    and returns a window containing this network. 
    \item
    an object of class \code{"lppm"}
    representing a fitted point process model on a linear network.
    In this case, \code{as.owin} extracts the linear network
    and returns a window containing this network. 
    \item
    A \code{data.frame} with exactly three columns. Each row of the
    data frame corresponds to one pixel. Each row contains the
    \eqn{x} and \eqn{y} coordinates of a pixel, and a logical value
    indicating whether the pixel lies inside the window.   
    \item
    A \code{data.frame} with exactly two columns. Each row of the
    data frame contains the \eqn{x} and \eqn{y} coordinates of a pixel
    that lies inside the window.
    \item
    an object of class \code{"distfun"}, \code{"nnfun"}
    or \code{"funxy"} representing a function of spatial location,
    defined on a spatial domain. The spatial domain of the function will be
    extracted.
    \item
    an object of class \code{"rmhmodel"} representing a
    point process model that can be simulated using \code{\link[spatstat.random]{rmh}}.
    The window (spatial domain) of the model will be extracted.
    The window may be \code{NULL} in some circumstances (indicating that the
    simulation window has not yet been determined). This is not treated
    as an error, because the argument \code{fatal} defaults to
    \code{FALSE} for this method.
    \item
    an object of class \code{"layered"} representing a
    list of spatial objects. See \code{\link[spatstat.geom]{layered}}.
    In this case, \code{as.owin} will be applied to each
    of the objects in the list, and the union of these windows
    will be returned.
    \item
    an object of some other suitable class from another package.
    For full details, see \code{vignette('shapefiles')}.
  }
  If the argument \code{W} is not in one of these formats
  and cannot be converted to a window, then an error will
  be generated (if \code{fatal=TRUE}) or a value of \code{NULL}
  will be returned (if \code{fatal=FALSE}).

  When \code{W} is a data frame, the argument \code{step}
  can be used to specify the pixel grid spacing; otherwise, the spacing
  will be guessed from the data.
}
\seealso{
  \code{\link[spatstat.geom]{as.owin}},
  \code{\link[spatstat.random]{as.owin.rmhmodel}},
  \code{\link[spatstat.linnet]{as.owin.lpp}}.
  
  \code{\link[spatstat.geom]{owin.object}},
  \code{\link[spatstat.geom]{owin}}.

  Additional methods for \code{as.owin} may be provided
  by other packages outside the \pkg{spatstat} family.
}
\examples{
  fit <- ppm(cells ~ 1)
  as.owin(fit)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
 
 
