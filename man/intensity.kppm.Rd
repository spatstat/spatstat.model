\name{intensity.kppm}  
\alias{intensity.kppm}
\title{
  Intensity of Fitted Point Process Model (Cluster Process or Cox Process)
}
\description{
  Computes the intensity of a fitted point process model
  (Neyman-Scott cluster process or Cox process).
}
\usage{
 \method{intensity}{kppm}(X, \dots)
}
\arguments{
  \item{X}{
    A fitted point process model (object of class \code{"kppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.ppm}} in some cases.
    See Details.
  }
}
\details{
  This is a method for the generic function \code{\link[spatstat.geom]{intensity}}
  for fitted point process models of class \code{"kppm"}. An object of
  this class is created by the model-fitting function \code{\link{kppm}} and
  represents a fitted Neyman-Scott cluster process or Cox
  process model (Baddeley, Rubak and Turner, 2015, chapter 12). 

  The intensity of a point process model is the expected
  number of random points per unit area. 

  The result of \code{intensity.kppm(X)} is a numerical value if \code{X}
  is a stationary point process, and a pixel image if \code{X}
  is non-stationary. In the latter case, the resolution of the pixel
  image is controlled by the arguments \code{\dots} which are passed
  to \code{\link{predict.ppm}}.
}
\value{
  A numeric value (if the model is stationary)
  or a pixel image.
}
\references{
  \baddrubaturnbook
}
\seealso{
  \code{\link[spatstat.geom]{intensity}},
  \code{\link[spatstat.model]{intensity.ppm}}
}
\examples{
   m <- kppm(redwood ~ 1, "LGCP")
   intensity(m)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
