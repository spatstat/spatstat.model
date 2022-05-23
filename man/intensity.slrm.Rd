\name{intensity.slrm}  
\alias{intensity.slrm}
\title{
  Intensity of Fitted Spatial Logistic Regression Model
}
\description{
  Computes the intensity of a fitted spatial logistic regression model,
  treated as a point process model.
}
\usage{
\method{intensity}{slrm}(X, \dots)
}
\arguments{
  \item{X}{
    A fitted spatial logistic regression model (object of class \code{"slrm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.slrm}} in some cases.
    See Details.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}
  for spatial logistic regression models (class \code{"slrm"}).

  The fitted spatial logistic regression model \code{X} is interpreted
  as a point process model. The intensity of a point process model is
  defined as the expected number of random points per unit area. The
  fitted probabilities of presence according to \code{X} are converted
  to intensity values.

  The result is a numerical value if \code{X}
  is stationary, and a pixel image if \code{X}
  is non-stationary. In the latter case, the resolution of the pixel
  image is controlled by the arguments \code{\dots} which are passed
  to \code{\link{predict.slrm}}.
}
\value{
  A numeric value (if the model is stationary)
  or a pixel image.
}
\references{
  Baddeley, A., Berman, M., Fisher, N.I., Hardegen, A., Milne, R.K.,
  Schuhmacher, D., Shah, R. and Turner, R. (2010)
  Spatial logistic regression and change-of-support
  for spatial Poisson point processes.
  \emph{Electronic Journal of Statistics}
  \bold{4}, 1151--1201.
  \code{DOI: 10.1214/10-EJS581}
}
\seealso{
  \code{\link{intensity}},
  \code{\link{intensity.ppm}}
}
\examples{
  fitS <- slrm(swedishpines ~ 1)
  intensity(fitS)
  fitX <- slrm(swedishpines ~ x)
  intensity(fitX)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
