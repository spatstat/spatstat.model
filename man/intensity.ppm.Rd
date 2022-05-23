\name{intensity.ppm}  
\alias{intensity.ppm}
\title{
  Intensity of Fitted Point Process Model
}
\description{
  Computes the intensity of a fitted point process model.
}
\usage{
 \method{intensity}{ppm}(X, \dots, approx=c("Poisson", "DPP"))
}
\arguments{
  \item{X}{
    A fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.ppm}} in some cases.
    See Details.
  }
  \item{approx}{
    Character string (partially matched) specifying the type of
    approximation to the intensity for a non-Poisson model.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}
  for fitted point process models (class \code{"ppm"}).

  The intensity of a point process model is the expected
  number of random points per unit area. 

  If \code{X} is a Poisson point process model, the intensity of the
  process is computed exactly.
  The result is a numerical value if \code{X}
  is a stationary Poisson point process, and a pixel image if \code{X}
  is non-stationary. (In the latter case, the resolution of the pixel
  image is controlled by the arguments \code{\dots} which are passed
  to \code{\link{predict.ppm}}.)

  If \code{X} is a Gibbs point process model that is not a Poisson model,
  the intensity is computed approximately:
  \itemize{
    \item if \code{approx="Poisson"} (the default),
    the intensity is computed  using the Poisson-saddlepoint approximation
    (Baddeley and Nair, 2012a, 2012b, 2017; Anderssen et al, 2014).
    This approximation is currently available for pairwise-interaction
    models (Baddeley and Nair, 2012a, 2012b)
    and for the area-interaction model and Geyer saturation model
    (Baddeley and Nair, 2017).

    If the model is non-stationary. the pseudostationary solution
    (Baddeley and Nair, 2012b; Anderssen et al, 2014) is used.
    The result is a pixel image,
    whose resolution is controlled by the arguments \code{\dots}
    which are passed to \code{\link{predict.ppm}}.
    
    \item if \code{approx="DPP"}, the intensity is calculated using
    the approximation of (Coeurjolly and Lavancier, 2018) based on a
    determinantal point process. This approximation is more accurate
    than the Poisson saddlepoint approximation, for inhibitory
    interactions. However the DPP approximation is only available
    for stationary pairwise interaction models.
  }
}
\value{
  A numeric value (if the model is stationary)
  or a pixel image.
}
\references{
  Anderssen, R.S., Baddeley, A., DeHoog, F.R. and Nair, G.M. (2014)
  Solution of an integral equation arising in spatial point process theory.
  \emph{Journal of Integral Equations and Applications} 
  \bold{26} (4) 437--453.

  Baddeley, A. and Nair, G. (2012a) 
  Fast approximation of the intensity of Gibbs point processes.
  \emph{Electronic Journal of Statistics} \bold{6} 1155--1169.

  Baddeley, A. and Nair, G. (2012b)
  Approximating the moments of a spatial point process.
  \emph{Stat} \bold{1}, 1, 18--30.
  \code{DOI: 10.1002/sta4.5}

  Baddeley, A. and Nair, G. (2017)
  Poisson-saddlepoint approximation for Gibbs point
  processes with infinite-order interaction: in memory
  of Peter Hall.
  \emph{Journal of Applied Probability} \bold{54}, 4, 1008--1026.

  Coeurjolly, J.-F. and Lavancier, F. (2018)
  Approximation intensity for pairwise interaction Gibbs point processes
  using determinantal point processes.
  \emph{Electronic Journal of Statistics} \bold{12} 3181--3203.
}
\seealso{
  \code{\link{intensity}},
  \code{\link{intensity.ppp}}
}
\examples{
  fitP <- ppm(swedishpines ~ 1)
  intensity(fitP)
  fitS <- ppm(swedishpines ~ 1, Strauss(9))
  intensity(fitS)
  intensity(fitS, approx="D")
  fitSx <- ppm(swedishpines ~ x, Strauss(9))
  lamSx <- intensity(fitSx)
  fitG <- ppm(swedishpines ~ 1, Geyer(9, 1))
  lamG <- intensity(fitG)
  fitA <- ppm(swedishpines ~ 1, AreaInter(7))
  lamA <- intensity(fitA)
}
\author{
  \adrian, Gopalan Nair, and \Frederic Lavancier.
}
\keyword{spatial}
\keyword{models}
