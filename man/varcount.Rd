\name{varcount}
\alias{varcount}
\title{
  Predicted Variance of the Number of Points 
}
\description{
  Given a fitted point process model, calculate the predicted variance
  of the number of points in a nominated set \code{B}.
}
\usage{
varcount(model, B=Window(model), \dots, dimyx = NULL, relative=FALSE)
}
\arguments{
  \item{model}{
    A fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).
  }
  \item{B}{
    A window (object of class \code{"owin"} specifying the region in
    which the points are counted.
    Alternatively a pixel image (object of class \code{"im"})
    or a function of spatial coordinates specifying a numerical weight
    for each random point.
    The default is the window of the original point pattern data to which
    the model was fitted.
  }
  \item{\dots}{
    Additional arguments passed to \code{B} when it is a function.
  }
  \item{dimyx}{
    Spatial resolution for the calculations.
    Argument passed to \code{\link[spatstat.geom]{as.mask}}.
  }
  \item{relative}{
    Logical value specifying whether to divide the variance
    by the mean value.
  }
}
\details{
  The function \code{varcount} calculates the variance of the number of points
  falling in a specified window \code{B} according to the \code{model}.
  It can also calculate the variance of a sum of weights attached
  to each random point.

  If \code{relative=FALSE} (the default), the result is the variance.
  If \code{relative=TRUE}, the result is the variance divided by the
  mean, which is the overdispersion index (equal to 1 if the number of
  points has a Poisson distribution).
  
  The \code{model} should be a fitted point process model
  (object of class \code{"ppm"}, \code{"kppm"} or \code{"dppm"}).

  \itemize{
    \item{
      If \code{B} is a window, \code{varcount} calculates the variance
      of the number of points falling in \code{B}, according to the
      fitted \code{model}.

      If the \code{model} depends on spatial covariates other than the
      Cartesian coordinates, then \code{B} should be a subset of the
      domain in which these covariates are defined. 
    }
    \item{
      If \code{B} is a pixel image,
      \code{varcount} calculates the variance of
      \eqn{T = \sum_i B(x_i)}{T = sum[i] B(x[i])},
      the sum of the values of \code{B} over all random points
      falling in the domain of the image.

      If the \code{model} depends on spatial covariates other than the
      Cartesian coordinates, then the domain of the pixel image,
      \code{as.owin(B)}, should be a subset of the domain in which these
      covariates are defined. 
    }
    \item{
      If \code{B} is a \code{function(x,y)} or \code{function(x,y,...)}
      then \code{varcount} calculates the variance of
      \eqn{T = \sum_i B(x_i)}{T = sum[i] B(x[i])},
      the sum of the values of \code{B} over all random points
      falling inside the window \code{W=as.owin(model)}, the window
      in which the original data were observed.
    }
  }

  The variance calculation involves the intensity and the
  pair correlation function of the model.
  The calculation is exact (up to discretisation error)
  for models of class \code{"kppm"} and \code{"dppm"},
  and for Poisson point process models of class \code{"ppm"}.
  For Gibbs point process models of class \code{"ppm"} the
  calculation depends on the Poisson-saddlepoint approximations
  to the intensity and pair correlation function, which are rough
  approximations. The approximation is not yet implemented
  for some Gibbs models.
}
\value{
  A single number.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{predict.ppm}},
  \code{\link{predict.kppm}},
  \code{\link{predict.dppm}}
}
\examples{
   fitT <- kppm(redwood ~ 1, "Thomas")
   B <- owin(c(0, 0.5), c(-0.5, 0))
   varcount(fitT, B)

   fitS <- ppm(swedishpines ~ 1, Strauss(9))
   BS <- square(50)
   varcount(fitS, BS)
}
\keyword{spatial}
\keyword{models}

