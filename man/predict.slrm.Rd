\name{predict.slrm}
\Rdversion{1.1}
\alias{predict.slrm}
\title{
  Predicted or Fitted Values from Spatial Logistic Regression
}
\description{
  Given a fitted Spatial Logistic Regression model,
  this function computes the fitted probabilities for each pixel,
  or the fitted point process intensity, or the values of the
  linear predictor in each pixel.
}
\usage{
 \method{predict}{slrm}(object, ..., type = "intensity",
                        newdata=NULL, window=NULL,
                        leaveoneout=FALSE, fast=TRUE)
}
\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[spatstat.geom]{pixellate}}
    determining the pixel resolution for the discretisation
    of the point pattern.
  }
  \item{type}{
    Character string (partially) matching one of 
    \code{"probabilities"}, \code{"intensity"} or \code{"link"}.
  }
  \item{newdata}{
    Optional.  List containing new covariate values for the prediction.
    See Details.
  }
  \item{window}{
    Optional.  New window in which to predict.
    An object of class \code{"owin"}.
  }
  \item{leaveoneout}{
    Logical value specifying whether to use leave-one-out prediction.
    Ignored if \code{newdata}, \code{window} or \code{\dots} are given.
  }
  \item{fast}{
    Logical value indicating whether leave-one-out prediction
    should be performed using the leverage (Taylor) approximation
    (\code{fast=TRUE}, the default) or using brute force computation
    (\code{fast=FALSE}). The leverage approximation is very much faster.
  }
}
\details{
  This is a method for \code{\link[stats]{predict}} for spatial logistic
  regression models (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  The argument \code{type} determines which quantity is computed.
  If \code{type="intensity"}), the value of the point process intensity
  is computed at each pixel. If \code{type="probabilities"}) the 
  probability of the presence of a random point in each pixel is
  computed. If \code{type="link"}, the value of the linear predictor is
  computed at each pixel. 

  If \code{newdata = NULL} (the default), the algorithm computes
  fitted values of the model (based on the data that was originally used
  to fit the model \code{object}). 

  If \code{leaveoneout=TRUE}, the fitted values are computed using the
  leave-one-out rule. To compute the fitted value at each pixel
  \code{j},
  \itemize{
    \item the model \code{object} is updated by re-fitting it to
    the dataset \emph{excluding} the data inside pixel \code{j};
    \item the predicted value at pixel \code{j} is computed from
    this updated model.
  }

  If \code{newdata} is given, the algorithm computes predicted values
  of the model, using the new values of the covariates provided by
  \code{newdata}. The argument \code{newdata} should be a list;
  names of entries in the list should correspond
  to variables appearing in the model formula of the \code{object}.
  Each list entry may be a pixel image or a single numeric
  value.

  The argument \code{leaveoneout} is ignored if \code{newdata}
  or \code{window} are given, or if arguments \code{\dots} specify
  a change of pixel resolution.
}
\value{
  A pixel image (object of class \code{"im"}) containing the predicted
  values for each pixel.
}
\seealso{
  \code{\link{slrm}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  plot(predict(fit))

  X <- copper$SouthPoints
  Y <- copper$SouthLines
  Z <- distmap(Y)
  fitc <- slrm(X ~ Z, eps=1)
  pc <- predict(fitc)
  pcminus <- predict(fitc, leaveoneout=TRUE)

  Znew <- distmap(copper$Lines)[copper$SouthWindow]
  pcnew <- predict(fitc, newdata=list(Z=Znew))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
