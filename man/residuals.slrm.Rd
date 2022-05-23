\name{residuals.slrm}
\alias{residuals.slrm}
\title{
  Residuals for Fitted Spatial Logistic Regression Model
}
\description{
  Given a spatial logistic regression model fitted to a point pattern,
  compute the residuals for each pixel.
}
\usage{
  \method{residuals}{slrm}(object,
                           type=c("raw", "deviance", "pearson", "working", 
                                  "response", "partial", "score"),
                           \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
    for which residuals should be calculated.
  }
  \item{type}{
    String (partially matched)
    indicating the type of residuals to be calculated.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A pixel image (if the residual values are scalar),
  or a list of pixel images (if the residual values are vectors).
}
\details{
  This function computes several kinds of residuals for the fit of
  a spatial logistic regression model to a spatial point pattern dataset.

  The argument \code{object} must be a fitted spatial logistic
  regression model (object of class \code{"slrm"}). Such objects are
  created by the fitting algorithm \code{\link{slrm}}.

  The residuals are computed for each pixel
  that was used to fit the original model. The residuals are
  returned as a pixel image (if the residual values are scalar),
  or a list of pixel images (if the residual values are vectors).
  
  The type of residual is chosen by the argument \code{type}.

  For a given pixel, suppose \eqn{p} is the fitted probability of
  presence of a point, and \eqn{y} is the presence indicator
  (equal to 1 if the pixel contains any data points, and equal to 0
  otherwise). Then 
  \itemize{
    \item \code{type="raw"} or \code{type="response"} specifies
    the response residual
    \deqn{r = y - p}
    \item \code{type="pearson"} is the Pearson residual
    \deqn{
      r_P = \frac{y - p}{\sqrt{p (1-p)}}
    }{
      rP = (y-p)/sqrt(p * (1-p))
    }
    \item \code{type="deviance"} is the deviance residual
    \deqn{
      r_D = (-1)^{y+1} \sqrt{-2(y log p + (1-y) log(1-p))}
    }{
      rD = (-1)^(y+1) (-2(y log p + (1-y) log(1-p)))^(1/2)
    }
    \item \code{type="score"} specifies the score residuals
    \deqn{
      r_S = (y-p) x
    }{
      rS = (y-p) x
    }
    where \code{x} is the vector of canonical covariate values
    for the pixel
    \item \code{type="working"} specifies the working residuals
    as defined in \code{\link[stats]{residuals.glm}}
    \item \code{type="partial"} specifies the partial residuals
    as defined in \code{\link[stats]{residuals.glm}}
  }
}
\seealso{
  \code{\link[stats]{residuals.glm}},
  \code{\link{residuals.ppm}}
}
\examples{
   d <- if(interactive()) 128 else 32
   H <- unmark(humberside)
   fit <- slrm(H ~ x + y, dimyx=d)

   plot(residuals(fit))

   plot(residuals(fit, type="score"))
}
\author{
  \adrian
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
