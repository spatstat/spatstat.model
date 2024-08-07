\name{vcov.ppm}
\alias{vcov.ppm}
\title{Variance-Covariance Matrix for a Fitted Point Process Model}
\description{
  Returns the variance-covariance matrix of the estimates of the
  parameters of a fitted point process model.
}
\usage{
  \method{vcov}{ppm}(object, \dots,
                    what = c("vcov", "corr", "fisher"),
                    verbose = TRUE,
                    fine=FALSE,
                    gam.action=c("warn", "fatal", "silent"),
                    matrix.action=c("warn", "fatal", "silent"),
                    logi.action=c("warn", "fatal", "silent"),
                    nacoef.action=c("warn", "fatal", "silent"),
                    hessian=FALSE)
}
\arguments{
  \item{object}{A fitted point process model (an object of class \code{"ppm"}.)}
  \item{\dots}{Ignored.}
  \item{what}{Character string (partially-matched)
    that specifies what matrix is returned.
    Options are \code{"vcov"} for the variance-covariance matrix,
    \code{"corr"} for the correlation matrix, and
    \code{"fisher"} or \code{"Fisher"}
    for the Fisher information matrix.
  }
  \item{fine}{
    Logical value indicating whether to use a quick estimate
    (\code{fine=FALSE}, the default) or a slower, more accurate
    estimate (\code{fine=TRUE}).
  }
  \item{verbose}{Logical. If \code{TRUE}, a message will be printed
    if various minor problems are encountered.
  }
  \item{gam.action}{String indicating what to do if \code{object} was
    fitted by \code{gam}. 
  }
  \item{matrix.action}{String indicating what to do if the matrix
    is ill-conditioned (so that its inverse cannot be calculated).
  }
  \item{logi.action}{String indicating what to do if \code{object} was
    fitted via the logistic regression approximation using a
    non-standard dummy point process.
  }
  \item{nacoef.action}{
    String indicating what to do if some of the fitted coefficients
    are \code{NA} (so that variance cannot be calculated).
  }
  \item{hessian}{
    Logical. Use the negative Hessian matrix
    of the log pseudolikelihood instead of the Fisher information.
  }
}
\details{
  This function computes the asymptotic variance-covariance
  matrix of the estimates of the canonical parameters in the
  point process model \code{object}. It is a method for the 
  generic function \code{\link{vcov}}.

  \code{object} should be an object of class \code{"ppm"}, typically
  produced by \code{\link{ppm}}.

  The canonical parameters of the fitted model \code{object}
  are the quantities returned by \code{coef.ppm(object)}.
  The function \code{vcov} calculates the variance-covariance matrix
  for these parameters.
  
  The argument \code{what} provides three options:
  \describe{
    \item{\code{what="vcov"}}{
      return the variance-covariance matrix of the parameter estimates
    }
    \item{\code{what="corr"}}{
      return the correlation matrix of the parameter estimates
    }
    \item{\code{what="fisher"}}{
      return the observed Fisher information matrix.
    }
  }
  In all three cases, the result is a square matrix.
  The rows and columns of the matrix correspond to the canonical
  parameters given by \code{\link{coef.ppm}(object)}. The row and column
  names of the matrix are also identical to the names in
  \code{\link{coef.ppm}(object)}.

  For models fitted by the Berman-Turner approximation (Berman and Turner, 1992;
  Baddeley and Turner, 2000) to the maximum pseudolikelihood (using the
  default \code{method="mpl"} in the call to \code{\link{ppm}}), the implementation works
  as follows.
  \itemize{
    \item
    If the fitted model \code{object} is a Poisson process,
    the calculations are based on standard asymptotic theory for the maximum
    likelihood estimator (Kutoyants, 1998).
    The observed Fisher information matrix of the fitted model
    \code{object} is first computed, by
    summing over the Berman-Turner quadrature points in the fitted model.
    The asymptotic variance-covariance matrix is calculated as the
    inverse of the
    observed Fisher information. The correlation matrix is then obtained
    by normalising.
    \item
    If the fitted model is not a Poisson process (i.e. it is some other
    Gibbs point process) then the calculations are based on
    Coeurjolly and Rubak (2012). A consistent estimator of the
    variance-covariance matrix is computed by summing terms over all
    pairs of data points. If required, the Fisher information is
    calculated as the inverse of the variance-covariance matrix.
  }

  For models fitted by the Huang-Ogata method (\code{method="ho"} in
  the call to \code{\link{ppm}}), the implementation uses the 
  Monte Carlo estimate of the Fisher information matrix that was
  computed when the original model was fitted. 

  For models fitted by the logistic regression approximation to the
  maximum pseudolikelihood (\code{method="logi"} in the call to
  \code{\link{ppm}}),
  \itemize{
    \item Calculations are based on Baddeley \emph{et al.} (2013).
    A consistent estimator of the variance-covariance matrix is
    computed by summing terms over all pairs of data points.
    If required, the Fisher information is calculated as the inverse of the
    variance-covariance matrix.
    \item The calculations depend on
    the type of dummy pattern used when the model was fitted:
    \itemize{
      \item currently only the dummy types
      \code{"stratrand"} (the default), \code{"binomial"} and \code{"poisson"} as
      generated by \code{\link[spatstat.geom]{quadscheme.logi}} are supported.
      \item For other dummy types the behavior depends on the argument
      \code{logi.action}. If \code{logi.action="fatal"} an error is
      produced. Otherwise, for dummy types
      \code{"grid"} and \code{"transgrid"} the formulas for
      \code{"stratrand"} are used which in many cases should be
      conservative. For an arbitrary, user-specified dummy pattern (type
      \code{"given"}), the formulas for \code{"poisson"} are used which in
      many cases should be conservative. If \code{logi.action="warn"} a
      warning is issued, otherwise the calculation proceeds without a
      warning.
    }
    \item The result of the calculation is \bold{random} (i.e. not
    deterministic) when dummy type is \code{"stratrand"} (the default)
    because some of the variance terms are estimated by random sampling.
    This can be avoided by specifying
    \code{dummytype='poisson'} or \code{dummytype='binomial'} in the 
    call to \code{\link{ppm}} when the model is fitted.
  }
  
  The argument \code{verbose} makes it possible to suppress some
  diagnostic messages.

  The asymptotic theory is not correct if the model was fitted using
  \code{gam} (by calling \code{\link{ppm}} with \code{use.gam=TRUE}).
  The argument \code{gam.action} determines what to do in this case.
  If \code{gam.action="fatal"}, an error is generated.
  If \code{gam.action="warn"}, a warning is issued and the calculation
  proceeds using the incorrect theory for the parametric case, which is
  probably a reasonable approximation in many applications.
  If \code{gam.action="silent"}, the calculation proceeds without a
  warning.
  
  If \code{hessian=TRUE} then the negative Hessian (second derivative)
  matrix of the log pseudolikelihood, and its inverse, will be computed.
  For non-Poisson models, this is not a valid estimate of variance,
  but is useful for other calculations.

  Note that standard errors and 95\% confidence intervals for
  the coefficients can also be obtained using
  \code{confint(object)} or \code{coef(summary(object))}.
}
\section{Error messages}{
  An error message that reports
  \emph{system is computationally singular} indicates that the
  determinant of the Fisher information matrix was either too large 
  or too small for reliable numerical calculation.

  If this message occurs, try repeating the calculation
  using \code{fine=TRUE}. 

  Singularity can occur because of numerical overflow or
  collinearity in the covariates. To check this, rescale the 
  coordinates of the data points and refit the model. See the Examples.

  In a Gibbs model, a singular matrix may also occur if the
  fitted model is a hard core process: this is a feature of the
  variance estimator. 
}
\value{
  A square matrix.
}
\examples{
  X <- rpoispp(42)
  fit <- ppm(X ~ x + y)
  vcov(fit)
  vcov(fit, what="Fish")

  # example of singular system
  m <- ppm(demopat ~polynom(x,y,2))
  \donttest{
    try(v <- vcov(m))
  }
  # rescale x, y coordinates to range [0,1] x [0,1] approximately
  demopatScale <- rescale(demopat, 10000)
  m <- ppm(demopatScale ~ polynom(x,y,2))
  v <- vcov(m)

  # Gibbs example
  fitS <- ppm(swedishpines ~1, Strauss(9))
  coef(fitS)
  sqrt(diag(vcov(fitS)))
}
\author{
  Original code for Poisson point process was written by
  \adrian and \rolf.
  New code for stationary Gibbs point processes was generously contributed by
  \ege and Jean-\Francois Coeurjolly.
  New code for generic Gibbs process written by \adrian.
  New code for logistic method written by \ege.
}
\seealso{
  \code{\link{vcov}} for the generic,

  \code{\link{ppm}} for information about fitted models,
  
  \code{\link[stats]{confint}} for confidence intervals.
}
\references{
  Baddeley, A., Coeurjolly, J.-F., Rubak, E. and Waagepetersen, R. (2014)
  Logistic regression for spatial Gibbs point processes.
  \emph{Biometrika} \bold{101} (2) 377--392.

  Coeurjolly, J.-F. and Rubak, E. (2013)
  Fast covariance estimation for innovations
  computed from a spatial Gibbs point process.
  Scandinavian Journal of Statistics \bold{40} 669--684.

  Kutoyants, Y.A. (1998) 
  \bold{Statistical Inference for Spatial Poisson Processes},
  Lecture Notes in Statistics 134. 
  New York: Springer 1998. 
}
\keyword{spatial}
\keyword{methods}
\keyword{models}

