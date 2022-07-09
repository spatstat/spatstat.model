\name{rhohat.ppm} 
\alias{rhohat.ppm}
\alias{rhohat.slrm}
\concept{Resource Selection Function}
\concept{Prospectivity}
\title{
  Nonparametric Estimate of Intensity as Function of a Covariate
}
\description{
  Computes a nonparametric estimate of the intensity of a point process,
  as a function of a (continuous) spatial covariate.
}
\usage{

\method{rhohat}{ppm}(object, covariate, ...,
       weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local", "decreasing", "increasing",
                  "mountain", "valley", "piecewise"),
       subset=NULL, 
       dimyx=NULL, eps=NULL,
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI, breaks=NULL)

\method{rhohat}{slrm}(object, covariate, ...,
       weights=NULL,
       method=c("ratio", "reweight", "transform"),
       horvitz=FALSE,
       smoother=c("kernel", "local", "decreasing", "increasing",
                  "mountain", "valley", "piecewise"),
       subset=NULL, 
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname, confidence=0.95, positiveCI, breaks=NULL)

}
\arguments{
  \item{object}{
    A point pattern (object of class \code{"ppp"} or \code{"lpp"}),
    a quadrature scheme (object of class \code{"quad"})
    or a fitted point process model (object of class \code{"ppm"},
    \code{"slrm"} or \code{"lppm"}).
  }
  \item{covariate}{
    Either a \code{function(x,y)} or a pixel image (object of
    class \code{"im"}) providing the values of the covariate at any
    location.
    Alternatively one of the strings \code{"x"} or \code{"y"}
    signifying the Cartesian coordinates.
  }
  \item{weights}{
    Optional weights attached to the data points.
    Either a numeric vector of weights for each data point,
    or a pixel image (object of class \code{"im"}) or
    a \code{function(x,y)} providing the weights.
  }
  \item{method}{
    Character string determining the estimation method. See Details.
  }
  \item{horvitz}{
    Logical value indicating whether to use Horvitz-Thompson weights.
    See Details.
  }
  \item{smoother}{
    Character string determining the smoothing algorithm
    and the type of curve that will be estimated. See Details.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, from which the estimate should
    be calculated.
  }
  \item{dimyx,eps}{
    Arguments controlling the pixel
    resolution at which the covariate will be evaluated.
    See Details.
  }
  \item{bw}{
    Smoothing bandwidth or bandwidth rule
    (passed to \code{\link{density.default}}).
  }
  \item{adjust}{
    Smoothing bandwidth adjustment factor
    (passed to \code{\link{density.default}}).
  }
  \item{n, from, to}{
    Arguments passed to \code{\link{density.default}} to
    control the number and range of values at which the function
    will be estimated.
  }
  \item{bwref}{
    Optional. An alternative value of \code{bw} to use when smoothing
    the reference density (the density of the covariate values
    observed at all locations in the window).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.default}}
    or \code{\link[locfit]{locfit}}.
  }
  \item{covname}{
    Optional. Character string to use as the name of the covariate.
  }
  \item{confidence}{
    Confidence level for confidence intervals.
    A number between 0 and 1.
  }
  \item{positiveCI}{
    Logical value.
    If \code{TRUE}, confidence limits are always positive numbers;
    if \code{FALSE}, the lower limit of the
    confidence interval may sometimes be negative.
    Default is \code{FALSE} if \code{smoother="kernel"}
    and \code{TRUE} if \code{smoother="local"}.
    See Details.
  }
  \item{breaks}{
    Breakpoints for the piecewise-constant function
    computed when \code{smoother='piecewise'}.
    Either a vector of numeric values specifying the breakpoints,
    or a single integer specifying the number of equally-spaced
    breakpoints. There is a sensible default.
  }
}
\details{
  This command estimates the relationship between
  point process intensity and a given spatial covariate.
  Such a relationship is sometimes called a
  \emph{resource selection function} (if the points are organisms
  and the covariate is a descriptor of habitat) or
  a \emph{prospectivity index} (if the points are mineral deposits
  and the covariate is a geological variable). 
  This command uses nonparametric methods which do not assume a
  particular form for the relationship.  

  If \code{object} is a point pattern, and \code{baseline} is missing or
  null, this command assumes that \code{object} is a realisation of a
  point process with intensity function
  \eqn{\lambda(u)}{lambda(u)} of the form
  \deqn{\lambda(u) = \rho(Z(u))}{lambda(u) = rho(Z(u))}
  where \eqn{Z} is the spatial
  covariate function given by \code{covariate}, and
  \eqn{\rho(z)}{rho(z)} is the resource selection function
  or prospectivity index.
  A nonparametric estimator of the function \eqn{\rho(z)}{rho(z)} is computed.

  If \code{object} is a point pattern, and \code{baseline} is given,
  then the intensity function is assumed to be
  \deqn{\lambda(u) = \rho(Z(u)) B(u)}{lambda(u) = rho(Z(u)) * B(u)}
  where \eqn{B(u)} is the baseline intensity at location \eqn{u}.
  A nonparametric estimator of the relative intensity  \eqn{\rho(z)}{rho(z)}
  is computed.

  If \code{object} is a fitted point process model, suppose \code{X} is
  the original data point pattern to which the model was fitted. Then
  this command assumes \code{X} is a realisation of a Poisson point
  process with intensity function of the form
  \deqn{
    \lambda(u) = \rho(Z(u)) \kappa(u)
  }{
    lambda(u) = rho(Z(u)) * kappa(u)
  }
  where \eqn{\kappa(u)}{kappa(u)} is the intensity of the fitted model
  \code{object}. A nonparametric estimator of
  the relative intensity \eqn{\rho(z)}{rho(z)} is computed.

  The nonparametric estimation procedure is controlled by the
  arguments \code{smoother}, \code{method} and \code{horvitz}.

  The argument \code{smoother} selects the type of estimation technique.
  \itemize{
    \item
    If \code{smoother="kernel"} (the default),
    the nonparametric estimator is a \emph{kernel smoothing estimator}
    of \eqn{\rho(z)}{rho(z)} (Guan, 2008; Baddeley et al, 2012).
    The estimated function \eqn{\rho(z)}{rho(z)} will be
    a smooth function of \eqn{z} which takes nonnegative values.
    Confidence bands are also computed, assuming a Poisson point process.
    See the section on \emph{Smooth estimates}.
    \item
    If \code{smoother="local"}, 
    the nonparametric estimator is a \emph{local regression estimator}
    of \eqn{\rho(z)}{rho(z)} (Baddeley et al, 2012) obtained using
    local likelihood.
    The estimated function \eqn{\rho(z)}{rho(z)} will be
    a smooth function of \eqn{z}.
    Confidence bands are also computed, assuming a Poisson point process.
    See the section on \emph{Smooth estimates}.
    \item
    If \code{smoother="increasing"}, we assume that
    \eqn{\rho(z)}{rho(z)} is an increasing function of \eqn{z},
    and use the \emph{nonparametric maximum likelihood estimator}
    of \eqn{\rho(z)}{rho(z)} described by Sager (1982).
    The estimated function will be a step function, that is increasing
    as a function of \eqn{z}. Confidence bands are not computed.
    See the section on \emph{Monotone estimates}.
    \item
    If \code{smoother="decreasing"}, we assume that
    \eqn{\rho(z)}{rho(z)} is a decreasing function of \eqn{z},
    and use the \emph{nonparametric maximum likelihood estimator}
    of \eqn{\rho(z)}{rho(z)} described by Sager (1982).
    The estimated function will be a step function, that is decreasing
    as a function of \eqn{z}. Confidence bands are not computed.
    See the section on \emph{Monotone estimates}.
    \item
    If \code{smoother="mountain"}, we assume that
    \eqn{\rho(z)}{rho(z)} is a function with an inverted U shape,
    with a single peak at a value \eqn{z_0}{z0}, so that
    \eqn{\rho(z)}{rho(z)} is an increasing function of \eqn{z} 
    for \eqn{z < z_0}{z < z0} and a decreasing function of \eqn{z}
    for \eqn{z > z_0}{z > z0}.
    We compute the \emph{nonparametric maximum likelihood estimator}.
    The estimated function will be a step function, which is 
    increasing and then decreasing as a function of \eqn{z}.
    Confidence bands are not computed.
    See the section on \emph{Unimodal estimates}.
    \item
    If \code{smoother="valley"}, we assume that
    \eqn{\rho(z)}{rho(z)} is a function with a U shape,
    with a single minimum at a value \eqn{z_0}{z0}, so that
    \eqn{\rho(z)}{rho(z)} is a decreasing function of \eqn{z} 
    for \eqn{z < z_0}{z < z0} and an increasing function of \eqn{z}
    for \eqn{z > z_0}{z > z0}.
    We compute the \emph{nonparametric maximum likelihood estimator}.
    The estimated function will be a step function, which is 
    decreasing and then increasing as a function of \eqn{z}.
    Confidence bands are not computed.
    See the section on \emph{Unimodal estimates}.
    \item
    If \code{smoother="piecewise"}, the estimate of
    \eqn{\rho(z)}{rho(z)} is piecewise constant.
    The range of covariate values is divided into several intervals
    (ranges or bands). The endpoints of these intervals are the
    breakpoints, which may be specified by the argument \code{breaks};
    there is a sensible default. The estimate of 
    \eqn{\rho(z)}{rho(z)} takes a constant value on each interval.
    The estimate of \eqn{\rho(z)}{rho(z)} in each interval of covariate
    values is simply the average intensity
    (number of points per unit area) in the relevant sub-region.
  }
  See Baddeley (2018) for a comparison of these estimation techniques
  (except for \code{"mountain"} and \code{"valley"}).
  
  If the argument \code{weights} is present, then the contribution
  from each data point \code{X[i]} to the estimate of \eqn{\rho}{rho} is
  multiplied by \code{weights[i]}.

  If the argument \code{subset} is present, then the calculations are
  performed using only the data inside this spatial region.

  This technique assumes that \code{covariate} has continuous values.
  It is not applicable to covariates with categorical (factor) values
  or discrete values such as small integers.
  For a categorical covariate, use
  \code{\link{intensity.quadratcount}} applied to the result of
  \code{\link{quadratcount}(X, tess=covariate)}.

  The argument \code{covariate} should be a pixel image, or a function,
  or one of the strings \code{"x"} or \code{"y"} signifying the
  cartesian coordinates. It will be evaluated on a fine grid of locations,
  with spatial resolution controlled by the arguments
  \code{dimyx,eps} which are passed to \code{\link{as.mask}}.
}
\section{Smooth estimates}{
  Smooth estimators of \eqn{\rho(z)}{rho(z)}
  were proposed by Baddeley and Turner (2005) and Baddeley et al (2012).
  Similar estimators were proposed by Guan (2008) and in the literature
  on relative distributions (Handcock and Morris, 1999).

  The estimated function \eqn{\rho(z)}{rho(z)} will be a smooth function
  of \eqn{z}.
  
  The smooth estimation procedure involves computing several density estimates
  and combining them. The algorithm used to compute density estimates is 
  determined by \code{smoother}:
  \itemize{
    \item
    If \code{smoother="kernel"},
    the smoothing procedure is based on
    fixed-bandwidth kernel density estimation,
    performed by \code{\link{density.default}}.
    \item
    If \code{smoother="local"}, the smoothing procedure
    is based on local likelihood density estimation, performed by
    \code{\link[locfit]{locfit}}.
  }
  The argument \code{method} determines how the density estimates will be
  combined to obtain an estimate of \eqn{\rho(z)}{rho(z)}:
  \itemize{
    \item
    If \code{method="ratio"}, then \eqn{\rho(z)}{rho(z)} is
    estimated by the ratio of two density estimates,
    The numerator is a (rescaled) density estimate obtained by
    smoothing the values \eqn{Z(y_i)}{Z(y[i])} of the covariate
    \eqn{Z} observed at the data points \eqn{y_i}{y[i]}. The denominator
    is a density estimate of the reference distribution of \eqn{Z}.
    See Baddeley et al (2012), equation (8). This is similar but not
    identical to an estimator proposed by Guan (2008).
    \item
    If \code{method="reweight"}, then \eqn{\rho(z)}{rho(z)} is
    estimated by applying density estimation to the 
    values \eqn{Z(y_i)}{Z(y[i])} of the covariate
    \eqn{Z} observed at the data points \eqn{y_i}{y[i]},
    with weights inversely proportional to the reference density of
    \eqn{Z}.
    See Baddeley et al (2012), equation (9).
    \item 
    If \code{method="transform"},
    the smoothing method is variable-bandwidth kernel
    smoothing, implemented by applying the Probability Integral Transform
    to the covariate values, yielding values in the range 0 to 1,
    then applying edge-corrected density estimation on the interval
    \eqn{[0,1]}, and back-transforming.
    See Baddeley et al (2012), equation (10).
  }
  If \code{horvitz=TRUE}, then the calculations described above
  are modified by using Horvitz-Thompson weighting.
  The contribution to the numerator from 
  each data point is weighted by the reciprocal of the
  baseline value or fitted intensity value at that data point;
  and a corresponding adjustment is made to the denominator.
  
  Pointwise confidence intervals for the true value of \eqn{\rho(z)}
  are also calculated for each \eqn{z},
  and will be plotted as grey shading.
  The confidence intervals are derived using the central limit theorem,
  based on variance calculations which assume a Poisson point process. 
  If \code{positiveCI=FALSE}, the lower limit of the confidence
  interval may sometimes be negative, because the confidence intervals
  are based on a normal approximation to the estimate of \eqn{\rho(z)}.
  If \code{positiveCI=TRUE}, the confidence limits are always
  positive, because the confidence interval is based on a normal
  approximation to the estimate of \eqn{\log(\rho(z))}{log(\rho(z))}.
  For consistency with earlier versions, the default is
  \code{positiveCI=FALSE} for \code{smoother="kernel"}
  and \code{positiveCI=TRUE} for \code{smoother="local"}.
}
\section{Monotone estimates}{
  The nonparametric maximum likelihood estimator
  of a monotone function \eqn{\rho(z)}{rho(z)} was described by Sager (1982).
  This method assumes that
  \eqn{\rho(z)}{rho(z)} is either an increasing
  function of \eqn{z}, or a decreasing function of \eqn{z}.
  The estimated function will be a step function,
  increasing or decreasing as a function of \eqn{z}.

  This estimator is chosen by specifying
  \code{smoother="increasing"} or \code{smoother="decreasing"}.
  The argument \code{method} is ignored this case.
  
  To compute the estimate of \eqn{\rho(z)}{rho(z)}, the algorithm first
  computes several primitive step-function estimates, and then takes
  the maximum of these primitive functions.
  
  If \code{smoother="decreasing"}, each primitive step function
  takes the form \eqn{\rho(z) = \lambda}{rho(z) = lambda} when \eqn{z \le t},
  and \eqn{\rho(z) = 0}{rho(z) = 0} when \eqn{z > t}, where
  and \eqn{\lambda}{lambda} is a primitive estimate of intensity
  based on the data for \eqn{Z \le t}{Z <= t}. The jump location \eqn{t}
  will be the value of the covariate \eqn{Z} at one of the
  data points. The primitive estimate \eqn{\lambda}{lambda}
  is the average intensity (number of points divided by area)
  for the region of space where the covariate value is less than
  or equal to \eqn{t}.

  If \code{horvitz=TRUE}, then the calculations described above
  are modified by using Horvitz-Thompson weighting.
  The contribution to the numerator from 
  each data point is weighted by the reciprocal of the
  baseline value or fitted intensity value at that data point;
  and a corresponding adjustment is made to the denominator.

  Confidence intervals are not available
  for the monotone estimators.
}
\section{Unimodal estimators}{
  If \code{smoother="valley"} then we estimate a U-shaped function.
  A function \eqn{\rho(z)}{rho(z)} is U-shaped if it is
  decreasing when \eqn{z < z_0}{z < z0} and
  increasing when \eqn{z > z_0}{z > z0}, where \eqn{z_0}{z0} is
  called the critical value. The nonparametric maximum likelihood
  estimate of such a function can be computed by profiling over \eqn{z_0}{z0}. 
  The algorithm considers all possible candidate values of the critical value 
  \eqn{z_0}{z0}, and estimates the function \eqn{\rho(z)}{rho(z)}
  separately on the left and right of \eqn{z_0}{z0} using the monotone
  estimators described above. These function estimates are combined into
  a single function, and the Poisson point process likelihood is
  computed. The optimal value of \eqn{z_0}{z0}
  is the one which maximises the Poisson point process likelihood.

  If \code{smoother="mountain"}  then we estimate a function which has
  an inverted U shape. A function \eqn{\rho(z)}{rho(z)} is
  inverted-U-shaped if it is
  increasing when \eqn{z < z_0}{z < z0} and
  decreasing when \eqn{z > z_0}{z > z0}. The nonparametric maximum likelihood
  estimate of such a function can be computed by profiling over
  \eqn{z_0}{z0} using the same technique \emph{mutatis mutandis}.

  Confidence intervals are not available for the unimodal estimators.
}
\value{
  A function value table (object of class \code{"fv"})
  containing the estimated values of \eqn{\rho}{rho}
  (and confidence limits) for a sequence of values of \eqn{Z}.
  Also belongs to the class \code{"rhohat"}
  which has special methods for \code{print}, \code{plot}
  and \code{predict}.
}
\references{
  Baddeley, A., Chang, Y.-M., Song, Y. and Turner, R. (2012)
  Nonparametric estimation of the dependence of a point
  process on spatial covariates.
  \emph{Statistics and Its Interface} \bold{5} (2), 221--236.
  
  Baddeley, A. and Turner, R. (2005)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.  

  Baddeley, A. (2018)
  A statistical commentary on mineral prospectivity analysis.
  Chapter 2, pages 25--65
  in \emph{Handbook of Mathematical Geosciences: Fifty Years of IAMG},
  edited by B.S. Daya Sagar, Q. Cheng and F.P. Agterberg.
  Springer, Berlin.

  Guan, Y. (2008) On consistent nonparametric intensity estimation
  for inhomogeneous spatial point processes.
  \emph{Journal of the American Statistical Association}
  \bold{103}, 1238--1247.

  Handcock, M.S. and Morris, M. (1999)
  \emph{Relative Distribution Methods in the Social Sciences}.
  Springer, New York.
  
  Sager, T.W. (1982) 
  Nonparametric maximum likelihood estimation of
  spatial patterns. \emph{Annals of Statistics} \bold{10}, 1125--1136.
}
\author{
  Smoothing algorithm by
  \adrian, 
  Ya-Mei Chang, Yong Song, 
  and \rolf.
  
  Nonparametric maximum likelihood algorithm by \adrian.
}
\seealso{
  \code{\link{rho2hat}},
  \code{\link{methods.rhohat}},
  \code{\link[spatstat.model]{parres}}.

  See \code{\link[spatstat.model]{ppm}} for a parametric method for the same problem.
}
\examples{
  X <-  rpoispp(function(x,y){exp(3+3*x)})

  fit <- ppm(X ~x)
  rr <- rhohat(fit, "y")

}
\keyword{spatial}
\keyword{models}
\keyword{nonparametric}
