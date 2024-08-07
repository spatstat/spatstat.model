\name{mincontrast}
\alias{mincontrast}
\title{Method of Minimum Contrast}
\description{
  A general low-level algorithm for fitting theoretical point process models
  to point pattern data by the Method of Minimum Contrast.
}
\usage{
mincontrast(observed, theoretical, startpar, \dots,
          ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=NULL),
          fvlab=list(label=NULL, desc="minimum contrast fit"),
          explain=list(dataname=NULL, modelname=NULL, fname=NULL),
          action.bad.values=c("warn", "stop", "silent"),
          control=list(), stabilize=TRUE,
          pspace=NULL)
}
\arguments{
  \item{observed}{
    Summary statistic, computed for the data.
    An object of class \code{"fv"}.
  }
  \item{theoretical}{
    An R language function that calculates the theoretical expected value
    of the summary statistic, given the model parameters.
    See Details.
  }
  \item{startpar}{
    Vector of initial values of the parameters of the
    point process model (passed to \code{theoretical}). 
  }
  \item{\dots}{
    Additional arguments passed to the function \code{theoretical}
    and to the optimisation algorithm \code{\link[stats]{optim}}.
  }
  \item{ctrl}{
    Optional. List of arguments controlling the optimisation. See Details.
  }
  \item{fvlab}{
    Optional. List containing some labels for the return value. See Details.
  }
  \item{explain}{
    Optional. List containing strings that give a human-readable description
    of the model, the data and the summary statistic.
  }
  \item{action.bad.values}{
    String (partially matched) specifying what to do if
    values of the summary statistic are \code{NA}, \code{NaN} or
    infinite. See Details.
  }
  \item{control}{
    Optional. Argument passed to \code{\link[stats]{optim}}.
    A list of parameters which control the behaviour of the optimization algorithm.
  }
  \item{stabilize}{
    Logical value specifying whether to numerically stabilize
    the optimization algorithm, by specifying suitable default values of
    \code{control$fnscale} and \code{control$parscale}.
  }
  \item{pspace}{
    For internal use by the package only.
  }
}
\details{
  This function is a general algorithm for fitting point process models
  by the Method of Minimum Contrast. If you want to fit the
  Thomas process, see \code{\link{thomas.estK}}.
  If you want to fit a log-Gaussian Cox process, see
  \code{\link{lgcp.estK}}.  If you want to fit the \Matern cluster
  process, see \code{\link{matclust.estK}}.

  The Method of Minimum Contrast (Pfanzagl, 1969; Diggle and Gratton, 1984)
  is a general technique for fitting
  a point process model to point pattern data. First a summary function
  (typically the \eqn{K} function) is computed from the data point
  pattern. Second, the theoretical expected
  value of this summary statistic under the point process model
  is derived (if possible, as an algebraic expression involving the
  parameters of the model) or estimated from simulations of the model. 
  Then the model is fitted by finding the optimal parameter values
  for the model to give the closest match between the theoretical
  and empirical curves. 

  The argument \code{observed} should be an object of class \code{"fv"}
  (see \code{\link[spatstat.explore]{fv.object}}) containing the values of a summary
  statistic computed from the data point pattern. Usually this is the
  function \eqn{K(r)} computed by \code{\link[spatstat.explore]{Kest}} or one of its relatives.
  
  The argument \code{theoretical} should be a user-supplied function
  that computes the theoretical expected value of the summary statistic.
  It must have an argument named \code{par} that will be the vector
  of parameter values for the model (the length and format of this
  vector are determined by the starting values in \code{startpar}).
  The function \code{theoretical} should also expect a second argument
  (the first argument other than \code{par})
  containing values of the distance \eqn{r} for which the theoretical
  value of the summary statistic \eqn{K(r)} should be computed.
  The value returned by \code{theoretical} should be a vector of the
  same length as the given vector of \eqn{r} values.

  The argument \code{ctrl} determines the contrast criterion
  (the objective function that will be minimised). 
  The algorithm minimises the criterion
  \deqn{
    D(\theta)=
    \int_{r_{\mbox{\scriptsize min}}}^{r_{\mbox{\scriptsize max}}}
    |\hat F(r)^q - F_\theta(r)^q|^p \, {\rm d}r
  }{
    D(theta) = integral from rmin to rmax of
    abs(Fhat(r)^q - F(theta,r)^q)^p
  }
  where \eqn{\theta}{theta} is the vector of parameters of the model,
  \eqn{\hat F(r)}{Fhat(r)} is the observed value of the summary statistic
  computed from the data, \eqn{F_\theta(r)}{F(theta,r)} is the
  theoretical expected value of the summary statistic,
  and \eqn{p,q} are two exponents. The default is \code{q = 1/4},
  \code{p=2} so that the contrast criterion is the integrated squared
  difference between the fourth roots of the two functions
  (Waagepetersen, 2007).

  The argument \code{action.bad.values} specifies what to do if
  some of the values of the summary statistic are \code{NA}, \code{NaN} or
  infinite. If \code{action.bad.values="stop"}, or if all of the values are bad,
  then a fatal error occurs. Otherwise, the domain of the summary
  function is shortened to avoid the bad values. The shortened domain is the
  longest interval on which the function values are finite 
  (provided this interval is at least half the length of the original
  domain). A warning is issued if \code{action.bad.values="warn"} (the default)
  and no warning is issued if \code{action.bad.values="silent"}.

  The other arguments just make things print nicely.
  The argument \code{fvlab} contains labels for the component
  \code{fit} of the return value.
  The argument \code{explain} contains human-readable strings
  describing the data, the model and the summary statistic.

  The \code{"..."} argument of \code{mincontrast} can be used to
  pass extra arguments to the function \code{theoretical}
  and/or to the optimisation function \code{\link[stats]{optim}}.
  In this case, the function \code{theoretical}
  should also have a \code{"..."} argument and should ignore it
  (so that it ignores arguments intended for \code{\link[stats]{optim}}).
}
\value{
  An object of class \code{"minconfit"}. There are methods for printing
  and plotting this object. It contains the following components:
  \item{par }{Vector of fitted parameter values.}
  \item{fit }{Function value table (object of class \code{"fv"})
    containing the observed values of the summary statistic
    (\code{observed}) and the theoretical values of the summary
    statistic computed from the fitted model parameters.
  }
  \item{opt }{The object returned from the optimizer \code{\link{optim}}.}
  \item{crtl }{List of parameters determining the contrast objective.}
  \item{info }{List of explanatory strings.}
}
\references{
  Diggle, P.J. and Gratton, R.J. (1984)
  Monte Carlo methods of inference for implicit statistical models.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{46}, 193 -- 212.

  \Moller, J. and Waagepetersen, R. (2003).
  Statistical Inference and Simulation for Spatial Point Processes.
  Chapman and Hall/CRC, Boca Raton.

  Pfanzagl, J. (1969).
  On the measurability and consistency of minimum contrast estimates.
  \emph{Metrika} \bold{14}, 249--276.

  Waagepetersen, R. (2007).
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}  
\author{
  \rasmus.
  Adapted for \pkg{spatstat} by
  \adrian.
}
\seealso{
  \code{\link{kppm}},
  \code{\link{lgcp.estK}},
  \code{\link{matclust.estK}},
  \code{\link{thomas.estK}},
}
\keyword{spatial}
\keyword{models}
