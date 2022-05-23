\name{ic.kppm}
\alias{ic}
\alias{ic.ppm}
\alias{ic.kppm}
\title{Model selection criteria for the intensity function of a point process}
\description{
  Information criteria for selecting the intensity function model
  of a Poisson, cluster or Cox point process.
}
\usage{
  ic(object)

  \method{ic}{ppm}(object)

  \method{ic}{kppm}(object)
}
\arguments{
  \item{object}{
    Fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"}).
  }
}
\details{
  This function returns information criteria for selecting the
  intensity function model of a Poisson, Cox or cluster point process
  fitted by first order composite likelihood
  (i.e. using the Poisson likelihood function). 
  
  Degrees of freedom \eqn{df} for the information criteria are given by
  the trace of \eqn{S^{-1} \Sigma} where \eqn{S} is the sensitivity matrix
  and \eqn{\Sigma} is the variance matrix for the log composite
  likelihood score function. In case of a Poisson process, \eqn{df} is
  the number of parameters in the model for the intensity function. 
  
  The composite Bayesian information criterion (cbic) is
  \eqn{-2\ell + \log(n) df}{-2 * ll + log(n) * df}
  where \eqn{\ell}{ll} is the maximal log first-order composite likelihood
  (Poisson loglikelihood for the intensity function) and
  \eqn{n} is the observed number of points.
  It reduces to the BIC criterion in case of a Poisson process.
  
  The composite information criterion (cic) is
  \eqn{-2\ell + 2 df}{- 2 * ll + 2 * df}
  and reduces to the AIC in case of a Poisson process.
  
  NOTE: the information criteria are for selecting the intensity
  function model (a set of covariates) within a given model class.
  They cannot be used to choose among different types of cluster
  or Cox point process models (e.g. can not be used to choose
  between Thomas and LGCP models).
}
\value{
  A list with entries \code{loglike}, \code{cbic}, \code{cic} and
  \code{df}. 
  Here \code{loglike} is the fitted log first-order composite likelihood,
  \code{cbic} is composite Bayesian information criterion,
  \code{cic} is is the composite likelihood criterion and
  \code{df} is the adjusted degrees of freedom for the
  fitted intensity function model.
 } 
\seealso{
  \code{\link{kppm}}
}
\references{
  Choiruddin, A., Coeurjolly, J.F. and Waagepetersen, R. (2020) 
  Information criteria for inhomogeneous spatial point processes.
  \emph{Australian and New Zealand Journal of Statistics}.
  To appear. 
}
\examples{
  if(interactive()) {

    # model with one covariate
    fit1 <- kppm(bei~elev,data=bei.extra)
    ic1  <- ic(fit1)

    # model with two covariates
    fit2 <- kppm(bei~elev+grad,data=bei.extra)
    ic2  <- ic(fit2)
  
    # smallest cbic for fit1 but smallest cic for fit2
  }
}
\author{
  Achmad Choiruddin, Jean-Francois Coeurjolly and Rasmus Waagepetersen.
}
\keyword{spatial}
\keyword{models}
\concept{point process model}
\concept{Cox point process}
\concept{cluster process}
\concept{Neyman-Scott cluster process}

