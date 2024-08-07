\name{simulate.ppm}
\alias{simulate.ppm}
\title{Simulate a Fitted Gibbs Point Process Model}
\description{
  Generates simulated realisations from a fitted Gibbs or Poisson
  point process model.
}
\usage{
  \method{simulate}{ppm}(object, nsim=1, ...,
                         singlerun = FALSE,
                         start = NULL,
                         control = default.rmhcontrol(object, w=w),
                         w = window, 
                         window = NULL, 
                         project=TRUE, new.coef=NULL,
                         verbose=FALSE, progress=(nsim > 1),
                         drop=FALSE)
}
\arguments{
  \item{object}{
    Fitted point process model.
    An object of class \code{"ppm"}.
  }
  \item{nsim}{
    Number of simulated realisations.
  }
  \item{singlerun}{
    Logical. Whether to generate the simulated realisations
    from a single long run of the Metropolis-Hastings algorithm
    (\code{singlerun=TRUE}) or from separate, independent runs of the
    algorithm (\code{singlerun=FALSE}, the default).
  }
  \item{start}{Data determining the initial state
    of the Metropolis-Hastings algorithm.  See
    \code{\link[spatstat.random]{rmhstart}} for description of these arguments.
    Defaults to \code{list(n.start=npoints(data.ppm(object)))}
    meaning that the initial state of the algorithm
    has the same number of points as the original dataset.
  }
  \item{control}{Data controlling the running of
    the Metropolis-Hastings algorithm.  See \code{\link[spatstat.random]{rmhcontrol}}
    for description of these arguments.
  }
  \item{w, window}{
    Optional. The window in which the model is defined.
    An object of class \code{"owin"}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link[spatstat.random]{rmhcontrol}},
    or to \code{\link[spatstat.random]{rmh.default}}, or to covariate functions in the model.
  }
  \item{project}{
    Logical flag indicating what to do if the fitted model is
    invalid (in the sense that the values of the fitted coefficients do not
    specify a valid point process).
    If \code{project=TRUE} the closest valid model will be simulated;
    if \code{project=FALSE} an error will occur.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    from \code{\link{rmh.ppm}}
    during the simulation of each point pattern.
  }
  \item{progress}{
    Logical flag indicating whether to print progress reports
    for the sequence of simulations.
  }
  \item{new.coef}{
    New values for the canonical parameters of the model.
    A numeric vector of the same length as \code{coef(object)}.
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE}, the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"ppm"} of fitted
  point process models.
  
  Simulations are performed by \code{\link{rmh.ppm}}.

  If \code{singlerun=FALSE} (the default), the simulated patterns are
  the results of independent runs of the Metropolis-Hastings
  algorithm. If \code{singlerun=TRUE}, a single long run of the
  algorithm is performed, and the state of the simulation is saved
  every \code{nsave} iterations to yield the simulated patterns.

  In the case of a single run, the behaviour is controlled
  by the parameters \code{nsave,nburn,nrep}. These 
  are described in \code{\link[spatstat.random]{rmhcontrol}}. They may be passed
  in the \code{\dots} arguments or included in \code{control}.
  It is sufficient to specify two
  of the three parameters \code{nsave,nburn,nrep}.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"ppp"}).
  It also belongs to the class \code{"solist"}, so that it can be
  plotted, and the class \code{"timed"}, so that the total computation
  time is recorded.
}
\examples{
  \testonly{op <- spatstat.options(rmh.nrep=10)}
  fit <- ppm(japanesepines, ~1, Strauss(0.1))
  simulate(fit, 2)
  simulate(fit, 2, singlerun=TRUE, nsave=1e4, nburn=1e4)
  \testonly{spatstat.options(op)}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{simulate.kppm}},
  \code{\link[stats]{simulate}}
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{models}
