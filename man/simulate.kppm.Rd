\name{simulate.kppm}
\alias{simulate.kppm}
\title{Simulate a Fitted Cluster Point Process Model}
\description{
  Generates simulated realisations from a fitted cluster point process model.
}
\usage{
  \method{simulate}{kppm}(object, nsim = 1, seed=NULL, ...,
         window=NULL, covariates=NULL,
         n.cond = NULL, w.cond = NULL,
         verbose=TRUE, retry=10,
         drop=FALSE)
}
\arguments{
  \item{object}{
    Fitted cluster point process model. An object of class \code{"kppm"}.
  }
  \item{nsim}{
    Number of simulated realisations.
  }
  \item{seed}{
    an object specifying whether and how to initialise
    the random number generator. Either \code{NULL} or an integer that will
    be used in a call to \code{\link[base:Random]{set.seed}}
    before simulating the point patterns. 
  }
  \item{\dots}{
    Additional arguments passed to the relevant random generator.
    See Details.
  }
  \item{window}{
    Optional. Window (object of class \code{"owin"}) in which the
    model should be simulated.
  }
  \item{covariates}{
    Optional. A named list containing new values for the covariates in the
    model. 
  }
  \item{n.cond}{
    Optional. Integer specifying a fixed number of points.
    See the section on \emph{Conditional Simulation}.
  }
  \item{w.cond}{
    Optional. Conditioning region. A window (object of class \code{"owin"})
    specifying the region which must contain exactly \code{n.cond} points.
    See the section on \emph{Conditional Simulation}.
  }
  \item{verbose}{
    Logical. Whether to print progress reports (when \code{nsim > 1}).
  }
  \item{retry}{
    Number of times to repeat the simulation if it fails
    (e.g. because of insufficient memory).
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE}, the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"kppm"} of fitted
  cluster point process models.
  
  Simulations are performed by
  \code{\link{rThomas}},
  \code{\link{rMatClust}},
  \code{\link{rCauchy}},
  \code{\link{rVarGamma}}
  or \code{\link{rLGCP}}
  depending on the model.

  Additional arguments \code{\dots} are passed to the relevant function
  performing the simulation.
  For example the argument \code{saveLambda} is recognised by all of the
  simulation functions.

  The return value is a list of point patterns.
  It also carries an attribute \code{"seed"} that
  captures the initial state of the random number generator.
  This follows the convention used in
  \code{simulate.lm} (see \code{\link[stats]{simulate}}).
  It can be used to force a sequence of simulations to be
  repeated exactly, as shown in the examples for \code{\link[stats]{simulate}}.
}
\section{Conditional Simulation}{
  If \code{n.cond} is specified, it should be a single integer.
  Simulation will be conditional on the event 
  that the pattern contains exactly \code{n.cond} points
  (or contains exactly \code{n.cond} points inside
  the region \code{w.cond} if it is given).

  Conditional simulation uses the rejection algorithm described
  in Section 6.2 of Moller, Syversveen and Waagepetersen (1998).
  There is a maximum number of proposals which will be attempted.
  Consequently the return value may contain fewer
  than \code{nsim} point patterns.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"ppp"}). (For conditional simulation,
  the length of the result may be shorter than \code{nsim}).

  The return value also carries an attribute \code{"seed"} that
  captures the initial state of the random number generator.
  See Details.
}
\examples{
  if(offline <- !interactive()) {
    spatstat.options(npixel=32, ndummy.min=16)
  }

  fit <- kppm(redwood ~x, "Thomas")
  simulate(fit, 2)

  simulate(fit, n.cond=60)

  if(offline) reset.spatstat.options()
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
  
  \Moller, J., Syversveen, A. and Waagepetersen, R. (1998)
  Log Gaussian Cox Processes.
  \emph{Scandinavian Journal of Statistics} \bold{25}, 451--482.
}
\seealso{
  \code{\link{kppm}},
  \code{\link{rThomas}},
  \code{\link{rMatClust}},
  \code{\link{rCauchy}},
  \code{\link{rVarGamma}},
  \code{\link{rLGCP}},
  \code{\link{simulate.ppm}},
  \code{\link[stats]{simulate}}
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{models}
