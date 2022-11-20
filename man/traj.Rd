\name{traj}
\alias{traj}
\title{
  Extract trajectory of function evaluations
}
\description{
  Extract the history of evaluations of the objective function
  performed when a cluster process model was fitted.
}
\usage{
traj(object)
}
\arguments{
  \item{object}{
    Fitted cluster point process model (object of class \code{"kppm"})
    or objective function surface (object of class \code{"objsurf"}).
  }
}
\details{
  Under appropriate circumstances, the fitted model \code{object}
  contains the history of evaluations of the objective function
  that were performed by the optimisation algorithm. This history is
  extracted by \code{traj}.

  The result is a data frame containing the input parameter values for
  the objective function, and the corresponding value of the objective
  function, that were considered by the optimisation algorithm. This
  data frame also belongs to the class \code{"traj"} which has
  methods for \code{plot}, \code{print} and other purposes.
}
\value{
  Either a data frame (belonging to class \code{"traj"})
  or \code{NULL}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{methods.traj}}
}
\examples{
  fit <- kppm(redwood, trajectory=TRUE)
  h <- traj(fit)
}

