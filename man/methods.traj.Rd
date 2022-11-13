\name{methods.traj}
\alias{methods.traj} %DoNotExport
\alias{lines.traj}
\alias{plot.traj}
\alias{print.traj}
\title{
  Methods for Trajectories of Function Evaluations
}
\description{
  Methods for objects of class \code{"traj"}.
}
\usage{
\method{print}{traj}(x, \dots)
\method{plot}{traj}(x, \dots, show.ends=TRUE, add=FALSE, xlab=NULL, ylab=NULL)
\method{lines}{traj}(x, \dots, directed=FALSE)
}
\arguments{
  \item{x}{
    Object of class \code{"traj"}.
  }
  \item{\dots}{
    Additional arguments passed to other methods.
  }
  \item{directed}{
    Logical value specifying whether to draw arrows instead of
    undirected lines.
  }
  \item{show.ends}{
    Logical value specifying whether to indicate the start and finish of
    the trajectory. The start is a blue circle; the finish is a red cross.
  }
  \item{add}{
    Logical value specifying whether to draw the trajectory on the
    existing plot (\code{add=TRUE}) or to start a new plot
    (\code{add=FALSE}, the default).
  }
  \item{xlab,ylab}{Optional labels for the horizontal and vertical axes.}
}
\details{
  An object of class \code{"traj"} represents
  the history of evaluations of the objective function
  performed when a cluster process model was fitted.
  It is a data frame containing the input parameter values for
  the objective function, and the corresponding value of the objective
  function, that were considered by the optimisation algorithm.

  These functions are methods for the generic
  \code{print}, \code{plot} and \code{lines}.
}
\value{
  Null.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{traj}}
}
\examples{
  fit <- kppm(redwood, pspace=list(save=TRUE))
  h <- traj(fit)
  h
  plot(h)
  lines(h)
}
