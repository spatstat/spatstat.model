\name{Kmodel.clusterprocess}
\alias{Kmodel.clusterprocess}
\alias{Lmodel.clusterprocess}
\alias{pcfmodel.clusterprocess}
\title{
  K Function or Pair Correlation Function of a Cluster Process Model
}
\description{
  Returns the theoretical \eqn{K} function, \eqn{L} function or
  pair correlation function of a cluster point process model
  with specified values of its parameters.
}
\usage{
 \method{Kmodel}{clusterprocess}(model, \dots)

 \method{Lmodel}{clusterprocess}(model, \dots)

 \method{pcfmodel}{clusterprocess}(model, \dots)
}
\arguments{
  \item{model}{
    Object of class \code{"clusterprocess"}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}
\value{
  A \code{function} in the \R language,
  which takes one argument \code{r}.
}
\details{
  For some point process models, it is possible to
  write down a mathematical expression for the \eqn{K} function
  or the pair correlation function of the model. 
  
  The commands \code{Kmodel}, \code{Lmodel} and \code{pcfmodel} are
  generic, with methods for the classes \code{"kppm"} (fitted cluster processes
  and Cox processes), \code{"ppm"} (fitted Gibbs processes),
  and \code{"clusterprocess"} (theoretical cluster processes).

  This page documents the methods for the class \code{"clusterprocess"}.
  An object of this class represents a
  Neyman-Scott-Cox cluster point process model with specified values
  of its parameters. 
  
  The return value is a \code{function} in the \R language,
  which takes one argument \code{r}.
  Evaluation of this function, on a numeric vector \code{r},
  yields values of the desired \eqn{K} function or pair correlation
  function at these distance values.
}
\seealso{
  \code{\link{clusterprocess}}.

  \code{\link[spatstat.explore]{Kest}}, \code{\link[spatstat.explore]{Lest}} or \code{\link[spatstat.explore]{pcf}}
  to estimate the \eqn{K} function, \eqn{L} function
  or pair correlation function nonparametrically from point pattern data.

  \code{\link{Kmodel.kppm}} for the method for fitted cluster processes
  and Cox processes.
  
  \code{\link{Kmodel.ppm}} for the method for fitted Gibbs processes.
}
\author{
  \adrian
}
\examples{
  m <- clusterprocess("Thomas", kappa=10, mu=5, scale=0.1)
  m2 <- clusterprocess("VarGamma", kappa=10, mu=10, scale=0.1, nu=0.7)
  m
  m2
  g <- pcfmodel(m)
  g(0.2)
  g2 <- pcfmodel(m2)
  g2(1)
}
\keyword{spatial}
\keyword{models}
