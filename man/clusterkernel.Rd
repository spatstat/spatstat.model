\name{clusterkernel}
\alias{clusterkernel.kppm}
\title{
  Extract Cluster Offspring Kernel
}
\description{
  Given a cluster point process model, this command
  returns the probability density of the cluster offspring.
}
\usage{
\method{clusterkernel}{kppm}(model, \dots)
}
\arguments{
  \item{model}{
    Cluster model. Either a fitted cluster or Cox model
    (object of class \code{"kppm"}), or a character string
    specifying the type of cluster model.
  }
  \item{\dots}{
    Parameter values for the model,
    when \code{model} is a character string.
  }
}
\details{
  Given a specification of a cluster point process model, this command
  returns a \code{function(x,y)} giving the two-dimensional
  probability density of the cluster offspring points assuming a cluster parent
  located at the origin.
}
\value{
  A function in the \R language with arguments \code{x,y,\dots}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{clusterfield}}, \code{\link[spatstat.core]{kppm}}
}
\examples{
  fit <- kppm(redwood ~ x, "MatClust")
  f <- clusterkernel(fit)
  f(0.05, 0.02)
}
\keyword{spatial}
