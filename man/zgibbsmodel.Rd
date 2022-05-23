\name{zgibbsmodel}
\alias{zgibbsmodel}
\title{
  Gibbs Model 
}
\description{
  Experimental code. Creates an object representing a Gibbs point
  process model. Typically used for theoretical calculations about
  such a model.
}
\usage{
zgibbsmodel(beta = 1, interaction = NULL, icoef = NULL)
}
\arguments{
  \item{beta}{
    First order trend term. A numeric value, numeric vector,
    pixel image, function, or a list of such objects.
  }
  \item{interaction}{
    Object of class \code{"interact"} specifying the interpoint
    interaction structure, or \code{NULL} representing the Poisson process.
  }
  \item{icoef}{
    Numeric vector of coefficients for the interpoint interaction.
  }
}
\details{
  Experimental.
}
\value{
  Object belonging to the experimental class \code{zgibbsmodel}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{methods.zgibbsmodel}}
}
\examples{
  m <- zgibbsmodel(10, Strauss(0.1), -0.5)
}
\keyword{spatial}
\keyword{models}
