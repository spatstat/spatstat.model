\name{as.fv.kppm}
\alias{as.fv.kppm}
\alias{as.fv.dppm}
\alias{as.fv.minconfit}
\title{Convert Fitted Model To Class fv}
\description{
  Converts fitted model into a function table (an object of class \code{"fv"}).
}
\usage{
  \method{as.fv}{kppm}(x)

  \method{as.fv}{dppm}(x)

  \method{as.fv}{minconfit}(x)
}
\arguments{
  \item{x}{A fitted model which will be converted into a function table}
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
}
\details{
  The generic command \code{as.fv} converts data \code{x}, that
  could be interpreted as the values of a function,
  into a function value table (object of the class \code{"fv"}
  as described in \code{\link{fv.object}}). This object can then
  be plotted easily using \code{\link{plot.fv}}.

  Objects of class \code{"kppm"} (and related classes) represent
  a model that has been fitted to a dataset 
  by computing a summary function of the dataset and
  matching it to the corresponding summary function of the model.
  The methods for \code{as.fv} for classes \code{"kppm"},
  \code{"dppm"} and \code{"minconfit"} extract this information:
  the result is a function table containing the
  observed summary function and the best fit summary function.
}
\examples{
  as.fv(kppm(redwood))
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
\concept{Model diagnostics}
