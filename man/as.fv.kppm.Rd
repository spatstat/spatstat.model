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
  An object of class \code{"fv"} (see \code{\link[spatstat.explore]{fv.object}}).
}
\details{
  The generic command \code{\link[spatstat.explore]{as.fv}} converts data \code{x}, that
  could be interpreted as the values of a function,
  into a function value table (object of the class \code{"fv"}
  as described in \code{\link[spatstat.explore]{fv.object}}). This object can then
  be plotted easily using \code{\link[spatstat.explore]{plot.fv}}.

  Objects of class \code{"kppm"} and \code{"dppm"} represent
  a point process model that has been fitted to a point pattern dataset.
  By default, the fitting is performed by the method of moments,
  by computing a summary function of the dataset and
  matching it to the corresponding summary function of the model.
  The methods for \code{as.fv} for classes \code{"kppm"},
  \code{"dppm"} and \code{"minconfit"} extract this information:
  the result is a function table containing the
  observed summary function and the theoretical summary function
  for the model of best fit.

  If the object \code{x} was created by 
  \code{\link{kppm}} or \code{\link{dppm}}
  using a non-default value of the argument \code{method},
  then the model was fitted in a different way (namely by maximising a
  surrogate for the likelihood). In this case \code{as.fv(x)}
  computes an empirical estimate of the pair correlation function
  using \code{\link[spatstat.explore]{pcf}} or \code{\link[spatstat.explore]{pcfinhom}},
  and computes the theoretical pair correlation function of the fitted
  model using \code{\link{pcfmodel.kppm}}. Again these two functions are
  returned in a function table.
}
\examples{
  as.fv(kppm(redwood3))
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
\concept{Model diagnostics}
