\name{methods.slrm}
\alias{methods.slrm} %DoNotExport
\alias{formula.slrm}
\alias{update.slrm}
\alias{print.slrm}
\alias{summary.slrm}
\alias{terms.slrm}
\alias{labels.slrm}
\alias{deviance.slrm}
\title{
  Methods for Spatial Logistic Regression Models
}
\description{
  These are methods for the class \code{"slrm"}. 
}
\usage{
\method{formula}{slrm}(x, \dots)
\method{print}{slrm}(x, ...)
\method{summary}{slrm}(object, ...)
\method{terms}{slrm}(x, \dots)
\method{labels}{slrm}(object, \dots)
\method{deviance}{slrm}(object, \dots)
\method{update}{slrm}(object, \dots, evaluate = TRUE, env = parent.frame())
}
\arguments{
  \item{x,object}{
    An object of class \code{"slrm"},
    representing a fitted spatial logistic regression model.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{evaluate}{
    Logical value. If \code{TRUE}, evaluate the updated call to
    \code{slrm}, so that the model is refitted; 
    if \code{FALSE}, simply return the updated call.
  }
  \item{env}{
    Optional environment in which the model should be updated.
  }
}
\details{
  These functions are methods for the generic commands
  \code{\link[stats]{formula}},
  \code{\link[stats]{update}},
  \code{\link{print}},
  \code{\link{summary}},
  \code{\link[stats]{terms}},
  \code{\link{labels}} and 
  \code{\link[stats]{deviance}}
  for the class \code{"slrm"}.

  An object of class \code{"slrm"} represents a fitted
  spatial logistic regression model. It is obtained from \code{\link{slrm}}.
}
\value{
  See the help files for the corresponding generic functions.
}
\author{
  \adrian
}
\seealso{
  \code{\link{slrm}},  \code{\link{plot.slrm}},
  \code{\link{predict.slrm}}, \code{\link{simulate.slrm}},
  \code{\link{vcov.slrm}},
  \code{\link{coef.slrm}}.
}
\examples{
  fit <- slrm(redwood ~ x)
  coef(fit)
  formula(fit)
  tf <- terms(fit)  
  labels(fit)
  deviance(fit) 
}
\keyword{spatial}
\keyword{methods}
