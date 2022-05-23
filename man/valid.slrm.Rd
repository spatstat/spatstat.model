\name{valid.slrm}
\alias{valid.slrm}
\title{
  Check Whether Spatial Logistic Regression Model is Valid
}
\description{
  Determines whether a fitted spatial logistic regression model
  is a well-defined model.
}
\usage{
  \method{valid}{slrm}(object, warn=TRUE, \dots)
}
\arguments{
  \item{object}{
    Fitted spatial logistic regression model (object of class \code{"slrm"}).
  }
  \item{warn}{
    Logical value indicating whether to issue a warning if the
    validity of the model cannot be checked (due to unavailability of
    the required code).
  }
  \item{\dots}{Ignored.}
}
\details{
  This is a method for the generic function \code{\link{valid}}
  for spatial logistic regression models (class \code{"slrm"}).
  
  In a model fitted by \code{\link{slrm}}, some of the fitted
  coefficients may be \code{NA} or infinite values.
  This can occur if the data are not
  adequate for estimation of the model parameters. The model is said
  to be \emph{unidentifiable} or \emph{confounded}.

  The function \code{valid.slrm} checks whether the fitted coefficients of
  \code{object} specify a well-defined model. It returns
  \code{TRUE} if the model is well-defined, and \code{FALSE} otherwise.

  Use the function \code{\link{emend.slrm}} to force the fitted model
  to be valid.
}
\value{
  A logical value, or \code{NA}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{slrm}},
  \code{\link{emend.slrm}}
}
\examples{
   fit1 <- slrm(cells ~ x)
   valid(fit1)
   fit2 <- slrm(cells ~ x + I(x))
   valid(fit2)
}
\keyword{spatial}
\keyword{models}
