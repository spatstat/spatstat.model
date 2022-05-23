\name{emend.slrm}
\alias{emend.slrm}
\title{
  Force Spatial Logistic Regression Model to be Valid
}
\description{
  Ensures that a fitted spatial logistic regression
  specifies a well-defined model.
}
\usage{
\method{emend}{slrm}(object, \dots, fatal=FALSE, trace=FALSE)
}
\arguments{
  \item{object}{
    Fitted point process model (object of class \code{"slrm"}).
  }
  \item{\dots}{Ignored.}
  \item{fatal}{
    Logical value indicating whether to generate an error
    if the model cannot be projected to a valid model.
  }
  \item{trace}{
    Logical value indicating whether to print a trace
    of the decision process. 
  }
}
\details{
  \code{emend.slrm} is a method for the generic \code{\link{emend}},

  The purpose of the function is to ensure that a fitted model
  is valid.
  
  The model-fitting function \code{\link{slrm}}
  fits spatial logistic regression models to point pattern data.

  In some circumstances, the fitted model returned by \code{\link{slrm}} may not
  specify a well-defined model, because
  some of the fitted coefficients of the model
  may be \code{NA} or infinite values. 
  This usually occurs when the data are insufficient to estimate
  all the parameters. The model is said to be
  \emph{unidentifiable} or \emph{confounded}.

  The function \code{emend.slrm} modifies the model \code{object}
  so that the model is valid. It 
  identifies the terms in the model \code{object}
  that are associated with illegal parameter values (i.e. parameter
  values which are either \code{NA}, infinite, or outside their permitted
  range). It considers all possible sub-models of \code{object}
  obtained by deleting one or more
  of these terms. It identifies which of these submodels are valid,
  and chooses the valid submodel with the largest pseudolikelihood. The result
  of \code{emend.slrm} or \code{project.slrm} is the
  true maximum pseudolikelihood fit to the data.

  For large datasets or complex models, the algorithm used in
  \code{emend.slrm} may be time-consuming, because it takes time to
  compute all the sub-models. A faster, approximate
  algorithm can be applied by setting
  \code{spatstat.options(project.fast=TRUE)}. This produces a
  valid submodel, which may not be the maximum likelihood submodel.

  Use the function \code{\link{valid.slrm}} to check whether a fitted model
  object specifies a well-defined model.
}
\value{
  Another point process model (object of class \code{"slrm"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{slrm}},
  \code{\link{valid.slrm}},
  \code{\link{emend}},
  \code{\link{spatstat.options}}
}
\examples{
   fit <- slrm(redwood ~ x + I(x))
   coef(fit)
   fit2 <- emend(fit)
   coef(fit2)
}
\keyword{spatial}
\keyword{models}
