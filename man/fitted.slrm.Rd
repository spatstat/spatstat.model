\name{fitted.slrm}
\Rdversion{1.1}
\alias{fitted.slrm}
\title{
  Fitted Probabilities for Spatial Logistic Regression
}
\description{
  Given a fitted Spatial Logistic Regression model,
  this function computes the fitted probabilities for each pixel,
  or the fitted probabilities at each original data point.
}
\usage{
  \method{fitted}{slrm}(object, \dots, type="probabilities",
                 dataonly=FALSE, leaveoneout=FALSE)
}
\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Ignored.
  }
  \item{type}{
    Character string (partially) matching one of 
    \code{"probabilities"}, \code{"intensity"} or \code{"link"}
    determining the quantity that should be predicted.
  }
  \item{dataonly}{
    Logical. If \code{TRUE}, then values will be computed
    at the points of the data point pattern. If \code{FALSE}, then
    values will be computed at the pixels used to fit the model.
  }
  \item{leaveoneout}{
    Logical value
    specifying whether to perform a leave-one-out calculation.
    See Details.
  }
}
\details{
  This is a method for the generic function
  \code{\link[stats:fitted.values]{fitted}}
  for spatial logistic regression models
  (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  By default, the fitted presence probabilities in each pixel are
  computed and returned as a pixel image.
  \itemize{
    \item If \code{dataonly=FALSE} (the default), a value is computed
    for each pixel in a fine grid.
    \item If \code{dataonly=TRUE}, a value is computed for each
    of the points in the original spatial point pattern.
  }
  If \code{leaveoneout=FALSE} (the default), fitted values are extracted
  from the fitted model \code{object}.
  If \code{leaveoneout=TRUE}, then a leave-one-out calculation is
  performed:
  \itemize{
    \item If \code{dataonly=FALSE} and \code{leaveoneout=TRUE},
    the predicted value in each \bold{pixel} is computed using a leave-one-out
    calculation. At each pixel \code{j}, the model is re-fitted using all the
    data \emph{except} the data inside pixel \code{j}, and this updated model is
    then predicted using the data inside pixel \code{j}.
    (\sQuote{Leave one pixel out})
    \item If \code{dataonly=TRUE} and  \code{leaveoneout=TRUE},
    the predicted value for each \bold{data point} \code{X[i]} is calculated
    by re-fitting the model to the data with \code{X[i]} removed,
    and then predicting this model at location \code{X[i]}.
    (\sQuote{Leave one point out})
  }
}
\value{
  A pixel image (object of class \code{"im"}) containing the
  fitted probability for each pixel,
  or a numeric vector containing the fitted probability at each data point.
}
\seealso{
  \code{\link{slrm}},
  \code{\link[stats:fitted.values]{fitted}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  plot(fitted(fit))
  fitted(fit, dataonly=TRUE)
}
\author{\spatstatAuthors.}
\keyword{spatial}
\keyword{models}
\keyword{methods}
