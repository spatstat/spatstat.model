\name{predict.mppm}
\alias{predict.mppm}
\title{Prediction for Fitted Multiple Point Process Model}
\description{
    Given a fitted multiple point process model obtained by \code{\link{mppm}},
    evaluate the spatial trend and/or the conditional intensity of the
    model. By default, predictions are evaluated over a grid of
    locations, yielding pixel images of the trend and conditional intensity.
    Alternatively predictions may be evaluated at specified
    locations with specified values of the covariates.
}
\usage{
\method{predict}{mppm}(object, ..., newdata = NULL, type = c("trend", "cif"),
             ngrid = 40, locations=NULL, verbose=FALSE)
}
\arguments{
  \item{object}{The fitted model. An object of class \code{"mppm"}
    obtained from \code{\link{mppm}}.
  }
  \item{\dots}{Ignored.}
  \item{newdata}{
    Optional.
    New values of the covariates, for which the predictions should be computed.
    See Details.
  }
  \item{type}{
    Type of predicted values required. A character string or vector of
    character strings. Options are \code{"trend"} for the spatial trend
    (first-order term) and \code{"cif"} or \code{"lambda"} for the
    conditional intensity.
    Alternatively \code{type="all"} selects all options.
  }
  \item{ngrid}{
    Dimensions of the grid of spatial locations at which prediction will be
    performed (if \code{locations=NULL}). An integer or a pair of integers.
  }
  \item{locations}{
    Optional. The locations at which
    predictions should be performed. A list of point patterns, with one entry
    for each row of \code{newdata}.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
  }
}
\details{
  This function computes the spatial trend and the conditional
  intensity of a spatial point process model that has been fitted to
  several spatial point patterns. See Chapter 16 of
  Baddeley, Rubak and Turner (2015) for explanation and examples.
  
  Note that by ``spatial trend'' we mean the (exponentiated) first
  order potential and not the intensity of the process. [For example
  if we fit the stationary Strauss process with parameters
  \eqn{\beta}{beta} and \eqn{\gamma}{gamma},
  then the spatial trend is constant and equal to \eqn{\beta}{beta}.]
  The conditional intensity \eqn{\lambda(u,X)}{lambda(u,X)} of the fitted
  model is evaluated at each required spatial location u, with respect
  to the response point pattern X.

  If \code{newdata=NULL}, predictions are computed for the original
  values of the covariates, to which the model was fitted.
  Otherwise \code{newdata} should be a hyperframe
  (see \code{\link[spatstat.geom]{hyperframe}}) containing columns of covariates
  as required by the model. If \code{type} includes \code{"cif"},
  then \code{newdata} must also include a column of spatial point
  pattern responses, in order to compute the conditional intensity.

  If \code{locations=NULL}, then predictions are performed
  at an \code{ngrid} by \code{ngrid} grid of locations in the window
  for each response point pattern. The result will be a hyperframe
  containing a column of images of the trend (if selected)
  and a column of images of the conditional intensity (if selected).
  The result can be plotted.

  If \code{locations} is given, then it should be a list of point
  patterns (objects of class \code{"ppp"}). Predictions are performed at these
  points, and the results are returned as mark values attached to the
  \code{locations}. The result is a hyperframe containing columns
  called \code{trend} and/or \code{cif}. The column called \code{trend}
  contains marked point patterns in which the point locations are
  the \code{locations} and the mark value is the predicted trend.
  The column called \code{cif}
  contains marked point patterns in which the point locations are
  the \code{locations} and the mark value is the predicted conditional
  intensity.
}
\section{Models that depend on row number}{
  The point process model that is described by an \code{mppm} object
  may be a different point process for each row of the original
  hyperframe of data. This occurs if
  the model formula includes the variable \code{id} (representing
  row number) or if the model has a different interpoint interaction on
  each row.

  If the point process model
  is different on each row of the original data,
  then either
  \itemize{
    \item \code{newdata} is missing. Predictions are computed for
    each row of the original data using the point process model
    that applies on each row.
    \item \code{newdata} must have the same number of rows
    as the original data. Each row of \code{newdata} is assumed
    to be a replacement for the corresponding row of the original data.
    The prediction for row \code{i} of \code{newdata}
    will be computed for the point process model that applies to row \code{i}
    of the original data.
    \item \code{newdata} must include a column called \code{id}
    specifying the row number, and therefore identifying which
    of the point process models should apply.
    The predictions for row \code{i}
    of \code{newdata} will be computed for the point process model that applies
    to row \code{k} of the original data, where \code{k = newdata$id[i]}.
  }
}
\value{
  A hyperframe with columns named \code{trend} and/or \code{cif}.

  If  \code{locations=NULL}, the entries of the hyperframe are
  pixel images.

  If \code{locations} is not null, the entries are
  marked point patterns constructed by attaching the predicted values
  to the \code{locations} point patterns.
}
\references{
  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Baddeley, A., Bischof, L., Sintorn, I.-M., Haggarty, S.,
  Bell, M. and Turner, R. 
  Analysis of a designed experiment where the response is a spatial
  point pattern. In preparation.

  \baddrubaturnbook
}
\author{
  \adrian, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  \spatstatAuthors.
}
\seealso{
  \code{\link{mppm}},
  \code{\link{fitted.mppm}},
  \code{\link[spatstat.geom]{hyperframe}}
}
\examples{
  h <- hyperframe(Bugs=waterstriders)
  fit <- mppm(Bugs ~ x, data=h, interaction=Strauss(7))
  # prediction on a grid
  p <- predict(fit)
  plot(p$trend)
  # prediction at specified locations
  loc <- with(h, runifpoint(20, Window(Bugs)))
  p2 <- predict(fit, locations=loc)
  plot(p2$trend)
}
\keyword{spatial}
\keyword{models}
