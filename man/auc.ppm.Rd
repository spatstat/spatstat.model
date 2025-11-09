\name{auc.ppm}
\alias{auc.ppm}
\alias{auc.kppm}
\alias{auc.slrm}
\title{
  Area Under ROC Curve
}
\description{
  Compute the AUC (area under the Receiver Operating Characteristic
  curve) for a fitted point process model.
}
\usage{
\method{auc}{ppm}(X, \dots, subset=NULL)

\method{auc}{kppm}(X, \dots, subset=NULL)

\method{auc}{slrm}(X, \dots, subset=NULL)
}
\arguments{
  \item{X}{
    Fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"lppm"})
    or fitted spatial logistic regression model
    (object of class \code{"slrm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.model]{roc.ppm}},
    and arguments passed to \code{\link[spatstat.geom]{as.mask}}
    controlling the pixel resolution for calculations,
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, for which the AUC should be
    calculated.
  }
}
\details{
  This command computes the AUC, the area under the Receiver Operating
  Characteristic curve. The ROC itself is computed by \code{\link[spatstat.explore]{roc}}.

  For a fitted point process model \code{X},
  the AUC measures the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  Suppose \eqn{\lambda(u)}{\lambda(u)} is the intensity function of the model.
  The AUC is the probability that
  \eqn{\lambda(x_i) > \lambda(U)}{\lambda(x[i]) > \lambda(U)}.
  That is, AUC is the probability that a randomly-selected data point
  has higher predicted intensity than does a randomly-selected spatial
  location.
  The AUC is \bold{not} a measure of the goodness-of-fit of the model
  (Lobo et al, 2007).

  (For spatial logistic regression models (class \code{"slrm"})
  replace \dQuote{intensity} by \dQuote{probability of presence}
  in the text above.)

  The algorithm also calculates the theoretically expected AUC value
  for this model, as described in \rocpapercite.
}
\value{
  Numeric vector of length 2 giving the AUC value 
  and the theoretically expected AUC value for this model.
}
\references{
  \rocpaper.
  
  Lobo, J.M.,
  \Jimenez-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in 
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S. 
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \Birkhauser, Basel.
}
\author{
  \rocketAuthors.
}
\seealso{
  \code{\link[spatstat.explore]{roc}},
  \code{\link[spatstat.model]{roc.ppm}},
  \code{\link[spatstat.explore]{youden}}.
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  auc(fit)
  auc(swedishpines, "x")
}
\keyword{spatial}
