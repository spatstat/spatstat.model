\name{auc.ppm}
\alias{auc.ppm}
\alias{auc.kppm}
\alias{auc.slrm}
\title{
  Area Under ROC Curve
}
\description{
  Compute the AUC (area under the Receiver Operating Characteristic
  curve) for
  a fitted point process model.
}
\usage{

\method{auc}{ppm}(X, \dots)

\method{auc}{kppm}(X, \dots)

\method{auc}{slrm}(X, \dots)

}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"} or \code{"lpp"})
    or fitted point process model (object of class \code{"ppm"},
    \code{"kppm"}, \code{"slrm"} or \code{"lppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} controlling the
    pixel resolution for calculations.
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
}
\value{
  Numeric.
  For \code{auc.ppm}, \code{auc.kppm} and \code{auc.lppm}, the result is a 
  numeric vector of length 2 giving the AUC value 
  and the theoretically expected AUC value for this model.
}
\references{
  Lobo, J.M.,
  \ifelse{latex}{\out{Jim{\'e}nez}}{Jimenez}-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in 
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S. 
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \ifelse{latex}{\out{Birkh{\"a}user}}{Birkhauser}, Basel.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.explore]{roc}}
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  auc(fit)
}
\keyword{spatial}
