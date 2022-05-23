\name{roc.ppm} 
\alias{roc.ppm}
\alias{roc.kppm}
\alias{roc.slrm}
\title{
  Receiver Operating Characteristic
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a point pattern or a fitted point process model.
}
\usage{

\method{roc}{ppm}(X, \dots)

\method{roc}{kppm}(X, \dots)

\method{roc}{slrm}(X, \dots)

}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"} or \code{"lpp"})
    or fitted point process model
    (object of class \code{"ppm"}, \code{"kppm"}, \code{"slrm"}
    or \code{"lppm"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} controlling the
    pixel resolution for calculations.
  }
}
\details{
  This command computes Receiver Operating
  Characteristic curve. The area under the ROC is computed by \code{\link{auc}}.

  For a point pattern \code{X} and a covariate \code{Z}, the
  ROC is a plot showing the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  For each possible threshold \eqn{z}, the algorithm calculates
  the fraction \eqn{a(z)} of area in the study region where the
  covariate takes a value greater than \eqn{z}, and the
  fraction \eqn{b(z)} of data points for which the covariate value
  is greater than \eqn{z}. The ROC is a plot of \eqn{b(z)} against
  \eqn{a(z)} for all thresholds \eqn{z}. 
  
  For a fitted point process model, 
  the ROC shows the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  The ROC is \bold{not} a diagnostic for the goodness-of-fit of the model
  (Lobo et al, 2007).

  (For spatial logistic regression models (class \code{"slrm"})
  replace \dQuote{intensity} by \dQuote{probability of presence}
  in the text above.)
}
\value{
  Function value table (object of class \code{"fv"})
  which can be plotted to show the ROC curve.
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
\code{\link{auc}}
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  plot(roc(fit))
}
\keyword{spatial}
