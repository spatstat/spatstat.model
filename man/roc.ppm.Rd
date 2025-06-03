\name{roc.ppm}
\alias{roc.ppm}
\alias{roc.slrm}
\alias{roc.kppm}
\title{
  Receiver Operating Characteristic For Fitted Point Process Model
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a fitted point process model.
}
\usage{
\method{roc}{ppm}(X, covariate=NULL,
                  \dots, baseline=NULL, high=TRUE,
                  method = "raw",
                  CI = "none", alpha=0.05,
                  leaveoneout=FALSE, subset=NULL)

\method{roc}{slrm}(X, covariate=NULL,
                  \dots, baseline=NULL, high=TRUE,
                  method = "raw",
                  CI = "none", alpha=0.05,
                  leaveoneout=FALSE, subset=NULL)

\method{roc}{kppm}(X, covariate=NULL,
                  \dots, baseline=NULL, high=TRUE,
                  method = "raw",
                  CI = "none", alpha=0.05,
                  leaveoneout=FALSE, subset=NULL)
}
\arguments{
  \item{X}{
    Fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"})
    or fitted spatial logistic regression model
    (object of class \code{"slrm"}).
  }
  \item{covariate}{
    Spatial covariate. Either a \code{function(x,y)},
    a pixel image (object of class \code{"im"}), or
    one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates.
    Traditionally omitted when \code{X} is a fitted model.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} controlling the
    pixel resolution for calculations.
  }
  \item{baseline}{
    Optional. A spatial object giving a baseline intensity.
    Usually a \code{function(x,y)} or
    a pixel image (object of class \code{"im"})
    giving the baseline intensity at
    any location within the observation window.
    Alternatively a point pattern (object of class \code{"ppp"})
    with the locations of the reference population.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
  \item{method}{
    The method or methods that should be used to estimate the ROC curve.
    A character vector: current choices are
    \code{"raw"}, \code{"monotonic"}, \code{"smooth"} and \code{"all"}.
    See Details.
  }
  \item{CI}{
    Character string (partially matched) specifying whether confidence
    intervals should be computed, and for which method.
    See Details.
  }
  \item{alpha}{
    Numeric value between 0 and 1. The confidence intervals will have
    confidence level \code{1-alpha}. The default gives 95\%
    confidence intervals.
  }
  \item{leaveoneout}{
    Logical value specifying (for \code{roc.ppm}, \code{roc.slrm},
    \code{roc.kppm} and \code{roc.lppm}) whether the fitted intensity
    of the model at each of the original data points should be computed
    by the leave-one-out procedure 
    (i.e. by removing the data point in question from the point pattern,
    re-fitting the model to the reduced point pattern, and computing the
    intensity of this modified model at the point in question)
    as described in \rocpapercite.
    It is also possible to specify \code{leaveoneout=c(TRUE,FALSE)}
    so that both versions are calculated.
  }
  \item{subset}{
    Optional. A spatial window (object of class \code{"owin"})
    specifying a subset of the data, from which the ROC should be
    calculated.
  }
}
\details{
  The generic function \code{roc}
  computes the Receiver Operating Characteristic (ROC)
  curve. The area under the ROC is computed by \code{\link[spatstat.explore]{auc}}.

  For a fitted model \code{X}
  and a spatial covariate \code{Z}, two ROC curves are computed and
  plotted together. The first curve is obtained by
  extracting the original point pattern data to which the model \code{X}
  was fitted, and computing the raw ROC curve as described above. 
  The second curve is the predicted value of the raw curve according to
  the model, as described in \rocpapercite.

  For a fitted model \code{X}, if \code{covariate}
  is missing or \code{NULL}, the default is to take the \code{covariate}
  to be the fitted probability of presence (if \code{X} is a spatial
  logistic regression model) or the fitted point process intensity
  (for other models).
  This is the standard use of the ROC curve in spatial ecology.
  The ROC shows the ability of the
  fitted model to separate the spatial domain
  into areas of high and low density of points.
  For a fitted spatial logistic regression model,
  similarly the ROC shows the ability of the fitted presence
  probabilities to segregate the spatial domain into pixels
  with high and low probability of presence.
  The ROC is \bold{not} a diagnostic for the goodness-of-fit of the model
  (Lobo et al, 2007).

  There are currently three methods to estimate the ROC curve:
  \describe{
    \item{\code{"raw"}}{
      uses the raw empirical spatial cummulative distribution function of the
      covariate.
    }
    \item{\code{"monotonic"}}{
      uses a monotonic regression to estimate the relation between the covariate
      and the point process intensity and then calculates the ROC from that.
      This corresponds to a either a convex minorant or a concave majorant of
      the raw ROC curve.
    }
    \item{\code{"smooth"}}{
      uses a smooth estimate of the relation between the covariate and the point
      process intensity and then calculates the ROC from that. See
      \code{\link[spatstat.explore]{roc.rhohat}} for details.
    }
    \item{\code{"all"}}{
      uses all of the above methods.
    }
  }

  If \code{CI} is one of the strings \code{'raw'},
  \code{'monotonic'} or \code{'smooth'}, then
  pointwise 95\% confidence intervals for the true ROC curve
  will be computed based on the \code{raw}, \code{monotonic} or
  \code{smooth} estimates, respectively.
  The confidence level is \code{1-alpha}, so that for example
  \code{alpha=0.01} would give 99\% confidence intervals.
  By default, confidence bands for the ROC curve are not computed.
}
\value{
  Function value table (object of class \code{"fv"})
  which can be plotted to show the ROC curve.
  Also belongs to class \code{"roc"}.
}
\references{
  \rocpaper.

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
  \rocketAuthors.
}
\seealso{
\code{\link[spatstat.explore]{roc}},
\code{\link[spatstat.linnet]{roc.lpp}},
\code{\link[spatstat.explore]{roc.rhohat}}.

\code{\link[spatstat.model]{addROC}},
\code{\link[spatstat.model]{dropROC}},
\code{\link[spatstat.model]{addapply}},
\code{\link[spatstat.model]{dropply}} for adding or removing variables from a model.

\code{\link[spatstat.explore]{auc}}
}
\examples{
  gold <- rescale(murchison$gold, 1000, "km")
  faults <- rescale(murchison$faults, 1000, "km")
  dfault <- distfun(faults)
  
  fit <- ppm(gold ~ dfault)
  plot(roc(fit, method="all"))

  ## compare raw estimates with and without using leave-one-out intensity
  plot(roc(fit, leaveoneout=c(FALSE, TRUE)))
}
\keyword{spatial}
\keyword{nonparametric}
