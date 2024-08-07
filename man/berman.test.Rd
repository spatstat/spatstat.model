\name{berman.test.ppm} 
\alias{berman.test.ppm}
\title{Berman's Tests for Point Process Model} 
\description{
  Tests the goodness-of-fit of a Poisson point process model
  using methods of Berman (1986).
}
\usage{

\method{berman.test}{ppm}(model, covariate,
                         which = c("Z1", "Z2"),
               alternative = c("two.sided", "less", "greater"), ...)

}
\arguments{
  \item{model}{
    A fitted point process model (object of class \code{"ppm"} or
    \code{"lppm"}).
  }
  \item{covariate}{
    The spatial covariate on which the test will be based.
    An image (object of class \code{"im"})
    or a function.
  }
  \item{which}{
    Character string specifying the choice of test.
  }
  \item{alternative}{
    Character string specifying the alternative hypothesis.
  }
  \item{\dots}{
    Additional arguments controlling the pixel resolution
    (arguments \code{dimyx}, \code{eps} and \code{rule.eps}
    passed to \code{\link[spatstat.geom]{as.mask}})
    or other undocumented features.
  }
}
\details{
  These functions perform a goodness-of-fit test of a Poisson point
  process model fitted to point pattern data. The observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same values under the model,
  are compared using either of two test statistics
  \eqn{Z_1}{Z[1]} and \eqn{Z_2}{Z[2]} proposed by Berman (1986).
  The \eqn{Z_1}{Z[1]} test is also known as the
  Lawson-Waller test.

  The function \code{berman.test} is generic, with methods for
  point patterns (\code{"ppp"} or \code{"lpp"})
  and point process models (\code{"ppm"} or \code{"lppm"}).
  \itemize{
    \item 
    If \code{X} is a point pattern dataset (object of class
    \code{"ppp"} or \code{"lpp"}), then
    \code{berman.test(X, ...)} performs a goodness-of-fit test of the
    uniform Poisson point process (Complete Spatial Randomness, CSR)
    for this dataset.
    \item
    If \code{model} is a fitted point process model
    (object of class \code{"ppm"} or \code{"lppm"})
    then \code{berman.test(model, ...)} performs
    a test of goodness-of-fit for this fitted model. In this case,
    \code{model} should be a Poisson point process.
  }
  
  The test is performed by comparing the observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same covariate under the model.
  Thus, you must nominate a spatial covariate for this test.
  
  The argument \code{covariate} should be either a \code{function(x,y)}
  or a pixel image (object of class \code{"im"} containing the values
  of a spatial function.
  If \code{covariate} is an image, it should have numeric values,
  and its domain should cover the observation window of the
  \code{model}. If \code{covariate} is a function, it should expect
  two arguments \code{x} and \code{y} which are vectors of coordinates,
  and it should return a numeric vector of the same length
  as \code{x} and \code{y}.  

  First the original data point pattern is extracted from \code{model}.
  The values of the \code{covariate} at these data points are
  collected.

  Next the values of the \code{covariate} at all locations in the
  observation window are evaluated. The point process intensity
  of the fitted model is also evaluated at all locations in the window.

  \itemize{
    \item If \code{which="Z1"},
    the test statistic \eqn{Z_1}{Z[1]} is computed as follows.
    The sum \eqn{S} of the covariate values at all data
    points is evaluated. The predicted mean \eqn{\mu}{\mu} and variance
    \eqn{\sigma^2}{\sigma^2} of \eqn{S} are computed
    from the values of the covariate at all locations in the window.
    Then we compute \eqn{Z_1 = (S-\mu)/\sigma}{Z[1]=(S-\mu)/\sigma}.
    Closely-related tests were proposed independently 
    by Waller et al (1993) and Lawson (1993)
    so this test is often termed the 
    Lawson-Waller test in epidemiological literature.
    \item If \code{which="Z2"},
    the test statistic \eqn{Z_2}{Z[2]} is computed as follows.
    The values of the \code{covariate} at all locations in the
    observation window, weighted by the point process intensity,
    are compiled into a cumulative distribution function \eqn{F}.
    The probability integral transformation is then applied:
    the values of the \code{covariate} at the original data points
    are transformed by the predicted cumulative distribution function
    \eqn{F} into numbers between 0 and 1. If the model is correct,
    these numbers are i.i.d. uniform random numbers.
    The standardised sample mean of these numbers is the
    statistic \eqn{Z_2}{Z[2]}. 
  }
  In both cases the null distribution of the test statistic
  is the standard normal distribution, approximately.

  The return value is an object of class \code{"htest"} containing the
  results of the hypothesis test. The print method for this class
  gives an informative summary of the test outcome.
}
\value{
  An object of class \code{"htest"} (hypothesis test)
  and also of class \code{"bermantest"},
  containing the results of the test. The return value can be
  plotted (by \code{\link[spatstat.explore]{plot.bermantest}}) or printed
  to give an informative summary of the test.
}
\section{Warning}{
  The meaning of a one-sided test must be carefully scrutinised: see
  the printed output.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.explore]{cdf.test}},
  \code{\link[spatstat.explore]{quadrat.test}},
  \code{\link[spatstat.model]{ppm}}
}
\references{
  Berman, M. (1986)
  Testing for spatial association between a point process
  and another stochastic process.
  \emph{Applied Statistics} \bold{35}, 54--62.

  Lawson, A.B. (1993)
  On the analysis of mortality events around a
  prespecified fixed point.
  \emph{Journal of the Royal Statistical Society, Series A}
  \bold{156} (3) 363--377.

  Waller, L., Turnbull, B., Clark, L.C. and Nasca, P. (1992)
  Chronic Disease Surveillance and testing of
  clustering of disease and exposure: Application to
  leukaemia incidence and TCE-contaminated dumpsites
  in upstate New York.
  \emph{Environmetrics} \bold{3}, 281--300.
}
\examples{
   # Berman's data
   X <- copper$SouthPoints
   L <- copper$SouthLines
   D <- distmap(L, eps=1)
   # test of fitted model
   fit <- ppm(X ~ x+y)
   berman.test(fit, D)
}
\keyword{htest}
\keyword{spatial}

\concept{Goodness-of-fit}
