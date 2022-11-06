\name{palmdiagnose}
\alias{palmdiagnose}
\title{
  Diagnostic based on Palm Intensity
}
\description{
  Given a fitted cluster process or Cox process model,
  calculate a diagnostic which compares nonparametric and parametric
  estimates of the Palm intensity.
}
\usage{
palmdiagnose(object, \dots, breaks = 30, trim = 30, rmax=Inf)
}
\arguments{
  \item{object}{
    Fitted model (object of class \code{"kppm"})
    or a list of fitted models.
  }
  \item{\dots}{
    Optional.
    Additional arguments which are fitted models of class \code{"kppm"}.
  }
  \item{breaks}{
    Optional argument passed to \code{\link[base]{cut.default}}
    determining the breakpoints of distance values 
    for the nonparametric estimate. Either an integer specifying the
    number of breakpoints, or a numeric vector of distance values to be
    used as the breakpoints.
  }
  \item{trim}{
    Optional. Maximum value of the translation edge correction weight.
  }
  \item{rmax}{
    Optional. Maximum interpoint distance \eqn{r} that should be
    considered. See Details.
  }
}
\details{
  This function computes the diagnostic proposed by
  Tanaka, Ogata and Stoyan (2008, Section 2.3)
  for assessing goodness-of-fit of a Neyman-Scott cluster process model
  to a point pattern dataset.

  The fitted model \code{object} should be an object of class \code{"kppm"}
  representing a Neyman-Scott cluster process model or a Cox process
  model. In the current implementation, the model must be stationary.

  The code computes parametric and non-parametric estimates of the
  Palm intensity \eqn{\lambda_0(r)}{\lambda[0](r)}, loosely speaking,
  the intensity of the point process given that there is a point at the origin.
  The parametric estimate is obtained from the fitted model by
  substituting the fitted parameter estimates into
  expressions for the pair correlation and the intensity.

  The non-parametric estimate is obtained by considering all pairs of
  data points, dividing the range of interpoint distances into
  several equally-spaced bands (determined by the argument
  \code{breaks}), counting the number of pairs of points whose
  interpoint distances fall in each band, and numerically adjusting for
  edge effects. Tanaka, Ogata and Stoyan (2008) used the
  periodic (toroidal) edge correction; our code uses the
  translation edge correction so that the method can be applied to
  data in any window.

  The result is a function value table (object of class \code{"fv"})
  containing the nonparametric and parametric estimates of the Palm
  intensity. The result also belongs to the class \code{"palmdiag"}
  which has a method for \code{plot}. The default behaviour of
  \code{\link{plot.palmdiag}} is to plot the model fit as a curve,
  and to display the nonparametric estimates as dots; this is the plot style
  proposed by Tanaka, Ogata and Stoyan (2008). Alternative display
  styles are also supported by \code{\link{plot.palmdiag}}.

  For computational efficiency, the argument \code{rmax}
  specifies the maximum value of interpoint distance \eqn{r}
  for which estimates of \eqn{\lambda_0(r)}{\lambda[0](r)} shall be computed.
  The default \code{rmax = Inf} implies
  there is no constraint on interpoint distance,
  and the resulting function object contains estimates of
  \eqn{\lambda_0(r)}{\lambda[0](r)} up to
  the maximum distance that would have been observable
  in the window containing the original point pattern data. 

  If there are additional arguments \code{\dots} which are fitted models
  of class \code{"kppm"}, or if \code{object} is a list of fitted models
  of class \code{"kppm"}, then the parametric estimates for each
  of the fitted models will be included in the resulting function object.
  If names are attached to these fitted models, the names will be used
  in the resulting function object.
}
\value{
  Function value table (object of class \code{"fv"}) containing the
  nonparametric and parametric estimates of the Palm intensity.
  Also belongs to the class \code{"palmdiag"} which has a \code{plot} method.
}
\author{
  \adrian.
}
\references{
  Tanaka, U., Ogata, Y. and Stoyan, D. (2008)
  Parameter estimation and model selection for Neyman-Scott Point
  Processes.
  \emph{Biometrical Journal} \bold{50}, 1, 43--57.
}
\examples{
   fitK <- kppm(redwood)
   R <- palmdiagnose(fitK)
   plot(R)

   fitg <- kppm(redwood, statistic="pcf")
   R2 <- palmdiagnose(A=fitK, B=fitg)
   plot(R2)
}
\seealso{
  \code{\link{plot.palmdiag}}
}
\keyword{spatial}
\keyword{models}
\keyword{nonparametric}
