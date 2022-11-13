\name{plot.palmdiag}
\alias{plot.palmdiag}
\title{
  Plot the Palm Intensity Diagnostic
}
\description{
  Plot the Palm intensity diagnostic for a fitted cluster process
  or Cox process model.
}
\usage{
\method{plot}{palmdiag}(x, \dots,
          style = c("intervals", "dots", "bands"),
          args.dots = list(pch = 16), args.intervals = list(),
          xlim=NULL, main)
}
\arguments{
  \item{x}{
    Object of class \code{"palmdiag"} produced by \code{\link{palmdiagnose}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[spatstat.explore]{plot.fv}}
    when the fitted curve is plotted.
  }
  \item{style}{
    Character string specifying the style of plot for the
    nonparametric estimates. See Details.
  }
  \item{args.dots}{
    Arguments passed to \code{\link[graphics]{points}}
    when \code{style="dots"}.
  }
  \item{args.intervals}{
    Arguments passed to \code{\link[graphics]{segments}}
    when \code{style="intervals"}.
  }
  \item{xlim}{Optional range of distances plotted along the horizontal
    axis. A numeric vector of length 2.}
  \item{main}{Optional main title for plot.}
}
\details{
  This function plots the diagnostic proposed by
  Tanaka, Ogata and Stoyan (2008, Section 2.3)
  for assessing goodness-of-fit of a Neyman-Scott cluster process model
  to a point pattern dataset. The diagnostic is computed by
  the function \code{\link{palmdiagnose}}.

  First the Palm intensity of the fitted model is plotted as a
  function of interpoint distance \eqn{r} using
  \code{\link[spatstat.explore]{plot.fv}}. Then the nonparametric estimates
  of the Palm intensity are plotted on the same graph as follows:
  \itemize{
    \item if \code{style="dots"}, the nonparametric
    estimate for each band of distances is plotted as a dot,
    with horizontal coordinate at the middle of the band.
    This is the style proposed by Tanaka et al (2008).
    \item if \code{style="intervals"} (the default), each nonparametric estimate
    is plotted as a dot, and a 95\% confidence interval is
    plotted as a vertical line segment, centred on the dot.
    The confidence interval is based on the Poisson approximation.
    \item if \code{style="bands"}, the nonparametric estimates
    are drawn as a continuous curve which is flat on each band of
    distances. The 95\% confidence intervals are drawn as grey shading.
  }
}
\value{
  Null.
}
\references{
  Tanaka, U., Ogata, Y. and Stoyan, D. (2008)
  Parameter estimation and model selection for Neyman-Scott Point
  Processes.
  \emph{Biometrical Journal} \bold{50}, 1, 43--57.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{palmdiagnose}}
}
\examples{
   fit <- kppm(redwood)
   R <- palmdiagnose(fit)
   plot(R, style="d")
   plot(R)
   plot(R, style="b")
}
\keyword{models}
\keyword{nonparametric}
