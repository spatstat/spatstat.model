\name{measureWeighted}
\alias{measureWeighted}
\title{
  Weighted Version of a Measure
}
\description{
  Given a measure \code{m} (object of class \code{"msr"})
  and a spatially-varying weight function, construct the
  weighted version of \code{m}.
}
\usage{
measureWeighted(m, w)
}
\arguments{
  \item{m}{
    A measure (object of class \code{"msr"}).
  }
  \item{w}{
    A pixel image (object of class \code{"im"})
    or a \code{function(x,y)} giving the numeric weight at each spatial location.
  }
}
\details{
  For any region of space \code{B},
  the weighted measure \code{wm} has the value
  \deqn{
    wm(B) = \int_B w(x) dm(x)
  }{
    wm(B) = integral[B] w(x) dm(x)
  }
  In any small region of space, the increment of the
  weighted measure \code{wm} is equal to the increment of the 
  original measure \code{m} multiplied by the weight \code{w} at that
  location.
}
\value{
  Another measure (object of class \code{"msr"})
  on the same spatial domain.
}
\references{
  Halmos, P.R. (1950) \emph{Measure Theory}. Van Nostrand.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{msr}}, \code{\link{with.msr}}, \code{\link{split.msr}},
  \code{\link{measurePositive}}
}
\examples{
   fit <- ppm(cells ~ x)
   res <- residuals(fit)
   measureWeighted(res, function(x,y){x})
}
\keyword{spatial}
\keyword{math}
