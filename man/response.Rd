\name{response}
\alias{response}
\alias{response.lm}
\alias{response.glm}
\alias{response.ppm}
\alias{response.dppm}
\alias{response.kppm}
\alias{response.slrm}
\alias{response.mppm}
\title{
  Extract the Values of the Response from a Fitted Model
}
\description{
  Given a fitted model (of any kind)
  extract the values of the response variable.
}
\usage{
response(object)

\method{response}{lm}(object)

\method{response}{glm}(object)

\method{response}{ppm}(object)

\method{response}{kppm}(object)

\method{response}{dppm}(object)

\method{response}{slrm}(object)

\method{response}{mppm}(object)
}
\arguments{
  \item{object}{
    A fitted model (object of class
    \code{"lm"},
    \code{"glm"},
    \code{"ppm"},
    \code{"kppm"},
    \code{"dppm"},
    \code{"slrm"} or
    \code{"mppm"}
    or some other class).
  }
}
\details{
  For fitted linear models of class \code{"lm"}
  and fitted generalized linear models of class \code{"glm"},
  the numerical values of the response variable are extracted
  if they are available, and otherwise \code{NULL} is returned.

  For fitted point process models of class
  \code{"ppm"},
  \code{"kppm"},
  \code{"dppm"},
  \code{"slrm"} or \code{"lppm"},
  the original data point pattern is extracted.

  For a fitted point process model of class \code{"mppm"},
  the list of original data point patterns is extracted.
}
\value{
  For \code{response.lm} and \code{response.glm},
  a numeric vector, or \code{NULL}.

  For
  \code{response.ppm},
  \code{response.kppm},
  \code{response.dppm}
  and \code{response.slrm}
  a two-dimensional spatial point pattern (class \code{"ppp"}).

  For \code{response.mppm}, a list of
  two-dimensional spatial point patterns (objects of class
  \code{"ppp"}).
  The list also belongs to classes \code{"solist"} and \code{"ppplist"}.
}
\author{
  \adrian.
}
\examples{
  fit <- ppm(cells ~ x)
  response(fit)  
}
\keyword{manip}
\keyword{models}
