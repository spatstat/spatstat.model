\name{is.poissonclusterprocess}
\alias{is.poissonclusterprocess}
\alias{is.poissonclusterprocess.kppm}
\alias{is.poissonclusterprocess.zclustermodel}
\alias{is.poissonclusterprocess.default}
\title{
  Recognise Poisson Cluster Process Models
}
\description{
  Given a point process model (either a model that has been fitted to
  data, or a model specified by its parameters), 
  determine whether the model is a Poisson cluster process.
}
\usage{
is.poissonclusterprocess(model)
\method{is.poissonclusterprocess}{kppm}(model)
\method{is.poissonclusterprocess}{zclustermodel}(model)
\method{is.poissonclusterprocess}{default}(model)
}
\arguments{
  \item{model}{
    Any kind of object representing a spatial point process model,
    either a model fitted to data, or a specification of a point process model.
  }
}
\details{
  The argument \code{model} represents a fitted spatial point process
  model (such as an object of class \code{"ppm"}, \code{"kppm"} or
  similar) or a specification of a point process model
  (such as an object of class \code{"zclustermodel"}).

  This function returns \code{TRUE} if the \code{model} is a Poisson cluster
  process, and \code{FALSE} otherwise.
  
  The function \code{is.poissonclusterprocess} is generic, with methods
  for classes \code{kppm} and \code{zclustermodel}, and a default
  method. 
}
\value{
  A logical value.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{kppm}}, \code{\link{zclustermodel}}.
}
\examples{
  fut <- kppm(redwood ~ 1, "Thomas")
  is.poissonclusterprocess(fut)
  fot <- slrm(cells ~ x, dimyx=16)
  is.poissonclusterprocess(fot)
}
\keyword{spatial}
\keyword{models}
