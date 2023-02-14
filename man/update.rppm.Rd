\name{update.rppm}
\alias{update.rppm}
\title{Update a Recursively Partitioned Point Process Model}  
\description{
  \code{update} method for class \code{"rppm"}.
}
\usage{
 \method{update}{rppm}(object, \dots, envir=environment(terms(object)))
}
\arguments{
  \item{object}{
    Fitted recursively partitioned point process model.
    An object of class \code{"rppm"},
    obtained from \code{\link{rppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{rppm}}.
  }
  \item{envir}{
    Environment in which to re-evaluate the call to \code{\link{rppm}}.
  }
}
\details{
  \code{object} should be a fitted recursively partitioned point process model,
  obtained from the model-fitting function \code{\link{rppm}}.
  
  The model will be updated according to the new arguments provided.
}
\value{
  Another fitted recursively partitioned point process model (object of
  class \code{"rppm"}.
}
\seealso{
  \code{\link{rppm}}.
}
\examples{
 fit <- rppm(nztrees ~ x)
 newfit <- update(fit, . ~ x + y)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
