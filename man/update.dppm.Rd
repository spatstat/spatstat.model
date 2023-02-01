\name{update.dppm}
\alias{update.dppm}
\title{Update a Fitted Determinantal Point Process Model}
\description{
  \code{update} method for class \code{"dppm"}.
}
\usage{
 \method{update}{dppm}(object, \dots, evaluate=TRUE,
                       envir=environment(terms(object)))
}
\arguments{
  \item{object}{
    Fitted determinantal point process model.
    An object of class \code{"dppm"},
    obtained from \code{\link{dppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{dppm}}.
  }
  \item{evaluate}{
    Logical value indicating whether to return the updated fitted model
    (\code{evaluate=TRUE}, the default) or just the updated call to \code{dppm}
    (\code{evaluate=FALSE}).
  }
  \item{envir}{
    Environment in which to re-evaluate the call to \code{\link{dppm}}.
  }
}
\details{
  \code{object} should be a fitted determinantal point process model,
  obtained from the model-fitting function \code{\link{dppm}}.
  The model will be updated according to the new arguments provided.

  If the argument \code{trend} is provided, it determines the
  intensity in the updated model. It should be an \R formula
  (with or without a left hand side). It may include the symbols
  \code{+} or \code{-} to specify addition or deletion of terms
  in the current model formula, as shown in the Examples below.
  The symbol \code{.} refers to the current contents of the
  formula.

  The intensity in the updated model is determined by the
  argument \code{trend} if it is provided, or otherwise by any unnamed
  argument that is a formula, or otherwise by the formula of the
  original model, \code{formula(object)}.

  The spatial point pattern data to which the new model is fitted
  is determined by the left hand side of the updated model formula,
  if this is present. Otherwise it is determined by the argument
  \code{X} if it is provided, or otherwise by any unnamed argument
  that is a point pattern or a quadrature scheme.

  The model is refitted using \code{\link{dppm}}.
}
\value{
  Another fitted cluster point process model (object of
  class \code{"dppm"}.
}
\seealso{
  \code{\link{dppm}},  \code{\link{plot.dppm}},
  \code{\link{predict.dppm}}, \code{\link{simulate.dppm}},
  \code{\link{methods.dppm}}.
}
\examples{
  fit <- dppm(swedishpines ~ x + y, dppGauss, method="c")
  fitx <- update(fit, ~x)
  fit2 <- update(fit, flipxy(swedishpines))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
