#' residuals for exactppm
#'      (in spatstat.model because it uses 'msr')
#'
#' $Revision: 1.3 $ $Date: 2025/11/24 04:13:23 $

residuals.exactppm <- function(object, type="raw", ...) {
  verifyclass(object, "exactppm")
  type <- pickoption("type", type,
                     c(inverse="inverse",
                       raw="raw",
                       pearson="pearson",
                       Pearson="pearson",
                       score="score"))
  typenames <- c(inverse="inverse-lambda residuals",
                 raw="raw residuals",
                 pearson="Pearson residuals",
                 score="score residuals")
  typename <- typenames[[type]]

  ## Make a quadrature scheme
  
  Q <- quadscheme(object$X, ...)
  Z <- is.data(Q) # indicator data/dummy
  U <- union.quad(Q)

  # Compute fitted conditional intensity at quadrature points
  lambda <- fitted(object, locations=U)

  # ..................... compute residuals .....................

  # indicator is 1 if lambda > 0
  # (adjusted for numerical behaviour of predict.glm)
  indicator <- (lambda > .Machine$double.eps)

  # Evaluate residual measure components

  discrete <- switch(type,
                     raw     = rep.int(1, sum(Z)), 
                     inverse = 1/lambda[Z],
                     pearson = 1/sqrt(lambda[Z]),
                     score   = rep.int(1, sum(Z))
                     )

  density <- switch(type,
                    raw     = -lambda,
                    inverse = -indicator,
                    pearson = -indicator * sqrt(lambda),
                    score   = -lambda)

  # Residual measure (return value)
  res <- msr(Q, discrete, density)

  # name the residuals
  attr(res, "type") <- type
  attr(res, "typename") <- typename

  return(res)
}

