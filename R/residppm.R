#
#  residppm.R
#
# computes residuals for fitted point process model
#
#
# $Revision: 1.27 $ $Date: 2023/02/02 02:37:29 $
#

residuals.ppm <-
  function(object, type="raw", ...,
           check=TRUE, drop=FALSE,
           fittedvalues = NULL,
           new.coef=NULL, dropcoef=FALSE,
           quad=NULL) {
  
  verifyclass(object, "ppm")
  trap.extra.arguments(..., .Context="In residuals.ppm")

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

  given.fitted <- !missing(fittedvalues) && !is.null(fittedvalues)

  # ................. determine fitted values .................

  NewCoef <- NULL
  if(is.null(new.coef) && is.null(quad)) {
    # use 'object' without modification
    # validate 'object'
    if(check && !given.fitted && damaged.ppm(object)) 
      stop("object format corrupted; try update(object, use.internal=TRUE)")
  } else {
    # determine a new set of model coefficients
    if(!is.null(new.coef)) {
      # use specified model parameters
      NewCoef <- new.coef
    } else {
      # estimate model parameters using a (presumably) denser set of dummy pts
      # Determine new quadrature scheme
      if(is.quad(quad))
        hi.res.quad <- quad
      else if(is.ppp(quad))
        hi.res.quad <- quadscheme(data=data.ppm(object), dummy=quad)
      else {
        # assume 'quad' is a list of arguments to 'quadscheme'
        hi.res.quad <- do.call(quadscheme,
                               append(list(data.ppm(object)),
                                      quad))
      }
      ## refit the model with new quadscheme
      e <- list2env(list(hi.res.quad=hi.res.quad), parent=object$callframe)
      hi.res.fit <- update(object, hi.res.quad, envir=e)
      NewCoef <- coef(hi.res.fit)
    }
  }
  #' now compute fitted values using new coefficients
  if(!given.fitted) 
    fittedvalues <- fitted(object, drop=drop, check=check,
                           new.coef=NewCoef, dropcoef=dropcoef)

  # ..................... compute residuals .....................

  # Extract quadrature points and weights
  Q <- quad.ppm(object, drop=drop, clip=drop)
#  U <- union.quad(Q) # quadrature points
  Z <- is.data(Q) # indicator data/dummy
#  W <- w.quad(Q) # quadrature weights

  # Compute fitted conditional intensity at quadrature points
  lambda <- fittedvalues

  # indicator is 1 if lambda > 0
  # (adjusted for numerical behaviour of predict.glm)
  indicator <- (lambda > .Machine$double.eps)

  if(type == "score") {
    # need the covariates
    X <- model.matrix(object)
    if(drop) {
      gs <- getglmsubset(object)
      ok <- !is.na(gs) & gs
      X <- X[ok, , drop=FALSE]
    }
  }
      
  # Evaluate residual measure components

  discrete <- switch(type,
                     raw     = rep.int(1, sum(Z)), 
                     inverse = 1/lambda[Z],
                     pearson = 1/sqrt(lambda[Z]),
                     score   = X[Z, , drop=FALSE]
                     )

  density <- switch(type,
                    raw     = -lambda,
                    inverse = -indicator,
                    pearson = -indicator * sqrt(lambda),
                    score   = -lambda * X)

  # Residual measure (return value)
  res <- msr(Q, discrete, density)

  # name the residuals
  attr(res, "type") <- type
  attr(res, "typename") <- typename

  return(res)
}

