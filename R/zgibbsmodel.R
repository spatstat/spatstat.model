#'
#'    zgibbsmodel.R
#'
#'    Experimental
#' 
#'    $Revision: 1.3 $ $Date: 2022/04/26 07:56:22 $

zgibbsmodel <- local({

  validclass <- function(x) { is.numeric(x) || is.im(x) || is.function(x) }

  validtrend <- function(x) { validclass(x) || (is.list(x) && all(sapply(x, validclass))) }

  missingxy <- function(f) { is.function(f) && !all(c("x", "y") %in% names(formals(f))) }

  zgibbsmodel <- function(beta=1, interaction=NULL, icoef=NULL) {
    ## validate trend
    if(!validtrend(beta)) 
      stop("beta should be a number, a numeric vector, a pixel image, a function, or a list of such things")
    if(missingxy(beta))
      stop("Function beta() should have arguments x,y (and possibly others)")
    if(is.list(beta) && any(sapply(beta, missingxy)))
      stop("Each function beta() should have arguments x,y (and possibly others)")
    ## validate interaction
    if(is.null(interaction)) {
      interaction <- Poisson()
      if(length(icoef))
        stop("interaction coefficients icoef should not be specified for the Poisson process")
      icoef <- numeric(0)
    } else if(inherits(interaction, "fii")) {
      if(is.null(icoef)) {
        icoef <- coef(interaction)
      } else if(length(icoef) != length(coef(interaction)))
        stop("supplied interaction coefficients icoef have the wrong length")
      interaction <- as.interact(interaction)
    } else if(!inherits(interaction, "interact"))
      stop("Argument 'interaction' should be an object of class 'interact' or 'fii'")
    ## 
    if(anyNA(interaction$par))
      stop("All irregular parameters of the interaction must be supplied")
    ## build
    out <- list(beta        = beta,
                interaction = interaction,
                icoef       = icoef)
    class(out) <- c("zgibbsmodel", class(out))
    return(out)
  }

  zgibbsmodel
})

is.poisson.zgibbsmodel <- function(x) { is.null(x$interaction$family) }

is.stationary.zgibbsmodel <- function(x) { is.numeric(x$beta) }

as.interact.zgibbsmodel <- function(object) { object$interaction }

as.isf.zgibbsmodel <- function(object) { object$interaction$family }

interactionorder.zgibbsmodel <- function(object) { interactionorder(as.interact(object)) }

print.zgibbsmodel <- function(x, ...) {
  splat(if(is.stationary(x)) "Stationary" else "Non-stationary",
        if(is.poisson(x)) "Poisson" else "Gibbs",
        "point process model")
  beta <- x$beta
  tname <- if(is.poisson(x)) "Intensity" else "Trend"
  tcolon <- paste0(tname, ":")
  if(is.numeric(beta)) {
    if(length(beta) == 1) {
      splat(tcolon, "numeric value =", beta)
    } else {
      splat(tcolon, "numeric vector =", paren(paste(beta, collapse=" "), "["))
    }
  } else if(is.function(beta)) {
    splat(tname, "= function:")
    print(beta)
  } else if(is.im(beta)) {
    splat(tname, "= pixel image:")
    print(beta)
  } else {
    splat(tname, "= list:")
    print(beta)
  }
  if(!is.poisson(x)) {
    print(as.interact(x))
    splat("Iinteraction coefficients:")
    print(x$icoef)
  }
  invisible(NULL) 
}

fakefii <- function(model) {
  ## create a 'fake' fii object from a zgibbsmodel
  stopifnot(inherits(model, "zgibbsmodel"))
  inte <- as.interact(model)
  if(is.multitype(inte)) stop("Not implemented for multitype interactions")
  ## determine dimension of potential, etc
  fakePOT <- inte$pot(d=matrix(, 0, 0), par=inte$par)
  IsOffset <- attr(fakePOT, "IsOffset")
  fakePOT <- ensure3Darray(fakePOT)
  Vnames <- dimnames(fakePOT)[[3]]
  p <- dim(fakePOT)[3]
  if(sum(nzchar(Vnames)) < p)
    Vnames <- if(p == 1) "Interaction" else paste0("Interaction.", 1:p)
  if(length(IsOffset) < p)
    IsOffset <- logical(p)
  ## determine interaction coefficients
  icoef <- model$icoef
  if(!any(nzchar(names(icoef)))) names(icoef) <- Vnames
  ## create fake object
  fii(inte, icoef, Vnames, IsOffset)
}

## contributed by Frederic Lavancier (hacked by Adrian)

intensity.zgibbsmodel <- function(X, ..., approx=c("Poisson","DPP")) {
  approx <- match.arg(approx)
  
  fint <- fakefii(X)
  beta <- X$beta

  lambda <- switch(approx,
                   Poisson = PoisSaddle(beta, fint),
                   DPP     = DPPSaddle(beta, fint))
  
  return(lambda)
}
