#'  poissonfitsbetter.R
#'
#'  poisson.fits.better()
#'  and underlying calculations
#'
#'  $Revision: 1.1 $  $Date: 2022/11/10 06:51:16 $
#' 
#'  Copyright (c) Adrian Baddeley 2022
#'  GNU Public Licence >= 2.0

poisson.fits.better <- function(object) {
  if(!is.null(ispo <- object$ispo)) return(isTRUE(ispo))
  y <- PoissonCompareCalc(object)
  if(is.null(y)) return(FALSE)
  answer <- with(y, if(maximising) (poisval >= optval) else (poisval <= optval))
  return(answer)
}

PoissonCompareCalc <- function(object) {
  stopifnot(is.kppm(object))
  if(!isTRUE(object$isPCP)) return(NULL)
  Fit <- object$Fit
  switch(Fit$method,
         mincon = {
           m <- Fit$mcfit
           canonical <- !is.null(m$par.canon)
           optpar  <- if(canonical) m$par.canon else m$par
           objfun  <- m$objfun
           objargs <- m$objargs
           maximising <- FALSE
         },
         palm = ,
         clik2 = {
           canonical <- !is.null(object$par.canon)
           optpar  <- if(canonical) object$par.canon else object$par
           objfun  <- Fit$objfun
           objargs <- Fit$objargs
           maximising <- TRUE
         },
         return(NULL)
         )
  ## optimised value
  optval <- objfun(optpar, objargs=objargs)
  ## value for Poisson
  poispar <- optpar
  if(canonical) {
    if(is.na(match("strength", names(optpar))))
      stop("Internal error: the canonical parameters do not include 'strength'")
    poispar[["strength"]] <- 0
  } else {
    if(is.na(match("kappa", names(optpar))))
      stop("Internal error: the parameters do not include 'kappa'")
    poispar[["kappa"]] <- Inf
  }
  poisval <- objfun(poispar, objargs=objargs)
  return(list(optval=optval, poisval=poisval, maximising=maximising))
}
  
