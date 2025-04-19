#'
#'   clusterprocess.R
#'
#'   $Revision: 1.1 $ $Date: 2025/04/19 05:17:27 $
#'

#' clusterprocess() is defined in spatstat.random

#' The following methods are for generics defined in spatstat.model

is.poissonclusterprocess.clusterprocess <- function(model) { TRUE }

pcfmodel.clusterprocess <- function(model, ...) {
  p <- model$rules$pcf
  mpar <- model$par.idio
  other <- model$other
  f <- function(r) {
    as.numeric(do.call(p, c(list(par=mpar, rvals=r), other)))
  }
  return(f)
}

Kmodel.clusterprocess <- function(model, ...) {
  K <- model$rules$K
  mpar <- model$par.idio
  other <- model$other
  f <- function(r) {
    as.numeric(do.call(K, c(list(par=mpar, rvals=r), other)))
  }
  return(f)
}
