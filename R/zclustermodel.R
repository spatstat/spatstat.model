#'
#'   zclustermodel.R
#'
#' Experimental
#' 
#'   $Revision: 1.10 $ $Date: 2022/02/21 02:24:07 $
#'

zclustermodel <- function(name="Thomas", ..., mu, kappa, scale) {
  if(missing(kappa)) stop("The parent intensity kappa must be given")
  if(missing(mu)) stop("The mean cluster size mu must be given")
  if(missing(scale)) stop("The cluster scale must be given")
  rules <- spatstatClusterModelInfo(name)
  par.std <- c(kappa=kappa, scale=scale)
  par.std <- rules$checkpar(par.std, native=FALSE)
  par.idio <- rules$checkpar(par.std, native=TRUE)
  other <- rules$resolveshape(...)
  clustargs <- rules$outputshape(other$margs)
  out <- list(name=name, rules=rules,
              par.std=par.std, par.idio=par.idio,
              mu=mu, clustargs=clustargs,
              other=other)
  class(out) <- "zclustermodel"
  return(out)
}

print.zclustermodel <- local({

  print.zclustermodel <- function(x, ...) {
    with(x, {
      splat(rules$printmodelname(list(clustargs=clustargs)))
      newpar <- rules$checkpar(par.std, native=FALSE)
      splat("Parent intensity kappa =", blurb("kappa", newpar["kappa"]))
      splat("Cluster scale = ", newpar["scale"])
      splat("Mean cluster size mu =", blurb("mu", mu))
      if(length(clustargs) > 0) {
        hdr <- paste("Cluster shape",
                     ngettext(length(clustargs), "parameter:", "parameters:"))
        if(is.list(clustargs) &&
           all(sapply(clustargs, is.numeric)) &&
           all(lengths(clustargs) == 1)) {
          splat(hdr,
                paste(names(clustargs), as.numeric(clustargs),
                      sep="=",
                      collapse=", "))
        } else {
          splat(hdr)
          print(clustargs)
        }
      }
    })
    return(invisible(NULL))
  }

  blurb <- function(name, value) {
    if(is.numeric(value)) as.character(value) else
    if(is.im(value)) "[image]" else "[unrecognized format]"
  }

  print.zclustermodel
})


                             
pcfmodel.zclustermodel <- function(model, ...) {
  p <- model$rules$pcf
  mpar <- model$par.idio
  other <- model$other
  f <- function(r) {
    do.call(p, c(list(par=mpar, rvals=r), other))
  }
  return(f)
}

Kmodel.zclustermodel <- function(model, ...) {
  K <- model$rules$K
  mpar <- model$par.idio
  other <- model$other
  f <- function(r) {
    as.numeric(do.call(K, c(list(par=mpar, rvals=r), other)))
  }
  return(f)
}

intensity.zclustermodel <- function(X, ...) {
  X$par[["kappa"]] * X$mu
}
  
predict.zclustermodel <- function(object, ...,
                                 locations,
                                 type="intensity",
                                 ngrid=NULL) {
  ## limited use!!!
  if(!identical(type, "intensity"))
    stop("Sorry, only type='intensity' is implemented")
  lambda <- object$par.std[["kappa"]] * object$mu
  if(is.numeric(lambda)) {
    if(is.ppp(locations))
      return(rep(lambda, npoints(locations)))
    W <- as.owin(locations)
    if(!is.mask(W))
      W <- as.mask(W, dimyx=ngrid, ...)
    return(as.im(lambda, W=W))
  }
  return(lambda[locations])
}

clusterradius.zclustermodel <- function(model, ...,
                                        thresh = NULL, precision = FALSE) {
  do.call(clusterradius.character,
          resolve.defaults(
            list(model = model$name, thresh = thresh, precision = precision),
            list(...),
            as.list(model$par.std), # sic
            model$clustargs)
          )
}

reach.zclustermodel <- function(x, ..., epsilon) {
  thresh <- if(missing(epsilon)) NULL else epsilon
  2 * clusterradius(x, ..., thresh=thresh)
}
