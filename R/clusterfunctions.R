## clusterfunctions.R
##
## Contains methods for the generic functions
##  - clusterkernel
##  - clusterfield
##  - clusterradius.
##
##   $Revision: 1.11 $  $Date: 2022/02/21 02:24:34 $
##

## The generic clusterkernel() is now in spatstat.random

clusterkernel.kppm <- function(model, ...) {
  kernelR <- Kpcf.kppm(model, what = "kernel")
  f <- function(x, y = 0, ...){
    kernelR(sqrt(x^2+y^2))
  }
  return(f)
}

## The generic clusterradius() is now in spatstat.random

clusterfield.kppm <- function(model, locations = NULL, ...) {
    f <- clusterkernel(model)
    if(is.null(locations)){
        if(!is.stationary(model))
            stop("The model is non-stationary. The argument ",
                 sQuote("locations"), " must be given.")
        locations <- centroid.owin(Window(model), as.ppp = TRUE)
    }
    clusterfield.function(f, locations, ..., mu = model$mu)
}

## The generic clusterradius is defined in spatstat.random

clusterradius.kppm <- function(model, ..., thresh = NULL, precision = FALSE){
    a <- list(model = model$clusters,
              thresh = thresh,
              precision = precision)
    a <- append(a, as.list(c(model$clustpar, model$clustargs)))
    do.call(clusterradius.character, a)
}

