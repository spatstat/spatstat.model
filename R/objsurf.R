#
#  objsurf.R
#
#  surface of the objective function for an M-estimator
#
#  $Revision: 1.33 $ $Date: 2022/11/13 07:08:49 $
#

objsurf <- function(x, ...) {
  UseMethod("objsurf")
}

objsurf.kppm <- objsurf.dppm <- function(x, ...,
                                         ngrid=32,
                                         xlim=NULL, ylim=NULL,
                                         enclose=FALSE,
                                         ratio=1.5,
                                         verbose=TRUE) {
  ## history of function evaluations
  h <- attr(x, "h") 
  if(!is.null(h)) {
    dotargs <- list(...)
    if(!is.null(parmap <- dotargs$parmap)) {
      ## transform to new parametrisation
      tran <- parmap[[1]]
      tranpars <- as.data.frame(t(apply(h, 1, tran)))
      h <- cbind(tranpars, value=h$value)
    }
    if(enclose) {
      ## determine xlim, ylim to enclose the history
      if(is.null(xlim)) xlim <- range(h[,1])
      if(is.null(ylim)) ylim <- range(h[,2])
    }
    class(h) <- unique(c("traj", class(h)))
  }
  ## objective function surface
  Fit <- x$Fit
  switch(Fit$method,
         mincon = {
           result <- objsurf(Fit$mcfit, ..., 
                             ngrid=ngrid, xlim=xlim, ylim=ylim, ratio=ratio,
                             verbose=verbose)
         },
         palm = ,
         clik2 = {
           optpar  <- x$par.canon %orifnull% x$par
           objfun  <- Fit$objfun
           objargs <- Fit$objargs
           result  <- objsurfEngine(objfun, optpar, objargs, ...,
                                    objname = "log composite likelihood",
                                    ngrid=ngrid, xlim=xlim, ylim=ylim,
                                    ratio=ratio, verbose=verbose)
         },
         stop(paste("Unrecognised fitting method", dQuote(Fit$method)),
              call.=FALSE)
         )
  ## history of function evaluations
  attr(result, "h") <- h
  return(result)
}

objsurf.minconfit <- function(x, ..., ngrid=32, xlim=NULL, ylim=NULL,
                              ratio=1.5, verbose=TRUE) {
  optpar  <- x$par.canon %orifnull% x$par
  objfun  <- x$objfun
  objargs <- x$objargs
  dotargs <- x$dotargs
  result <- objsurfEngine(objfun, optpar, objargs, ...,
                          objname = "contrast",
                          dotargs=dotargs,
                          ngrid=ngrid, xlim=xlim, ylim=ylim, ratio=ratio,
                          verbose=verbose)
  return(result)
}

objsurfEngine <- function(objfun, optpar, objargs, 
                          ...,
                          dotargs=list(),
                          objname="objective", 
                          new.objargs=list(),
                          parmap = NULL,
                          ngrid=32,
                          xlim=NULL, ylim=NULL,
                          ratio=1.5, verbose=TRUE) {
  trap.extra.arguments(...)
  if(!is.function(objfun))
    stop("Object is in an outdated format and needs to be re-fitted")
  if(length(new.objargs))
    objargs <- resolve.defaults(new.objargs, objargs)
  npar    <- length(optpar)
  if(npar != 2)
    stop("Only implemented for functions of 2 arguments")
  ## create grid of values of (possibly transformed) parameters  
  ratio <- ensure2vector(ratio)
  ngrid <- ensure2vector(ngrid)
  stopifnot(all(ratio > 1))
  if(is.null(parmap)) {
    ## use original parameters
    if(is.null(xlim)) xlim <- optpar[1] * c(1/ratio[1], ratio[1])
    if(is.null(ylim)) ylim <- optpar[2] * c(1/ratio[2], ratio[2])
    xgrid <- seq(xlim[1], xlim[2], length=ngrid[1])
    ygrid <- seq(ylim[1], ylim[2], length=ngrid[2])
    pargrid <- expand.grid(xgrid, ygrid)
    colnames(pargrid) <- names(optpar)
  } else {
    if(!(length(parmap) == 2 && all(sapply(parmap, is.function))))
      stop("parmap should be a list of 2 functions")
    tran <- parmap[[1]]
    invtran <- parmap[[2]]
    ## transformed parameters
    Toptpar <- tran(optpar) 
    if(is.null(xlim)) xlim <- Toptpar[1] * c(1/ratio[1], ratio[1])
    if(is.null(ylim)) ylim <- Toptpar[2] * c(1/ratio[2], ratio[2])
    xgrid <- seq(xlim[1], xlim[2], length=ngrid[1])
    ygrid <- seq(ylim[1], ylim[2], length=ngrid[2])
    Tpargrid <- expand.grid(xgrid, ygrid)
    colnames(Tpargrid) <- names(Toptpar)
    ## inverse transform
    pargrid <- t(apply(Tpargrid, 1, invtran))
    colnames(pargrid) <- names(optpar)
    ## finally overwrite optimal values
    optpar <- Toptpar
  }
  # evaluate objective function
  if(verbose) cat(paste("Evaluating", nrow(pargrid), "function values..."))
  values <- do.call(apply,
                    append(list(pargrid, 1, objfun, objargs=objargs), dotargs))
  if(verbose) cat("Done.\n")
  result <- list(x=xgrid, y=ygrid, z=matrix(values, ngrid[1], ngrid[2]))
  attr(result, "optpar") <- optpar
  attr(result, "objname") <- objname
  class(result) <- "objsurf"
  return(result)
}

print.objsurf <- function(x, ...) {
  cat("Objective function surface\n")
  optpar <- attr(x, "optpar")
  objname <- attr(x, "objname")
  nama <- names(optpar)
  cat(paste("\tFunction value:", objname, "\n"))
  cat("Parameter limits:\n")
  cat(paste("\t", paste0(nama[1L], ":"), prange(signif(range(x$x), 4)), "\n"))
  cat(paste("\t", paste0(nama[2L], ":"), prange(signif(range(x$y), 4)), "\n"))
  if(!is.null(attr(x, "h")))
    splat("[Includes history of evaluations of objective function]")
  invisible(NULL)
}

summary.objsurf <- function(object, ...) {
  y <- list(xrange=range(object$x),
            yrange=range(object$y),
            objrange=range(object$z),
            optpar=as.list(attr(object, "optpar")),
            objname=attr(object, "objname")
            , trace=attr(object, "h") # may be NULL
            )
  class(y) <- c("summary.objsurf", class(y))
  return(y)
}

print.summary.objsurf <- function(x, ...) {
  with(x, {
    cat("Objective function surface\n")
    cat(paste("\tFunction value:", objname, "\n"))
    cat(paste("\tRange of values:", prange(objrange), "\n"))
    cat("Parameter limits (xrange, yrange):\n")
    nama <- names(optpar)
    cat(paste("\t", paste0(nama[1L], ":"), prange(xrange), "\n"))
    cat(paste("\t", paste0(nama[2L], ":"), prange(yrange), "\n"))
    cat("Selected parameter values (optpar):\n")
    cat(paste("\t", paste(nama, "=", optpar, collapse=", "), "\n"))
    if(!is.null(trace))
      splat("[Includes history of evaluations of objective function]")
  })
  return(invisible(NULL))
}


image.objsurf <- plot.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- attr(x, "optpar")
  nama <- names(optpar)
  xx <- unclass(x)
  dont.complain.about(xx)
  do.call(image,
          resolve.defaults(list(x=quote(xx)), 
                           list(...),
                           list(xlab=nama[1L], ylab=nama[2L], main=xname)))
  abline(v=optpar[1L], lty=3)
  abline(h=optpar[2L], lty=3)
  return(invisible(NULL))
}

lines.objsurf <- function(x, ..., directed=FALSE) {
  sx <- summary(x)
  if(!is.null(h <- sx[["trace"]])) {
    lines(h, ..., directed=directed)
  } else message("No trajectory data")
  return(invisible(NULL))
}

contour.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- summary(x)[["optpar"]]
  nama <- names(optpar)
  xx <- unclass(x)
  dont.complain.about(xx)
  do.call(contour,
          resolve.defaults(list(x=quote(xx)), 
                           list(...),
                           list(xlab=nama[1], ylab=nama[2], main=xname)))
  abline(v=optpar[1], lty=3)
  abline(h=optpar[2], lty=3)
  return(invisible(NULL))
}

  
persp.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- attr(x, "optpar")
  objname <- attr(x, "objname")
  nama <- names(optpar)
  xx <- x$x
  yy <- x$y
  zz <- x$z
  dont.complain.about(xx, yy, zz)
  r <- do.call(persp,
               resolve.defaults(list(x=quote(xx), y=quote(yy), z=quote(zz)),
                                list(...),
                                list(xlab=nama[1], ylab=nama[2],
                                     zlab=objname, main=xname)))
  return(invisible(r))
}


