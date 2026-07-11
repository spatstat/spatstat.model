#
#  objsurf.R
#
#  surface of the objective function for an M-estimator
#
#  $Revision: 1.42 $ $Date: 2026/07/11 06:44:29 $
#

objsurf <- function(x, ...) {
  UseMethod("objsurf")
}

objsurf.kppm <- objsurf.dppm <- function(x, ...,
                                         ngrid=32,
                                         par1lim=NULL, par2lim=NULL,
                                         enclose=FALSE,
                                         ratio=1.5,
                                         verbose=TRUE) {
  ## history of function evaluations
  h <- attr(x, "h") 
  if(!is.null(h)) {
    if(enclose) {
      ## determine par1lim, par2lim to enclose the history
      if(is.null(par1lim)) par1lim <- range(h[,1L])
      if(is.null(par2lim)) par2lim <- range(h[,2L])
    }
    parmap <- resolve.parmap(...)
    if(!is.null(parmap)) {
      ## transform original parameter values to new scale/parametrisation
      tran1 <- parmap[[1L]]
      tran2 <- parmap[[3L]]
      h[,1L] <- suppressWarnings(tran1(h[,1L]))
      h[,2L] <- suppressWarnings(tran2(h[,2L]))
    }
    class(h) <- unique(c("traj", class(h)))
  }
  ## objective function surface
  Fit <- x$Fit
  switch(Fit$method,
         mincon = {
           result <- objsurf(Fit$mcfit, ..., 
                             ngrid=ngrid,
                             par1lim=par1lim, par2lim=par2lim,
                             ratio=ratio,
                             verbose=verbose)
         },
         palm = ,
         waag = ,
         clik2 = {
           optpar  <- x$par.canon %orifnull% x$par
           objfun  <- Fit$objfun
           objargs <- Fit$objargs
           result  <- objsurfEngine(objfun, optpar, objargs, ...,
                                    objname = "log composite likelihood",
                                    ngrid=ngrid,
                                    par1lim=par1lim, par2lim=par2lim,
                                    ratio=ratio,
                                    verbose=verbose)
         },
         adapcl = {
           stop(paste("Sorry, objective function is not available for",
                      "adaptive composite likelihood method"),
                call.=FALSE)
         },
         stop(paste("Unrecognised fitting method", dQuote(Fit$method)),
              call.=FALSE)
         )
  ## history of function evaluations
  attr(result, "trace") <- h
  return(result)
}

objsurf.minconfit <- function(x, ..., ngrid=32, par1lim=NULL, par2lim=NULL,
                              ratio=1.5, verbose=TRUE) {
  optpar  <- x$par.canon %orifnull% x$par
  objfun  <- x$objfun
  objargs <- x$objargs
  dotargs <- x$dotargs
  result <- objsurfEngine(objfun, optpar, objargs, ...,
                          objname = "contrast",
                          dotargs=dotargs,
                          ngrid=ngrid,
                          par1lim=par1lim, par2lim=par2lim,
                          ratio=ratio,
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
                          par1lim=NULL, par2lim=NULL,
                          ratio=1.5, verbose=TRUE) {
  trap.extra.arguments(...)
  if(!is.function(objfun))
    stop("Object is in an outdated format and needs to be re-fitted")
  if(length(new.objargs))
    objargs <- resolve.defaults(new.objargs, objargs)
  optpar <- unlist(optpar)
  npar   <- length(optpar)
  if(npar != 2)
    stop("Only implemented for functions of 2 arguments")
  ## specified range of parameter values
  if(!is.null(par1lim)) check.range(par1lim)
  if(!is.null(par2lim)) check.range(par2lim)
  ## possible transformation of parameters
  parmap <- resolve.parmap(parmap=parmap)
  ## create grid of values of (possibly transformed) parameters
  ratio <- ensure2vector(ratio)
  ngrid <- ensure2vector(ngrid)
  stopifnot(all(ratio > 1))
  if(untransformed <- is.null(parmap)) {
    ## parameter values will be evenly spaced on original parameter scale
    if(is.null(par1lim)) par1lim <- optpar[1L] * c(1/ratio[1L], ratio[1L])
    if(is.null(par2lim)) par2lim <- optpar[2L] * c(1/ratio[2L], ratio[2L])
    xgrid <- seq(par1lim[1L], par1lim[2L], length=ngrid[1L])
    ygrid <- seq(par2lim[1L], par2lim[2L], length=ngrid[2L])
    pargrid <- expand.grid(xgrid, ygrid)
    colnames(pargrid) <- names(optpar)
    optxy <- optpar
    par1 <- xgrid
    par2 <- ygrid
  } else {
    ## parameter values will be evenly spaced on transformed scale
    ## e.g. parmap=list(log, exp) 
    tran1 <- parmap[[1L]]  # transformation of first parameter
    inv1  <- parmap[[2L]]
    tran2 <- parmap[[3L]]  # transformation of second parameter
    inv2  <- parmap[[4L]]
    ## optimal parameters on transformed scale (x,y)
    optxy <- c(tran1(optpar[1L]), tran2(optpar[2L]))
    ## range of arguments on transformed scale (x,y), in neighbourhood of optxy
    if(is.null(par1lim)) {
      xlim <- sort(optxy[1L] * c(1/ratio[1L], ratio[1L]))
      par1lim  <- sort(inv1(xlim))
    } else {
      xlim <- sort(tran1(par1lim))
    }
    if(is.null(par2lim)) {
      ylim <- sort(optxy[2L] * c(1/ratio[2L], ratio[2L]))
      par2lim  <- sort(inv2(ylim))
    } else {
      ylim <- sort(tran2(par2lim))
    }
    ## equally spaced arguments on transformed scale (x,y)
    xgrid <- seq(xlim[1L], xlim[2L], length=ngrid[1L])
    ygrid <- seq(ylim[1L], ylim[2L], length=ngrid[2L])
    xygrid <- expand.grid(xgrid, ygrid)
    ## parameter values for evaluation
    par1 <- inv1(xgrid)
    par2 <- inv2(ygrid)
    pargrid <- cbind(inv1(xygrid[,1L]), inv2(xygrid[,2L]))
    colnames(pargrid) <- names(optpar)
  }
  # evaluate objective function
  if(verbose) cat(paste("Evaluating", nrow(pargrid), "function values..."))
  values <- do.call(apply,
                    append(list(pargrid, 1, objfun, objargs=objargs), dotargs))
  if(verbose) cat("Done.\n")
  ## save as objsurf object
  result <- objsurfObject(
    z       = matrix(values, ngrid[1], ngrid[2]),
    par1    = par1,
    par2    = par2,
    optpar  = optpar,
    optxy   = optxy,
    par1lim = par1lim,
    par2lim = par2lim,
    parmap  = parmap,
    x       = xgrid,
    y       = ygrid)
  return(result)
}

objsurfObject <- function(par1, par2, z,
                          ..., 
                          optpar  = NULL,
                          optxy   = NULL,
                          par1lim = NULL,
                          par2lim = NULL,
                          parmap  = NULL,
                          x       = NULL,
                          y       = NULL,
                          objname = "objective") {
  #' par1:
  #'    vector of length nrow(z) containing values of first parameter
  #' par2:
  #'    vector of length ncol(z) containing values of second parameter
  #' z:
  #'    a matrix containing the results of evaluating the objective function
  #'    at each combination of parameter values (par1[i], par2[j])
  #' optpar:
  #'    vector of length 2 giving the optimal pair of parameter values
  #' parmap:
  #'    either NULL or list of functions.
  #'    Transformations of the parameter values to a scale x, y
  #'    on which par1, par2 are evenly spaced. This enables plotting
  #'    using image.default, contour.default etc.
  #'    parmap = list(f, g) or list(f1, g1, f2, g2)
  #'    where f is the transformation theta -> x
  #'          g is the inverse x -> theta
  #'    For equal spacing on log scale, set parmap = list(log, exp)
  #' x, y:
  #'    optional vectors equivalent to par1, par2 after transformation
  #'    Used only for display
  stopifnot(is.matrix(z))
  stopifnot(length(par1) == nrow(z))
  stopifnot(length(par2) == ncol(z))
  parmap <- resolve.parmap(parmap=parmap)
  if(transformed <- !is.null(parmap)) {
    ## transformation of parameters
    tran1 <- parmap[[1L]]     # transformation of first parameter theta1 -> x
    inv1  <- parmap[[2L]]  # inverse x -> theta1
    tran2 <- parmap[[3L]]     # transformation of second parameter theta2 -> y
    inv2  <- parmap[[4L]]  # inverse y -> theta2
  } else {
    tran1 <- tran2 <- inv1 <- inv2 <- I
  }
  if(is.null(x)) {
    x <- tran1(par1)
  } else {
    stopifnot(length(x) == nrow(z))
  }
  if(is.null(y)) {
    y <- tran2(par2)
  } else {
    stopifnot(length(y) == ncol(z))
  }
  ## create object used for image display
  a <- list(x=x, y=y, z=z)
  ## x and y are evenly spaced (as required for image.default)
  ## but may be the result of transforming parameter values e.g.  x = log(theta)
  ## optpar are the optimal parameter values (untransformed scale)
  optpar <- unlist(optpar)
  stopifnot(length(optpar) == 2)
  if(is.null(optxy)) {
    ## optxy are the optimal parameter values (transformed scale, for plots)
    if(!transformed) {
      optxy <- optpar
    } else {
      optxy <- c(tran1(optpar[1L]), tran2(optpar[2L]))
    }
  }
  if(is.null(par1lim)) 
    par1lim <- if(!transformed) range(x) else range(inv1(x))
  if(is.null(par2lim)) 
    par2lim <- if(!transformed) range(y) else range(inv2(y))
  attr(a, "optpar")  <- optpar # optimal param values (original scale)
  attr(a, "optxy")   <- optxy  # optimal param values (transformed scale)
  attr(a, "objname") <- objname
  attr(a, "par1lim") <- par1lim   # original scale
  attr(a, "par2lim") <- par2lim
  attr(a, "parmap")  <- parmap  # list of functions, or NULL
  class(a) <- "objsurf"
  return(a)
}


print.objsurf <- function(x, ...) {
  cat("Objective function surface\n")
  optpar  <- attr(x, "optpar")
  objname <- attr(x, "objname")
  par1lim <- attr(x, "par1lim")
  par2lim <- attr(x, "par2lim") 
  if(is.null(par1lim) || is.null(par2lim)) {
    message("objsurf object is in an older format -- recovering")
    par1lim <- range(x$x) # older format did not support transformations
    par2lim <- range(x$y)
  }
  nama <- names(optpar)
  cat(paste("\tFunction value:", objname, "\n"))
  cat("Parameter limits:\n")
  cat(paste("\t", paste0(nama[1L], ":"), prange(signif(par1lim, 4)), "\n"))
  cat(paste("\t", paste0(nama[2L], ":"), prange(signif(par2lim, 4)), "\n"))
  if(!is.null(attr(x, "h")))
    splat("[Includes history of evaluations of objective function]")
  invisible(NULL)
}

summary.objsurf <- function(object, ...) {
  y <- list(par1lim   = attr(object, "par1lim"),
            par2lim   = attr(object, "par2lim"),
            objrange  = range(object$z),
            optpar    = as.list(attr(object, "optpar")),
            objname   = attr(object, "objname"),
            parmap    = attr(object, "parmap"), # could be NULL
            trace     = attr(object, "h") # may be NULL
            )
  if(is.null(y$par1lim) || is.null(y$par2lim)) {
    message("objsurf object is in an older format --- recovering")
    y$par1lim <- range(object$x) # older format did not support transforms
    y$par2lim <- range(object$y)
  }
  class(y) <- c("summary.objsurf", class(y))
  return(y)
}

print.summary.objsurf <- function(x, ...) {
  with(x, {
    cat("Objective function surface\n")
    cat(paste("\tFunction value:", objname, "\n"))
    cat(paste("\tRange of values:", prange(objrange), "\n"))
    cat("Parameter limits:\n")
    nama <- names(optpar)
    cat(paste("\t", paste0(nama[1L], ":"), prange(par1lim), "\n"))
    cat(paste("\t", paste0(nama[2L], ":"), prange(par2lim), "\n"))
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
  parmap <- attr(x, "parmap")
  if(untransformed <- is.null(parmap)) {
    xaxistype <- yaxistype <- "normal"
  } else {
    xaxistype <- if(identical(parmap[[1L]], log)) "log" else "handmade"
    yaxistype <- if(identical(parmap[[3L]], log)) "log" else "handmade"
  }
  xx <- unclass(x)
  dont.complain.about(xx)
  do.call(image,
          resolve.defaults(list(x=quote(xx), axes=untransformed),
                           list(...),
                           list(xlab=nama[1L], ylab=nama[2L], main=xname)))
  switch(xaxistype,
         normal = { },
         log = {
           rx10 <- range(x$x)/log(10)
           px <- axisTicks(rx10, log=TRUE)
           axis(1, at=log(px), labels=px)
         },
         handmade = {
           ## annotate axes on original scale
           tran1 <- parmap[[1L]]
           inv1 <- parmap[[2L]]
           px <- pretty(range(inv1(x$x)))
           axis(1, at=tran1(px), labels=px)
         })
  switch(yaxistype,
         normal = { },
         log = {
           ry10 <- range(x$y)/log(10)
           py <- axisTicks(ry10, log=TRUE)
           axis(2, at=log(py), labels=py)
         },
         handmade = {
           ## annotate axes on original scale
           tran2 <- parmap[[3L]]
           inv2 <- parmap[[4L]]
           py <- pretty(range(inv2(x$y)))
           axis(2, at=tran2(py), labels=py)
         })
  ## plot coordinates of optimal parameter values in transformed scale
  optxy <- attr(x, "optxy") %orifnull% optpar
  abline(v=optxy[1L], lty=3)
  abline(h=optxy[2L], lty=3)
  return(invisible(NULL))
}

lines.objsurf <- function(x, ..., directed=FALSE) {
  if(!is.null(h <- attr(x, "trace"))) {
    lines(h, ..., directed=directed) # trace is already in transformed coords
  } else message("No trajectory data")
  return(invisible(NULL))
}

contour.objsurf <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  optpar <- summary(x)[["optpar"]]
  parmap <- attr(x, "parmap")
  if(untransformed <- is.null(parmap)) {
    xaxistype <- yaxistype <- "normal"
  } else {
    xaxistype <- if(identical(parmap[[1L]], log)) "log" else "handmade"
    yaxistype <- if(identical(parmap[[3L]], log)) "log" else "handmade"
  }
  nama <- names(optpar)
  xx <- unclass(x)
  dont.complain.about(xx)
  do.call(contour,
          resolve.defaults(list(x=quote(xx), axes=untransformed),
                           list(...),
                           list(xlab=nama[1], ylab=nama[2], main=xname)))
  switch(xaxistype,
         normal = { },
         log = {
           rx10 <- range(x$x)/log(10)
           px <- axisTicks(rx10, log=TRUE)
           axis(1, at=log(px), labels=px)
         },
         handmade = {
           ## annotate axes on original scale
           tran1 <- parmap[[1L]]
           inv1 <- parmap[[2L]]
           px <- pretty(range(inv1(x$x)))
           axis(1, at=tran1(px), labels=px)
         })
  switch(yaxistype,
         normal = { },
         log = {
           ry10 <- range(x$y)/log(10)
           py <- axisTicks(ry10, log=TRUE)
           axis(2, at=log(py), labels=py)
         },
         handmade = {
           ## annotate axes on original scale
           tran2 <- parmap[[3L]]
           inv2 <- parmap[[4L]]
           py <- pretty(range(inv2(x$y)))
           axis(2, at=tran2(py), labels=py)
         })
  ## plot coordinates of optimal parameter values in transformed scale
  optxy <- attr(x, "optxy") %orifnull% optpar
  abline(v=optxy[1L], lty=3)
  abline(h=optxy[2L], lty=3)
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


resolve.parmap <- function(..., parmap=NULL) {
  if(is.null(parmap)) return(NULL)
  if(!(is.list(parmap) &&
       length(parmap) %in% c(2,4) &&
       all(sapply(parmap, is.function))))
    stop("Argument parmap should be a list of 2 or 4 functions", call.=FALSE)
  if(length(parmap) == 2)
    parmap <- rep(parmap, 2)
  return(parmap)
}
