#'  kppmComLik.R
#'
#'      Fitting algorithm for kppm(method = 'clik2')
#' 
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#'      C  o  m  p  o  s  i  t  e    L  i  k  e  l  i  h  o  o  d
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#'
#' Implementation of Guan (2006) second order composite likelihood
#' 
#'  Guan, Y. (2006) A composite likelihood approach
#'  in fitting spatial point process models.
#'  Journal of the American Statistical Association 101, 1502-1512.
#'
#'  $Revision: 1.3 $ $Date: 2026/01/29 01:08:50 $
#'
#'  Copyright (c) 2001-2025 Adrian Baddeley, Rolf Turner, Ege Rubak,
#'                Abdollah Jalilian and Rasmus Plenge Waagepetersen
#'  GNU Public Licence (>= 2.0)

kppmComLik <- function(X, Xname, po, clusters, control=list(), stabilize=TRUE,
                       weightfun, rmax, algorithm="Nelder-Mead",
                       DPP=NULL, ..., 
                       pspace=NULL) {
  pspace <- do.call(make.pspace,
                    resolve.defaults(
                      list(fitmethod="clik2", clusters=clusters),
                      list(...), ## ellipsis arguments override pspace
                      as.list(pspace),
                      .MatchNull=FALSE))
  W <- as.owin(X)
  if(is.null(rmax))
    rmax <- rmax.rule("K", W, intensity(X))
  ## identify unordered pairs of points that contribute
  cl <- closepairs(X, rmax, what="ijd", twice=FALSE, neat=FALSE)
  dIJ <- cl$d
  # compute weights for unordered pairs of points
  if(is.function(weightfun)) {
    wIJ <- weightfun(dIJ)
    sumweight <- safePositiveValue(sum(wIJ))
  } else {
    npairs <- length(dIJ)
    wIJ <- rep.int(1, npairs)
    sumweight <- npairs
  }
  # convert window to mask, saving other arguments for later
  dcm <- do.call.matched(as.mask,
                         append(list(w=W), list(...)),
                         sieve=TRUE)
  M         <- dcm$result
  otherargs <- dcm$otherargs

  ## Detect DPP usage
  isDPP <- inherits(clusters, "detpointprocfamily")

  # compute intensity at pairs of data points
  # and c.d.f. of interpoint distance in window
  if(stationary <- is.stationary(po)) {
    # stationary unmarked Poisson process
    lambda <- intensity(X)
#    lambdaIJ <- lambda^2
    # compute cdf of distance between two uniform random points in W
    g <- distcdf(W, delta=rmax/4096)
    # scaling constant is (area * intensity)^2
    gscale <- npoints(X)^2  
  } else {
    # compute fitted intensity at data points and in window
#    lambdaX <- fitted(po, dataonly=TRUE)
    lambda <- lambdaM <- predict(po, locations=M)
    # lambda(x_i) * lambda(x_j)
#    lambdaIJ <- lambdaX[I] * lambdaX[J]
    # compute cdf of distance between two random points in W
    # with density proportional to intensity function
    g <- distcdf(M, dW=lambdaM, delta=rmax/4096)
    # scaling constant is (integral of intensity)^2
    gscale <- safePositiveValue(integral.im(lambdaM)^2, default=npoints(X)^2)
  }

  # Detect DPP model and change clusters and intensity correspondingly
  isDPP <- !is.null(DPP)
  if(isDPP){
    tmp <- dppmFixIntensity(DPP, lambda, po)
    clusters <- tmp$clusters
    lambda <- tmp$lambda
    po <- tmp$po
  }

  # trim 'g' to [0, rmax] 
  g <- g[with(g, .x) <= rmax,]
  # get pair correlation function (etc) for model
  info <- spatstatClusterModelInfo(clusters)
  pcfun        <- info$pcf
  selfstart    <- info$selfstart
  isPCP        <- info$isPCP
  resolveshape <- info$resolveshape
  modelname    <- info$modelname
  # Assemble information required for computing pair correlation
  if(is.function(resolveshape)) {
    # Additional 'shape' parameters of cluster model are required.
    # These may be given as individual arguments,
    # or in a list called 'covmodel'
    clustargs <- if("covmodel" %in% names(otherargs))
                 otherargs[["covmodel"]] else otherargs
    shapemodel <- do.call(resolveshape, clustargs)$covmodel
  } else shapemodel <- NULL
  pcfunargs <- shapemodel
  margs <- pcfunargs$margs
  # determine starting parameter values
  startpar <- selfstart(X)
  
  pspace.given <- pspace
  #' ............ permit negative parameter values ...................
  strict <- !isFALSE(pspace$strict)
  if(!strict) pcfunargs <- append(pcfunargs, list(strict=FALSE))
  #' ............ parameter corresponding to Poisson process .........
  if(isDPP) {
    poispar <- NULL
  } else if(isPCP) {
    if(!("kappa" %in% names(startpar)))
      stop("Internal error: startpar does not include 'kappa'")
    poispar <- startpar
    poispar[["kappa"]] <- Inf
  } else {
    ## LGCP
    if(!("sigma2" %in% names(startpar)))
      stop("Internal error: startpar does not include 'sigma2'")
    poispar <- startpar
    poispar[["sigma2"]] <- .Machine$double.eps # i.e. 0
  }
  #' ............ use canonical parameters .........................
  usecanonical <- isTRUE(pspace$canonical)
  if(usecanonical) {
     tocanonical <- info$tocanonical
     tohuman <- info$tohuman
     if(is.null(tocanonical) || is.null(tohuman)) {
       warning("Canonical parameters are not yet supported for this model")
       usecanonical <- FALSE
     }
  }
  startpar.human <- startpar
  poispar.human <- poispar
  if(usecanonical) {
    pcftheo <- pcfun
    startpar <- tocanonical(startpar, margs=margs)
    if(!is.null(poispar)) poispar <- tocanonical(poispar, margs=margs)
    pcfun <- function(par, ...) { pcftheo(tohuman(par, ...), ...) }
  }
  #' ............ penalty ......................................
  penalty    <- pspace$penalty
  penal.args <- pspace$penal.args
  tau        <- pspace$tau %orifnull% 1
  if(is.function(penalty)) {
    ## penalised optimisation
    if(usecanonical) {
      penalty.human <- penalty
      penalty <- function(par, ...) { penalty.human(tohuman(par, ...), ...) }
    }
    ## data-dependent arguments in penalty
    if(is.function(penal.args)) 
      penal.args <- penal.args(X)
    ## exchange rate (defer evaluation if it is a function)
    if(!is.function(tau)) check.1.real(tau)
    ## reinsert in 'pspace' for insurance
    pspace$penalty <- penalty
    pspace$penal.args <- penal.args
    pspace$tau <- tau
  }
  #' ............ debugger .....................................
  TRACE <- isTRUE(pspace$trace)
  if(SAVE <- isTRUE(pspace$save)) {
    saveplace <- new.env()
    assign("h", NULL, envir=saveplace)
  } else saveplace <- NULL
  
  # .....................................................
  # create local function to evaluate pair correlation
  #  (with additional parameters 'pcfunargs' in its environment)
  paco <- function(d, par) {
    do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
  }
  #' ..........  define objective function ......................
  if(!is.function(weightfun)) {
    # pack up necessary information
    objargs <- list(dIJ=dIJ, sumweight=sumweight, g=g, gscale=gscale, 
                    envir=environment(paco),
                    ## The following variables are temporarily omitted
                    ## in order to calculate the objective function
                    ## without using them, or their side effects.
                    penalty=NULL,   # updated below
                    penal.args=NULL,   # updated below
                    tau=NULL,   # updated below
                    TRACE=FALSE, # updated below
                    saveplace=NULL, # updated below
                    BIGVALUE=1, # updated below
                    SMALLVALUE=.Machine$double.eps)
    # define objective function (with 'paco' in its environment)
    # This is the log composite likelihood minus the constant term 
    #       2 * (sum(log(lambdaIJ)) - npairs * log(gscale))
    obj <- function(par, objargs) {
      with(objargs, {
        logprod <- sum(log(safePositiveValue(paco(dIJ, par))))
        integ <- unlist(stieltjes(paco, g, par=par))
        integ <- pmax(SMALLVALUE, integ)
        logcl <- 2*(logprod - sumweight * log(integ))
        logcl <- safeFiniteValue(logcl, default=-BIGVALUE)
        ## penalty
        if(hasPenalty <- is.function(penalty)) {
          straf <- do.call(penalty, append(list(par), penal.args))
          logclPEN <- logcl - tau * straf
        }
        ## debugger
        if(isTRUE(TRACE)) {
          cat("Parameters:", fill=TRUE)
          print(par)
          splat("\tlogprod:", logprod)
          splat("\tinteg:", integ)
          splat("log composite likelihood:", logcl)
          if(hasPenalty)
            splat("penalised log composite likelihood:", logclPEN)
        }
        ## save state
        if(is.environment(saveplace)) {
          h <- get("h", envir=saveplace)
          value <- list(logcl=logcl)
          if(hasPenalty)
            value <- append(value, list(logclPEN=logclPEN))
          hplus <- as.data.frame(append(par, value))
          h <- rbind(h, hplus)
          assign("h", h, envir=saveplace)
        }
        return(if(hasPenalty) logclPEN else logcl)
      },
      enclos=objargs$envir)
    }
    ## Determine the values of some parameters
    ## (1) Determine a suitable large number to replace Inf
    objargs$BIGVALUE <- bigvaluerule(obj, objargs, startpar)
    ## (2) Evaluate exchange rate 'tau'
    if(is.function(penalty) && is.function(tau)) {
      ## data-dependent exchange rate 'tau': evaluate now
      if("poisval" %in% names(formals(tau))) {
        ## New style: requires value of (unpenalised) objective function at Poisson process
        poisval <- obj(poispar, objargs)
        tau <- tau(X, poisval=poisval)
      } else {
        tau <- tau(X)
      }
      check.1.real(tau)
    }
    ## Now insert the penalty, etc
    objargs <- resolve.defaults(list(penalty    = penalty,
                                     penal.args = penal.args,
                                     tau        = tau,
                                     saveplace  = saveplace,
                                     TRACE      = TRACE),
                                objargs)
  } else {
    # create local function to evaluate  pair correlation(d) * weight(d)
    #  (with additional parameters 'pcfunargs', 'weightfun' in its environment)
    force(weightfun)
    wpaco <- function(d, par) {
      y <- do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
      w <- weightfun(d)
      return(y * w)
    }
    # pack up necessary information
    objargs <- list(dIJ=dIJ, wIJ=wIJ, sumweight=sumweight, g=g, gscale=gscale, 
                    envir=environment(wpaco),
                    penalty=NULL,   # updated below
                    penal.args=NULL,   # updated below
                    tau=NULL,   # updated below
                    TRACE=FALSE, # updated below
                    saveplace=NULL, # updated below
                    BIGVALUE=1, # updated below
                    SMALLVALUE=.Machine$double.eps)
    # define objective function (with 'paco', 'wpaco' in its environment)
    # This is the log composite likelihood minus the constant term 
    #       2 * (sum(wIJ * log(lambdaIJ)) - sumweight * log(gscale))
    obj <- function(par, objargs) {
      with(objargs,
      {
        integ <- unlist(stieltjes(wpaco, g, par=par))
        integ <- pmax(SMALLVALUE, integ)
        logcl <- safeFiniteValue(
          2*(sum(wIJ * log(safePositiveValue(paco(dIJ, par))))
            - sumweight * log(integ)),
          default=-BIGVALUE)
        ## penalty
        if(hasPenalty <- is.function(penalty)) {
          straf <- do.call(penalty, append(list(par), penal.args))
          logclPEN <- logcl - tau * straf
        }
        ## debugger
        if(isTRUE(TRACE)) {
          cat("Parameters:", fill=TRUE)
          print(par)
          splat("\tinteg:", integ)
          splat("log composite likelihood:", logcl)
          if(hasPenalty)
            splat("penalised log composite likelihood:", logclPEN)
        }
        if(is.environment(saveplace)) {
          h <- get("h", envir=saveplace)
          value <- list(logcl=logcl)
          if(hasPenalty)
            value <- append(value, list(logclPEN=logclPEN))
          hplus <- as.data.frame(append(par, value))
          h <- rbind(h, hplus)
          assign("h", h, envir=saveplace)
        }
        return(if(hasPenalty) logclPEN else logcl)
      },
      enclos=objargs$envir)
    }
    ## Determine the values of some parameters
    ## (1) Determine a suitable large number to replace Inf
    objargs$BIGVALUE <- bigvaluerule(obj, objargs, startpar)
    ## (2) Evaluate exchange rate 'tau'
    if(is.function(penalty) && is.function(tau)) {
      ## data-dependent exchange rate 'tau': evaluate now
      if("poisval" %in% names(formals(tau))) {
        ## New style: requires value of (unpenalised) objective function at Poisson process
        poisval <- obj(poispar, objargs)
        tau <- tau(X, poisval=poisval)
      } else {
        tau <- tau(X)
      }
      check.1.real(tau)
    }
    ## Now insert the penalty, etc
    objargs <- resolve.defaults(list(penalty    = penalty,
                                     penal.args = penal.args,
                                     tau        = tau,
                                     saveplace  = saveplace,
                                     TRACE      = TRACE),
                                objargs)
  }

  ## ......................  Optimization settings  ........................

  if(stabilize) {
    ## Numerical stabilisation 
    ## evaluate objective at starting state
    startval <- obj(startpar, objargs)
    ## use to determine appropriate global scale
    smallscale <- sqrt(.Machine$double.eps)
    fnscale <- -max(abs(startval), smallscale)
    parscale <- pmax(abs(startpar), smallscale)
    scaling <- list(fnscale=fnscale, parscale=parscale)
  } else {
    scaling <- list(fnscale=-1)
  }

  ## Update list of algorithm control arguments
  control.updated <- resolve.defaults(control, scaling, list(trace=0))

  ## Initialise list of all arguments to 'optim'
  optargs <- list(par=startpar, fn=obj, objargs=objargs,
                  control=control.updated, method=algorithm)

  ## DPP case: check startpar and modify algorithm
  changealgorithm <- length(startpar)==1 && algorithm=="Nelder-Mead"
  if(isDPP){
    alg <- dppmFixAlgorithm(algorithm, changealgorithm, clusters,
                            startpar.human)
    algorithm <- optargs$method <- alg$algorithm
    if(algorithm=="Brent" && changealgorithm){
      optargs$lower <- alg$lower
      optargs$upper <- alg$upper
    }
  }
  
  if(isTRUE(pspace$debug)) {
    splat("About to optimize... Objective function arguments:")
    print(objargs)
  }
  
  ## ..........   optimize it ..............................
  opt <- do.call(optim, optargs)

  ## raise warning/error if something went wrong
  signalStatus(optimStatus(opt), errors.only=TRUE)
  
  ## .......... extract fitted parameters .....................
  if(!usecanonical) {
    optpar.canon <- NULL
    optpar.human <- opt$par
    names(optpar.human) <- names(startpar.human)
  } else {
    optpar.canon <- opt$par
    names(optpar.canon) <- names(startpar)
    optpar.human <- tohuman(optpar.canon, margs=margs)
    names(optpar.human) <- names(startpar.human)
  }
  opt$par             <- optpar.human
  opt$par.canon       <- optpar.canon
  ## save starting values in 'opt' for consistency with mincontrast()
  opt$startpar        <- startpar.human

  ## Finish in DPP case
  if(!is.null(DPP)){
    ## all info that depends on the fitting method:
    Fit <- list(method       = "clik2",
                clfit        = opt,
                weightfun    = weightfun,
                rmax         = rmax,
                objfun       = obj,
                objargs      = objargs,
                maxlogcl     = opt$value,
                pspace.given = pspace.given,
                pspace.used  = pspace)
    # pack up
    clusters <- update(clusters, as.list(opt$par))
    result <- list(Xname      = Xname,
                   X          = X,
                   stationary = stationary,
                   fitted     = clusters,
                   modelname  = modelname,
                   po         = po,
                   lambda     = lambda,
                   Fit        = Fit)
    if(SAVE) {
      h <- get("h", envir=saveplace)
      if(!is.null(h)) class(h) <- unique(c("traj", class(h)))
      attr(result, "h") <- h
    }
    return(result)
  }
  
  ## meaningful model parameters
  modelpar <- info$interpret(optpar.human, lambda)
  # infer parameter 'mu'
  if(isPCP) {
    # Poisson cluster process: extract parent intensity kappa
    kappa <- optpar.human[["kappa"]]
    # mu = mean cluster size
    mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
  } else {
    # LGCP: extract variance parameter sigma2
    sigma2 <- optpar.human[["sigma2"]]
    # mu = mean of log intensity 
    mu <- if(stationary) log(lambda) - sigma2/2 else
          eval.im(log(lambda) - sigma2/2)    
  }
  # all info that depends on the fitting method:
  Fit <- list(method       = "clik2",
              clfit        = opt,
              weightfun    = weightfun,
              rmax         = rmax,
              objfun       = obj,
              objargs      = objargs,
              maxlogcl     = opt$value,
              pspace.given = pspace.given,
              pspace.used  = pspace)
  # pack up
  result <- list(Xname      = Xname,
                 X          = X,
                 stationary = stationary,
                 clusters   = clusters,
                 modelname  = modelname,
                 isPCP      = isPCP,
                 po         = po,
                 lambda     = lambda,
                 mu         = mu,
                 par        = optpar.human,
                 par.canon  = optpar.canon,
                 clustpar   = info$checkpar(par=optpar.human, native=FALSE, strict=strict),
                 clustargs  = info$outputshape(shapemodel$margs),
                 modelpar   = modelpar,
                 covmodel   = shapemodel,
                 Fit        = Fit)
  
  if(SAVE) {
    h <- get("h", envir=saveplace)
    if(!is.null(h)) class(h) <- unique(c("traj", class(h)))
    attr(result, "h") <- h
  }

  return(result)
}

