#'
#'   kppmWaagLik.R
#'
#'       Fitting algorithm for kppm(method = 'waag')
#' 
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#'               W  a  a  g  e  p  e  t  e  r  s  e  n
#'                 L  i  k  e  l  i  h  o  o  d              
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#'
#' Implementation of Waagepetersen (2007) second order composite likelihood
#'
#' Waagepetersen, R.P. (2007)
#' An estimating function approach to inference for inhomogeneous
#' Neyman-Scott processes. Biometrics 63, 252-258.
#'
#'  $Revision: 1.2 $ $Date: 2025/12/05 02:01:24 $
#'
#'  Copyright (c) 2001-2025 Adrian Baddeley, Rolf Turner, Ege Rubak
#'                and Bethany McDonald
#'             
#'  GNU Public Licence (>= 2.0)

kppmWaagLik <- function(X, Xname, po, clusters, control=list(),
                        stabilize=TRUE, weightfun, rmax,
                        algorithm="Nelder-Mead", DPP=NULL, ...,
                        pspace=NULL) {
  pspace <- do.call(make.pspace,
                    resolve.defaults(
                      list(fitmethod="waag", clusters=clusters),
                      list(...), ## ellipsis arguments override pspace
                      as.list(pspace),
                      .MatchNull=FALSE))
  W <- as.owin(X)
  if(is.null(rmax))
    rmax <- rmax.rule("K", W, intensity(X))
  
  #' identify unordered pairs of points { i, j} that contribute
  #' listing each pair only once, with i < j
  cl <- closepairs(X, rmax, twice=FALSE, neat=FALSE)
  dIJ <- cl$d
  npairs <- length(dIJ)
  #' compute weights for unordered pairs of points. Must be symmetric.
  wIJ <- if(is.function(weightfun)) weightfun(dIJ) else rep.int(1, npairs)
  
  #' convert window to mask, saving other arguments for later
  dcm <- do.call.matched(as.mask,
                         append(list(w=W), list(...)),
                         sieve=TRUE)
  M         <- dcm$result
  otherargs <- dcm$otherargs

  #' Detect DPP usage
  isDPP <- inherits(clusters, "detpointprocfamily")

  #' compute intensity at data points, intensity in window,
  #' and c.d.f. of interpoint distance in window (for Borel's method)
  if(stationary <- is.stationary(po)) {
    #' constant intensity
    lambda <- intensity(X)
    #' compute cdf of distance between two uniformly random points in W
    PairDistCDF <- distcdf(M, delta=rmax/4096)
    #' scaling constant is (intensity * area)^2
    PairDistCDFScale <- (lambda * area(M))^2
  } else {
    #' spatially-varying intensity -- compute at data points and in window
    lambdaX <- fitted(po, dataonly=TRUE)
    lambda <- lambdaM <- predict(po, locations=M)
    #' extract intensity at close pairs (listed only with i < j)
    lambdaI <- lambdaX[cl$i] 
    lambdaJ <- lambdaX[cl$j] 
    #' compute cdf of distance between two random points in W,
    #' i.i.d. with density proportional to intensity function
    PairDistCDF <- distcdf(M, dW=lambdaM, delta=rmax/4096)
   #' scaling constant is (integral of intensity)^2
    PairDistCDFScale <- safePositiveValue(integral.im(lambdaM)^2,
                                    default=npoints(X)^2)
  }

  #' Detect DPP model and change clusters and intensity correspondingly
  isDPP <- !is.null(DPP)
  if(isDPP){
    tmp <- dppmFixIntensity(DPP, lambda, po)
    clusters <- tmp$clusters
    lambda <- tmp$lambda
    po <- tmp$po
  }

  #' trim 'PairDistCDF' to [0, rmax] 
  PairDistCDF <- PairDistCDF[with(PairDistCDF, .x) <= rmax,]
  #' get pair correlation function (etc) for model
  info <- spatstatClusterModelInfo(clusters)
  pcfun        <- info$pcf
  selfstart    <- info$selfstart
  isPCP        <- info$isPCP
  resolveshape <- info$resolveshape
  modelname    <- info$modelname
  #' Assemble information required for computing pair correlation
  if(is.function(resolveshape)) {
    #' Additional 'shape' parameters of cluster model are required.
    #' These may be given as individual arguments,
    #' or in a list called 'covmodel'
    clustargs <- if("covmodel" %in% names(otherargs))
                 otherargs[["covmodel"]] else otherargs
    shapemodel <- do.call(resolveshape, clustargs)$covmodel
  } else shapemodel <- NULL
  pcfunargs <- shapemodel
  margs <- pcfunargs$margs
  #' determine starting parameter values
  startpar <- selfstart(X)
  
  pspace.given <- pspace
  #' ............ permit negative parameter values ................
  strict <- !isFALSE(pspace$strict)
  if(!strict) pcfunargs <- append(pcfunargs, list(strict=strict))
  #' ............ parameter corresponding to Poisson process ......
  if(isDPP) {
    poispar <- NULL
  } else if(isPCP) {
    if(!("kappa" %in% names(startpar)))
      stop("Internal error: startpar does not include 'kappa'")
    poispar <- startpar
    poispar[["kappa"]] <- Inf
  } else {
    #' LGCP
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
  #' ............ penalty .......................................
  penalty <- pspace$penalty
  penal.args <- pspace$penal.args
  tau <- pspace$tau %orifnull% 1
  if(is.function(penalty)) {
    #' penalised optimisation
    if(usecanonical) {
      penalty.human <- penalty
      penalty <- function(par, ...) { penalty.human(tohuman(par, ...), ...) }
    }
    #' data-dependent arguments in penalty
    if(is.function(penal.args)) 
      penal.args <- penal.args(X)
    #' exchange rate (defer evaluation if it is a function)
    if(!is.function(tau)) check.1.real(tau)
    #' reinsert in 'pspace' for insurance
    pspace$penalty <- penalty
    pspace$penal.args <- penal.args
    pspace$tau <- tau
  }
  #' ............ debugger ......................................
  TRACE <- isTRUE(pspace$trace)
  if(SAVE <- isTRUE(pspace$save)) {
    saveplace <- new.env()
    assign("h", NULL, envir=saveplace)
  } else saveplace <- NULL

  #' >>>>>>>>>>>>>>>>> CREATE OBJECTIVE FUNCTION <<<<<<<<<<<<<<<<<<<<<<<
  #'
  #' create local function to evaluate pair correlation
  #'  (with additional parameters 'pcfunargs' in its environment)
  paco <- function(d, par) {
    do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
  }
  #'
  if(!is.function(weightfun)) {
    #' Unweighted case -----------------------------------------------------
    #' First calculate the term which is not dependent on cluster parameters
    #' \sum\sum_{i \neq j} C(x_i, x_j) log(lambda(x_i) * lambda(x_j))
    #' where C(u, v) = 1 iff ||u - v|| <= rmax.
    #' Factor 2 is present below because close pairs were listed only once
    #' (with I < J)
    sumloglamlam <- 2 * 
      if(stationary) {
        npairs * safeFiniteValue(2 * log(lambda))
      } else {
        safeFiniteValue(sum(log(safePositiveValue(lambdaI * lambdaJ))))
      }
    #' pack up necessary information for objective function
    objargs <- list(dIJ              = dIJ,
                    PairDistCDF      = PairDistCDF,
                    PairDistCDFScale = PairDistCDFScale,
                    sumloglamlam     = sumloglamlam,
                    envir            = environment(paco),
                    #' The following variables are temporarily omitted
                    #' in order to calculate the objective function
                    #' without using them, or their side effects.
                    penalty          = NULL,   # updated below
                    penal.args       = NULL,   # updated below
                    tau              = NULL,   # updated below
                    TRACE            = FALSE, # updated below
                    saveplace        = NULL, # updated below
                    BIGVALUE         = 1, # updated below
                    SMALLVALUE       = .Machine$double.eps)
    #' Define OBJECTIVE FUNCTION (with 'paco' in its environment)
    #' This evaluates the log Waagepetersen composite likelihood
    obj <- function(par, objargs) {
      with(objargs, {
        #' calculate \sum\sum_{i \neq j} C(x_i, x_j) log pcf(x_i, x_j)
        #' Factor 2 is present because close pairs were listed only once in dIJ
        sumlogpaco <- 2 * sum(log(safePositiveValue(paco(dIJ, par)))) 
        #' calculate \int_W \int_W pcf(u - v) du dv using Borel's method
        integ <- unlist(stieltjes(paco, PairDistCDF, par=par))
        integ <- PairDistCDFScale * pmax(SMALLVALUE, integ)
        #' evaluate surrogate likelihood
        logwlik <- safeFiniteValue(
          sumloglamlam + sumlogpaco - integ,
          default=-BIGVALUE)
        #' penalty
        if(hasPenalty <- is.function(penalty)) {
          straf <- do.call(penalty, append(list(par), penal.args))
          logwlikPEN <- logwlik - tau * straf
        }
        #' debugger
        if(isTRUE(TRACE)) {
          cat("Parameters:", fill=TRUE)
          print(par)
          splat("integral:", integ)
          splat("log Waagepetersen likelihood:", logwlik)
          if(hasPenalty)
            splat("penalised log Waagepetersen likelihood:", logwlikPEN)
        }
        if(is.environment(saveplace)) {
          h <- get("h", envir=saveplace)
          value <- list(logwlik=logwlik)
          if(hasPenalty)
            value <- append(value, list(logwlikPEN=logwlikPEN))
          hplus <- as.data.frame(append(par, value))
          h <- rbind(h, hplus)
          assign("h", h, envir=saveplace)
        }
        return(if(hasPenalty) logwlikPEN else logwlik)
      },
      enclos=objargs$envir)
    }
    #' Determine the values of some parameters
    #' (1) Determine a suitable large number to replace Inf
    objargs$BIGVALUE <- bigvaluerule(obj, objargs, startpar)
    #' (2) Evaluate exchange rate 'tau'
    if(is.function(penalty) && is.function(tau)) {
      #' data-dependent exchange rate 'tau': evaluate now
      if("poisval" %in% names(formals(tau))) {
        #' New style: requires value of (unpenalised) objective function at Poisson process
        poisval <- obj(poispar, objargs)
        tau <- tau(X, poisval=poisval)
      } else {
        tau <- tau(X)
      }
      check.1.real(tau)
    }
    objargs <- resolve.defaults(list(penalty    = penalty,
                                     penal.args = penal.args,
                                     tau        = tau,
                                     saveplace  = saveplace,
                                     TRACE      = TRACE),
                                objargs)
  } else {
    #' -------------- Weighted case  -----------------------------------------
    #' Create local function to evaluate  pair correlation(d) * weight(d)
    #'  (with additional parameters 'pcfunargs', 'weightfun' in its environment)
    force(weightfun)
    wpaco <- function(d, par) {
      y <- do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
      w <- weightfun(d)
      return(y * w)
    }
    #' Calculate conserved quantity
    #' \sum \sum_{i \neq j} C(x_i, x_j) w(x_i, x_j) log(lambda(x_i) * lambda(x_j))
    #' where C(u, v) = 1 iff ||u - v|| <= rmax.
    #' Factor 2 is present below because close pairs were listed only once (with I < J)
    sumwloglamlam <- 2 *
      if(stationary) {
        safeFiniteValue(sum(wIJ) * 2 * log(safePositiveValue(lambda)))
      } else {
        safeFiniteValue(sum(wIJ * log(safePositiveValue(lambdaI * lambdaJ))))
      }
    #' pack up information necessary for objective function
    objargs <- list(dIJ              = dIJ,
                    wIJ              = wIJ,
                    PairDistCDF      = PairDistCDF,
                    PairDistCDFScale = PairDistCDFScale,
                    sumwloglamlam    = sumwloglamlam,
                    envir            = environment(wpaco),
                    penalty          = NULL,   # updated below
                    penal.args       = NULL,   # updated below
                    tau              = NULL,     # updated below
                    TRACE            = FALSE, # updated below
                    saveplace        = NULL, # updated below
                    BIGVALUE         = 1, # updated below
                    SMALLVALUE       = .Machine$double.eps)
    #' define objective function (with 'paco', 'wpaco' in its environment)
    #' This is the log Waagepetersen likelihood
    obj <- function(par, objargs) {
      with(objargs, {
        #' calculate \sum\sum_{i \neq j} w(x_i, x_j) log pcf(x_i, x_j)
        #' Factor 2 is present below because close pairs
        #' were listed only once (with I < J)
        sumwlogpaco <- 2 * safeFiniteValue(
                       sum(wIJ * log(safePositiveValue(paco(dIJ, par)))))
        #' calculate \int_W \int_W w(u,v) pcf(u, v) \lambda(u) \lambda(v) du dv
        #' using Borel's method        
        integ <- unlist(stieltjes(wpaco, PairDistCDF, par=par))
        integ <- PairDistCDFScale * pmax(SMALLVALUE, integ)
        #' evaluate surrogate likelihood
        logwlik <- safeFiniteValue(sumwloglamlam + sumwlogpaco - integ,
                                   default=-BIGVALUE)
        #' penalty
        if(hasPenalty <- is.function(penalty)) {
          straf <- do.call(penalty, append(list(par), penal.args))
          logwlikPEN <- logwlik - tau * straf
        }
        #' debugger
        if(isTRUE(TRACE)) {
          cat("Parameters:", fill=TRUE)
          print(par)
          splat("integral:", integ)
          splat("log Waagepetersen likelihood:", logwlik)
          if(hasPenalty)
            splat("penalised log Waagepetersen likelihood:", logwlikPEN)
        }
        if(is.environment(saveplace)) {
          h <- get("h", envir=saveplace)
          value <- list(logwlik=logwlik)
          if(hasPenalty)
            value <- append(value, list(logwlikPEN=logwlikPEN))
          hplus <- as.data.frame(append(par, value))
          h <- rbind(h, hplus)
          assign("h", h, envir=saveplace)
        }
        return(if(hasPenalty) logwlikPEN else logwlik)        
      },
      enclos=objargs$envir)
    }
    #' Determine the values of some parameters
    #' (1) Determine a suitable large number to replace Inf
    objargs$BIGVALUE <- bigvaluerule(obj, objargs, startpar)
    #' (2) Evaluate exchange rate 'tau'
    if(is.function(penalty) && is.function(tau)) {
      #' data-dependent exchange rate 'tau': evaluate now
      if("poisval" %in% names(formals(tau))) {
        #' New style: requires value of (unpenalised) objective function at Poisson process
        poisval <- obj(poispar, objargs)
        tau <- tau(X, poisval=poisval)
      } else {
        tau <- tau(X)
      }
      check.1.real(tau)
    }
    #' Now insert penalty, etc.
    objargs <- resolve.defaults(list(penalty    = penalty,
                                     penal.args = penal.args,
                                     tau        = tau,
                                     saveplace  = saveplace,
                                     TRACE      = TRACE),
                                objargs)
  }

  #' ......................  Optimization settings  ........................

  if(stabilize) {
    #' Numerical stabilisation 
    #' evaluate objective at starting state
    startval <- obj(startpar, objargs)
    #' use to determine appropriate global scale
    smallscale <- sqrt(.Machine$double.eps)
    fnscale <- -max(abs(startval), smallscale)
    parscale <- pmax(abs(startpar), smallscale)
    scaling <- list(fnscale=fnscale, parscale=parscale)
  } else {
    scaling <- list(fnscale=-1)
  }

  #' Update list of algorithm control arguments
  control.updated <- resolve.defaults(control, scaling, list(trace=0))

  #' Initialise list of all arguments to 'optim'
  optargs <- list(par=startpar, fn=obj, objargs=objargs,
                  control=control.updated, method=algorithm)

  #' DPP case: check startpar and modify algorithm
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

  #' .......................................................................
  
  if(isTRUE(pspace$debug)) {
    splat("About to optimize... Objective function arguments:")
    print(objargs)
  }

  #' optimize it
  opt <- do.call(optim, optargs)
  #' raise warning/error if something went wrong
  signalStatus(optimStatus(opt), errors.only=TRUE)
  
  #' Extract optimal values of parameters
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
  opt$par            <- optpar.human
  opt$par.canon      <- optpar.canon
  #' save starting values in 'opt' for consistency with minconfit()
  opt$startpar       <- startpar.human

  #' Finish in DPP case
  if(!is.null(DPP)){
    #' all info that depends on the fitting method:
    Fit <- list(method       = "waag",
                clfit        = opt,
                weightfun    = weightfun,
                rmax         = rmax,
                objfun       = obj,
                objargs      = objargs,
                maxlogcl     = opt$value,
                pspace.given = pspace.given,
                pspace.used  = pspace)
    #' pack up
    clusters <- update(clusters, as.list(optpar.human))
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
  #' meaningful model parameters
  modelpar <- info$interpret(optpar.human, lambda)
  #' infer parameter 'mu'
  if(isPCP) {
    #' Poisson cluster process: extract parent intensity kappa
    kappa <- optpar.human[["kappa"]]
    #' mu = mean cluster size
    mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
  } else {
    #' LGCP: extract variance parameter sigma2
    sigma2 <- optpar.human[["sigma2"]]
    #' mu = mean of log intensity 
    mu <- if(stationary) log(lambda) - sigma2/2 else
          eval.im(log(lambda) - sigma2/2)    
  }
  #' all info that depends on the fitting method:
  Fit <- list(method       = "waag",
              clfit        = opt,
              weightfun    = weightfun,
              rmax         = rmax,
              objfun       = obj,
              objargs      = objargs,
              maxlogcl     = opt$value,
              pspace.given = pspace.given,
              pspace.used  = pspace)
  #' pack up
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

