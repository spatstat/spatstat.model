#'  kppmMinCon.R
#'
#'      Fitting algorithm for kppm(method = 'mincon'), the default
#' 
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#'                M  i  n  i  m  u  m       C  o  n  t  r  a  s  t
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#' 
#'  $Revision: 1.2 $ $Date: 2026/01/21 06:26:39 $
#'
#'  Copyright (c) 2001-2025 Adrian Baddeley, Rolf Turner, Ege Rubak,
#'  GNU Public Licence (>= 2.0)

kppmMinCon <- function(X, Xname, po, clusters, control=list(), stabilize=TRUE, statistic, statargs,
                       algorithm="Nelder-Mead", DPP=NULL, ...,
                       pspace=NULL) {
  # Minimum contrast fit
  stationary <- is.stationary(po)
  pspace <- do.call(make.pspace,
                    resolve.defaults(
                      list(fitmethod="mincon", clusters=clusters),
                      list(...), ## ellipsis arguments override pspace
                      as.list(pspace),
                      .MatchNull=FALSE))
  # compute intensity
  if(stationary) {
    lambda <- summary(po)$trend$value
  } else {
    # compute intensity at high resolution if available
    w <- as.owin(po, from="covariates")
    if(!is.mask(w)) w <- NULL
    lambda <- predict(po, locations=w)
  }
  # Detect DPP model and change clusters and intensity correspondingly
  if(!is.null(DPP)){
    tmp <- dppmFixIntensity(DPP, lambda, po)
    clusters <- tmp$clusters
    lambda <- tmp$lambda
    po <- tmp$po
  }
  mcfit <- clusterfit(X, clusters, lambda = lambda,
                      dataname = Xname, control = control,  stabilize=stabilize,
                      statistic = statistic, statargs = statargs,
                      algorithm=algorithm, pspace=pspace, ...)
  fitinfo <- attr(mcfit, "info")
  attr(mcfit, "info") <- NULL
  # all info that depends on the fitting method:
  Fit <- list(method       = "mincon",
              statistic    = statistic,
              Stat         = fitinfo$Stat,
              StatFun      = fitinfo$StatFun,
              StatName     = fitinfo$StatName,
              FitFun       = fitinfo$FitFun,
              statargs     = statargs,
              pspace.given = pspace,
              pspace.used  = fitinfo$pspace.used, 
              mcfit        = mcfit,
              maxlogcl     = NULL)
  # results
  if(!is.null(DPP)){
    clusters <- update(clusters, as.list(mcfit$par))
    out <- list(Xname      = Xname,
                X          = X,
                stationary = stationary,
                fitted     = clusters,
                po         = po,
                Fit        = Fit)
  } else{
    out <- list(Xname      = Xname,
                X          = X,
                stationary = stationary,
                clusters   = clusters,
                modelname  = fitinfo$modelname,
                isPCP      = fitinfo$isPCP,
                po         = po,
                lambda     = lambda,
                mu         = mcfit$mu,
                par        = mcfit$par,
                par.canon  = mcfit$par.canon,
                clustpar   = mcfit$clustpar,
                clustargs  = mcfit$clustargs,
                modelpar   = mcfit$modelpar,
                covmodel   = mcfit$covmodel,
                Fit        = Fit)
  }
  h <- attr(mcfit, "h")
  if(!is.null(h)) class(h) <- unique(c("traj", class(h)))
  attr(out, "h") <- h
  return(out)
}

clusterfit <- function(X, clusters, lambda = NULL, startpar = NULL,
                       ...,
                       q=1/4, p=2, rmin=NULL, rmax=NULL, 
                       ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
                       statistic = NULL, statargs = NULL,
                       algorithm="Nelder-Mead", verbose=FALSE,
                       pspace=NULL){
  if(verbose) splat("Fitting cluster model")
  ## If possible get dataname from dots
  dataname <- list(...)$dataname
  ## Cluster info:
  info <- spatstatClusterModelInfo(clusters)
  if(verbose) splat("Retrieved cluster model information")
  ## Determine model type
  isPCP <- isTRUE(info$isPCP)
  isDPP <- inherits(clusters, "detpointprocfamily")

  ## resolve algorithm parameters
  default.ctrl <- list(q=if(isDPP) 1/2 else 1/4,
                       p=2,
                       rmin=NULL,
                       rmax=NULL)
  given.ctrl <- if(missing(ctrl)) list() else ctrl[names(default.ctrl)]
  given.args <- c(if(missing(q)) NULL else list(q=q),
                  if(missing(p)) NULL else list(p=p),
                  if(missing(rmin)) NULL else list(rmin=rmin),
                  if(missing(rmax)) NULL else list(rmax=rmax))
  ctrl <- resolve.defaults(given.args, given.ctrl, default.ctrl)
  if(verbose) {
    splat("Algorithm parameters:")
    print(ctrl)
  }
  ##
  if(inherits(X, "ppp")){
      if(verbose) 
        splat("Using point pattern data")
      if(is.null(dataname))
         dataname <- getdataname(short.deparse(substitute(X), 20), ...)
      if(is.null(statistic))
          statistic <- "K"
      # Startpar:
      if(is.null(startpar))
          startpar <- info$selfstart(X)
      stationary <- is.null(lambda) || (is.numeric(lambda) && length(lambda)==1)
      if(verbose) {
        splat("Starting parameters:")
        print(startpar)
        cat("Calculating summary function...")
      }
      # compute summary function
      if(stationary) {
          if(is.null(lambda)) lambda <- intensity(X)
          StatFun <- if(statistic == "K") "Kest" else "pcf"
          StatName <-
              if(statistic == "K") "K-function" else "pair correlation function"
          Stat <- do.call(StatFun,
                          resolve.defaults(list(X=quote(X)),
                                           statargs,
                                           list(correction="best")))
      } else {
          StatFun <- if(statistic == "K") "Kinhom" else "pcfinhom"
          StatName <- if(statistic == "K") "inhomogeneous K-function" else
          "inhomogeneous pair correlation function"
          Stat <- do.call(StatFun,
                          resolve.defaults(list(X=quote(X), lambda=lambda),
                                           statargs,
                                           list(correction="best")))
      }
      if(verbose) splat("Done.")
  } else if(inherits(X, "fv")){
      if(verbose) 
        splat("Using the given summary function")
      Stat <- X
      ## Get statistic type
      stattype <- attr(Stat, "fname")
      StatFun <- paste0(stattype)
      StatName <- NULL
      if(is.null(statistic)){
          if(is.null(stattype) || !is.element(stattype[1L], c("K", "pcf")))
              stop("Cannot infer the type of summary statistic from argument ",
                   sQuote("X"), " please specify this via argument ",
                   sQuote("statistic"))
          statistic <- stattype[1L]
      }
      if(stattype[1L]!=statistic)
          stop("Statistic inferred from ", sQuote("X"),
               " not equal to supplied argument ",
               sQuote("statistic"))
      # Startpar:
      if(is.null(startpar)){
          if(isDPP)
              stop("No rule for starting parameters in this case. Please set ",
                   sQuote("startpar"), " explicitly.")
          startpar <- info$checkpar(startpar, native=FALSE)
          startpar[["scale"]] <- mean(range(Stat[[fvnames(Stat, ".x")]]))
      }
  } else{
      stop("Unrecognised format for argument X")
  }
  
  ## avoid using g(0) as it may be infinite
  if(statistic=="pcf"){
      if(verbose) splat("Checking g(0)")
      argu <- fvnames(Stat, ".x")
      rvals <- Stat[[argu]]
      if(rvals[1L] == 0 && (is.null(rmin) || rmin == 0)) {
          if(verbose) splat("Ignoring g(0)")
          rmin <- rvals[2L]
      }
  }

  ## DPP resolving algorithm and checking startpar
  changealgorithm <- length(startpar)==1 && algorithm=="Nelder-Mead"
  if(isDPP){
    if(verbose) splat("Invoking dppmFixAlgorithm")
    alg <- dppmFixAlgorithm(algorithm, changealgorithm, clusters, startpar)
    algorithm <- alg$algorithm
  }

  #' determine initial values of parameters
  startpar <- info$checkpar(startpar, native=TRUE)
  #' code to compute the theoretical summary function of the model
  theoret <- info[[statistic]]
  #' explanatory text
  desc <- paste("minimum contrast fit of", info$descname)

  #' determine shape parameters if any
  dots <- info$resolveshape(...)
  margs <- dots$margs
  
  #' ............ permit negative parameter values .........................
  strict <- !isFALSE(pspace$strict)
  sargs <- if(strict) list() else list(strict=FALSE)
  
  #' ............ adjust K function or pair correlation function ...........
  do.adjust <- isTRUE(pspace$adjusted)
  if(do.adjust) {
    if(verbose) splat("Applying kppm adjustment")
    W <- Window(X)
    delta <- if(is.null(rmax)) NULL else rmax/4096
    ## pack up precomputed information needed for adjustment
    adjdata <- list(paircorr = info[["pcf"]],
                    pairWcdf = distcdf(W, delta=delta),
		    tohuman  = NULL)
    adjfun <- function(theo, par, auxdata, ..., margs=NULL) {
      with(auxdata, {
        if(!is.null(tohuman))
	  par <- tohuman(par, ..., margs=margs)
        a <- as.numeric(stieltjes(paircorr, pairWcdf, par=par, ..., margs=margs))
	return(theo/a)
      })
    }
    pspace$adjustment <- list(fun=adjfun, auxdata=adjdata)
  } 

  #' parameter vector corresponding to Poisson process
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
    if(verbose) splat("Converting to canonical parameters")
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
    htheo <- theoret
    startpar <- tocanonical(startpar, margs=margs)
    if(!is.null(poispar)) poispar <- tocanonical(poispar, margs=margs)
    theoret <- function(par, ...) { htheo(tohuman(par, ...), ...) }
    if(do.adjust)
      pspace$adjustment$auxdata$tohuman <- tohuman
  }
  #' ............ penalty .........................
  penalty    <- pspace$penalty
  penal.args <- pspace$penal.args
  tau        <- pspace$tau %orifnull% 1
  if(is.function(penalty)) {
    # penalised optimisation
    if(usecanonical) {
      penalty.human <- penalty
      penalty <- function(par, ...) { penalty.human(tohuman(par, ...), ...) }
    }
    ## data-dependent arguments in penalty
    if(is.function(penal.args)) 
      penal.args <- penal.args(X)
    ## exchange rate (defer evaluation if it is a function)
    if(!is.function(tau)) check.1.real(tau)
    ## reinsert in 'pspace' to pass to 'mincontrast'
    pspace$penalty <- penalty
    pspace$penal.args <- penal.args
    pspace$tau <- tau
    ## unpenalised version
    pspace.unpen <- pspace
    pspace.unpen[c("penalty", "penal.args", "tau")] <- NULL
  }
  #' ...................................................

  #'
  mcargs <- resolve.defaults(list(observed=Stat,
                                  theoretical=theoret,
                                  startpar=startpar,
                                  ctrl=ctrl,
                                  method=algorithm,
                                  fvlab=list(label="%s[fit](r)",
                                      desc=desc),
                                  explain=list(dataname=dataname,
                                      fname=statistic,
                                      modelname=info$modelname),
                                  margs=dots$margs,
                                  model=dots$model,
                                  pspace=pspace), ## As modified above
                             list(...),
                             sargs)


  if(isDPP && algorithm=="Brent" && changealgorithm)
    mcargs <- resolve.defaults(mcargs, list(lower=alg$lower, upper=alg$upper))

  if(is.function(penalty) && is.function(tau)) {
    ## data-dependent exchange rate 'tau': evaluate now
    if("poisval" %in% names(formals(tau))) {
      ## New style: requires value of (unpenalised) objective function at Poisson process
      mcargs.unpen <- mcargs
      mcargs.unpen$pspace <- pspace.unpen
      ## Evaluate using undocumented argument 'evalpar' to mincontrast
      mcargs.unpen$evalpar <- poispar
      poisval <- do.call(mincontrast, mcargs.unpen)
      tau <- tau(X, poisval=poisval)
    } else {
      ## old style
      tau <- tau(X)
    }
    check.1.real(tau)
    ## update 'tau' in argument list
    pspace$tau <- tau
    mcargs$pspace$tau <- tau
  }
  
  ## .............. FIT .......................
  if(verbose) splat("Starting minimum contrast fit")
  mcfit <- do.call(mincontrast, mcargs)
  if(verbose) splat("Returned from minimum contrast fit")
  ## ..........................................

  ## extract fitted parameters and reshape
  if(!usecanonical) {
    optpar.canon <- NULL
    optpar.human <- mcfit$par
    names(optpar.human) <- names(startpar.human)
  } else {
    optpar.canon <- mcfit$par
    names(optpar.canon) <- names(startpar)
    optpar.human <- tohuman(optpar.canon, margs=margs)
    names(optpar.human) <- names(startpar.human)
  }
  mcfit$par       <- optpar.human
  mcfit$par.canon <- optpar.canon
  ## extract fitted parameters
  optpar        <- mcfit$par
  names(optpar) <- names(startpar)
  mcfit$par     <- optpar
  
  # Return results for DPPs
  if(isDPP){
    extra <- list(Stat      = Stat,
                  StatFun   = StatFun,
                  StatName  = StatName,
                  modelname  = info$modelabbrev,
                  lambda     = lambda)
    attr(mcfit, "info") <- extra
    if(verbose) splat("Returning from clusterfit (DPP case)")
    return(mcfit)
  }
  ## Extra stuff for ordinary cluster/lgcp models
  ## imbue with meaning
  ## infer model parameters
  mcfit$modelpar <- info$interpret(optpar.human, lambda)
  mcfit$internal <- list(model=ifelse(isPCP, clusters, "lgcp"))
  mcfit$covmodel <- dots$covmodel
  
  if(isPCP) {
    # Poisson cluster process: extract parent intensity kappa
    kappa <- mcfit$par[["kappa"]]
    # mu = mean cluster size
    mu <- lambda/kappa
  } else {
    # LGCP: extract variance parameter sigma2
    sigma2 <- mcfit$par[["sigma2"]]
    # mu = mean of log intensity 
    mu <- log(lambda) - sigma2/2
  }
  ## Parameter values (new format)
  mcfit$mu <- mu
  mcfit$clustpar <- info$checkpar(mcfit$par, native=FALSE, strict=strict)
  mcfit$clustargs <- info$outputshape(dots$margs)

  ## The old fit fun that would have been used (DO WE NEED THIS?)
  FitFun <- paste0(tolower(clusters), ".est", statistic)

  extra <- list(FitFun       = FitFun,
                Stat         = Stat,
                StatFun      = StatFun,
                StatName     = StatName,
                modelname    = info$modelabbrev,
                isPCP        = isPCP,
                lambda       = lambda,
                pspace.used  = pspace) # Modified from call to 'clusterfit'
  attr(mcfit, "info") <- extra
  if(verbose) splat("Returning from clusterfit")
  return(mcfit)
}

