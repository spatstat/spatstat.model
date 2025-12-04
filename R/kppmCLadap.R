#'  kppmCLadap.R
#'
#'       Fitting algorithm for kppm(method = 'adapcl')
#'
#' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#'  A  d  a  p  t  i  v  e     C  o  m  p  o  s  i  t  e
#'               L  i  k  e  l  i  h  o  o  d
#' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#'
#'  Contributed by Chiara Fend 
#'  needs nonlinear equation solver nleqslv
#'
#'  $Revision: 1.3 $ $Date: 2025/12/04 05:32:38 $
#' 
#'  Copyright (c) 2020-2025 Chiara Fend, Ege Rubak,
#'                          Adrian Baddeley and Rolf Turner
#'  GNU Public Licence (>= 2.0)


kppmCLadap <- function(X, Xname, po, clusters, control, weightfun, 
                       rmax=NULL, epsilon=0.01, DPP=NULL,
                       algorithm="Broyden", ..., 
                       startpar=NULL, globStrat="dbldog") {

  if(!requireNamespace("nleqslv", quietly=TRUE))
    stop(paste("The package", sQuote("nleqslv"), "is required"),
           call.=FALSE)
  
  W <- as.owin(X)

  if(is.null(rmax)) # specified for numerical stability
    rmax <- shortside(Frame(W))

  # identify pairs of points that might contribute
  cl <- closepairs(X, rmax)
  dIJ <- cl$d #pairwise distances
  Rmin <- min(dIJ)
  indexmin <- which(dIJ==Rmin) #for later use
  
  # convert window to mask, saving other arguments for later
  dcm <- do.call.matched(as.mask,
                         append(list(w=W), list(...)),
                         sieve=TRUE)
  M         <- dcm$result
  otherargs <- dcm$otherargs

  #' penalised estimation is not supported
  if(any(m <- (c("penalised", "pspace") %in% names(otherargs)))) 
    warning(paste(ngettext(sum(m), "Argument", "Arguments"),
                  commasep(sQuote(c("penalised", "pspace")[m])),
                  ngettext(sum(m), "is", "are"),
                  "not supported for adaptive composite likelihood"),
            call.=FALSE)
    
  # compute intensity at pairs of data points
  # and c.d.f. of interpoint distance in window
  if(stationary <- is.stationary(po)) {
    # stationary unmarked Poisson process
    lambda <- intensity(X)
    # compute cdf of distance between two uniform random points in W
    g <- distcdf(W, delta=rmax/4096)
    # scaling constant is (area * intensity)^2
    gscale <- npoints(X)^2  
  } else {
    # compute fitted intensity at data points and in window
    #    lambdaX <- fitted(po, dataonly=TRUE)
    lambda <- lambdaM <- predict(po, locations=M)
    # compute cdf of distance between two random points in W
    # with density proportional to intensity function
    g <- distcdf(M, dW=lambdaM, delta=rmax/4096)
    # scaling constant is (integral of intensity)^2
    gscale <- safePositiveValue(integral.im(lambdaM)^2, default=npoints(X)^2)
  }
  
  isDPP <- !is.null(DPP)
  if(isDPP){
    tmp <- dppmFixIntensity(DPP, lambda, po)
    clusters <- tmp$clusters
    lambda <- tmp$lambda
    po <- tmp$po
  }
  
  # get pair correlation function (etc) for model
  info <- spatstatClusterModelInfo(clusters)
  pcfun        <- info$pcf
  dpcfun       <- info$Dpcf
  selfstart    <- info$selfstart
  isPCP        <- info$isPCP
  resolveshape <- info$resolveshape
  modelname    <- info$modelname
  # Assemble information required for computing pair correlation
  if(is.function(resolveshape)) {
    # Additional parameters of cluster model are required.
    # These may be given as individual arguments,
    # or in a list called 'covmodel'
    clustargs <- if("covmodel" %in% names(otherargs))
      otherargs[["covmodel"]] else otherargs
    shapemodel <- do.call(resolveshape, clustargs)$covmodel
  } else shapemodel <- NULL
  pcfunargs <- shapemodel
  
  ## determine starting parameter values
  if(is.null(startpar)) {
    startpar <- selfstart(X)
  } else if(!isDPP){
    startpar <- info$checkpar(startpar, native=TRUE)
  }
  ## optimization will be over the logarithms of the parameters
  startparLog <- log(startpar)
  pcfunLog <- function(par, ...) { pcfun(exp(par), ...) }
  dpcfunLog <- function(par, ...) { dpcfun(exp(par), ...) }
  
  # create local functions to evaluate pair correlation and its gradient
  #  (with additional parameters 'pcfunargs' in its environment)
  paco <- function(d, par) {
    do.call(pcfunLog, append(list(par=par, rvals=d), pcfunargs))
  }
  
  dpaco <- function(d, par) {
    do.call(dpcfunLog, append(list(par=par, rvals=d), pcfunargs))
  }
  
  # trim 'g' to [0, rmax] 
  g <- g[with(g, .x) <= rmax,]
  
  
  #' ..........  define objective function ......................
  # create local function to evaluate  weight(epsilon*M/(pcf(d)-1))
  weight <- function(d, par) {
    y <- paco(d=d, par=par)
    # calculate M (only needs to be calculated for cluster models)
    M <- 1
    if(!isDPP){
      M <- abs(paco(d=0, par=par)-1)
    }
    return(weightfun(epsilon*M/(y-1)))
  }

  wlogcl2score <- function(par, paco, dpaco, dIJ, gscale, epsilon, cdf=g){
    p <- length(par)
    temp <- rep(0, p)
    
    # check if current parameter is valid, if not return inf
    if(isDPP){
      if(length(par)==1 && is.null(names(par)))
        names(par) <- clusters$freepar
      mod <- update(clusters, as.list(exp(par)))
      if(!valid(mod)){
        return(rep(Inf, p))
      }
    }
    
    # everything can be computed
    wdIJ <- weight(d=dIJ, par=par)
    index <- unique(c(which(wdIJ!=0), indexmin))
    dIJcurrent <- dIJ[index]
    
    for(i in 1:p){
      parname <- names(par)[i]
      # weighted derivatives wrt log of parname
      dpcfweighted <- function(d, par){ 
        y <- dpaco(d = d, par = par)[parname,]*exp(par[i])
        return(y*weight(d = d, par = par))
      } 
      temp[i] <- sum(dpcfweighted(d = dIJcurrent, par=par)/paco(d = dIJcurrent, par = par)) - gscale * stieltjes(dpcfweighted,cdf, par=par)$f
    }
    return(temp)
  }
  
  ## .................   optimize it ..............................
  opt <- nleqslv::nleqslv(x = startparLog, fn = wlogcl2score, 
                          method = algorithm,
                          global = globStrat, control = control, 
                          paco=paco, dpaco=dpaco, 
                          dIJ=dIJ, gscale=gscale, epsilon=epsilon)
    
  ## .......... extract fitted parameters on original scale ...............
  optpar        <- exp(opt$x)
  names(optpar) <- names(startpar)
  ## insert entries expected in 'opt'
  opt$par      <- optpar
  opt$startpar <- startpar

  ## Finish in DPP case
  if(isDPP){
    # all info that depends on the fitting method:
    Fit <- list(method    = "adapcl",
                cladapfit = opt,
                weightfun = weightfun,
                rmax      = rmax,
                epsilon   = epsilon,
                objfun    = wlogcl2score,
                objargs   = control,
                estfunc  = opt$fvec)
    # pack up
    clusters <- update(clusters, as.list(exp(opt$x)))
    result <- list(Xname      = Xname,
                   X          = X,
                   stationary = stationary,
                   fitted     = clusters,
                   modelname  = modelname,
                   po         = po,
                   lambda     = lambda,
                   Fit        = Fit)
    return(result)
  }
  
  ## meaningful model parameters
  modelpar <- info$interpret(optpar, lambda)
  # infer parameter 'mu'
  if(isPCP) {
    # Poisson cluster process: extract parent intensity kappa
    kappa <- optpar[["kappa"]]
    # mu = mean cluster size
    mu <- if(stationary) lambda/kappa else eval.im(lambda/kappa)
  } else {
    # LGCP: extract variance parameter sigma2
    sigma2 <- optpar[["sigma2"]]
    # mu = mean of log intensity 
    mu <- if(stationary) log(lambda) - sigma2/2 else
      eval.im(log(lambda) - sigma2/2)    
  }
  # all info that depends on the fitting method:
  Fit <- list(method    = "adapcl",
              cladapfit = opt,
              weightfun = weightfun,
              rmax      = rmax,
              epsilon   = epsilon,
              objfun    = wlogcl2score,
              objargs   = control,
              estfunc  = opt$fvec)
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
                 par        = optpar,
                 clustpar   = info$checkpar(par=optpar, native=FALSE),
                 clustargs  = info$outputshape(shapemodel$margs),
                 modelpar   = modelpar,
                 covmodel   = shapemodel,
                 Fit        = Fit)
  return(result)
}

