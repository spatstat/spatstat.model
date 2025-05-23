#'
#'    simulatekppm.R
#'
#'    simulate.kppm
#'
#'    $Revision: 1.14 $ $Date: 2025/05/22 05:25:49 $

simulate.kppm <- function(object, nsim=1, seed=NULL, ...,
                          window=NULL, covariates=NULL,
                          n.cond=NULL, w.cond=NULL,
                          verbose=TRUE, retry=10,
                          drop=FALSE) {
  starttime <- proc.time()
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  if(nsim == 0) return(simulationresult(list()))
  verbose <- verbose && (nsim > 1)
  check.1.real(retry)
  # .... copied from simulate.lm ....
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  ## ..................................
  ## determine window for simulation results
  if(!is.null(window)) {
    stopifnot(is.owin(window))
    win <- window
  } else {
    win <- as.owin(object)
  }
  ## ..................................
  ## conditional simulation 
  if(!is.null(n.cond)) {
    ## fixed number of points
    out <- CondSimCox(object, nsim=nsim, seed=NULL, ..., 
                      win=win, covariates=covariates, 
                      n.cond=n.cond, w.cond=w.cond,
                      verbose=verbose, retry=retry, drop=drop)
    out <- timed(out, starttime=starttime)
    attr(out, "seed") <- RNGstate
    return(out)
  }

  ## ..................................
  # determine parameters
  mp <- as.list(object$modelpar)

  # parameter 'mu'
  # = parent intensity of cluster process
  # = mean log intensity of log-Gaussian Cox process
  
  if(is.null(covariates) && (object$stationary || is.null(window))) {
    # use existing 'mu' (scalar or image)
    mu <- object$mu
  } else {
    # recompute 'mu' using new data
    switch(object$clusters,
           Cauchy=,
           VarGamma=,
           Thomas=,
           MatClust={
             # Poisson cluster process
             kappa <- mp$kappa
             lambda <- predict(object, window=win, covariates=covariates)
             mu <- eval.im(lambda/kappa)
           },
           LGCP={
             # log-Gaussian Cox process
             sigma2 <- mp$sigma2
             lambda <- predict(object, window=win, covariates=covariates)
             mu <- eval.im(log(lambda) - sigma2/2)
           },
           stop(paste("Simulation of", sQuote(object$clusters),
                      "processes is not yet implemented"))
           )
  }

  ## simulate
  switch(object$clusters,
         Thomas={
           kappa <- mp$kappa
           sigma <- mp$sigma
           out <- rThomas(kappa=kappa, scale=sigma, mu=mu, win=win,
                          ..., nsim=nsim, drop=FALSE)
         },
         MatClust={
           kappa <- mp$kappa
           r     <- mp$R
           out   <- rMatClust(kappa=kappa, scale=r, mu=mu, win=win,
                              ..., nsim=nsim, drop=FALSE)
         },
         Cauchy = {
           kappa <- mp$kappa
           omega <- mp$omega
           out   <- rCauchy(kappa = kappa, scale=omega, mu=mu, win=win,
                            ..., nsim=nsim, drop=FALSE)
         },
         VarGamma = {
           kappa  <- mp$kappa
           omega  <- mp$omega
           nu.ker <- object$covmodel$margs$nu.ker
           out    <- rVarGamma(kappa=kappa, scale=omega, mu=mu,
                               nu=nu.ker, win=win,
                               ..., nsim=nsim, drop=FALSE)
         },
         LGCP={
           sigma2 <- mp$sigma2
           alpha  <- mp$alpha
           cm <- object$covmodel
           model <- cm$model
           margs <- cm$margs
           param <- append(list(var=sigma2, scale=alpha), margs)
           #'
           if(!is.im(mu)) {
             #' model will be simulated in 'win'
             out <- rLGCP(model=model, mu=mu, param=param,
                          win=win,
                          ..., nsim=nsim, drop=FALSE)
           } else {
             #' model will be simulated in as.owin(mu), then change window
             outwin <- rLGCP(model=model, mu=mu, param=param,
                             ..., nsim=nsim, drop=FALSE)
             #' clip to 'win'
             out <- solapply(outwin, "[", i=win)
             #' reinstate Lambda images
             if(isTRUE(list(...)$saveLambda)) {
               for(j in seq_along(out))
                 attr(out[[j]], "Lambda") <- attr(outwin[[j]], "Lambda")
             }
           }
         })
  #' pack up
  out <- simulationresult(out, nsim, drop)
  out <- timed(out, starttime=starttime)
  attr(out, "seed") <- RNGstate
  return(out)
}

