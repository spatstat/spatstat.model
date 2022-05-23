#'
#'    simulatekppm.R
#'
#'    simulate.kppm
#'
#'    $Revision: 1.9 $ $Date: 2022/04/06 08:51:41 $

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
    out <- condSimCox(object, nsim=nsim, seed=NULL, ..., 
                      window=win, covariates=covariates, 
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

  # prepare data for execution
  out <- list()
  switch(object$clusters,
         Thomas={
           kappa <- mp$kappa
           sigma <- mp$sigma
           cmd <- expression(rThomas(kappa,sigma,mu,win, ...))
           dont.complain.about(kappa, sigma, mu)
         },
         MatClust={
           kappa <- mp$kappa
           r     <- mp$R
           cmd   <- expression(rMatClust(kappa,r,mu,win, ...))
           dont.complain.about(kappa, r)
         },
         Cauchy = {
           kappa <- mp$kappa
           omega <- mp$omega
           cmd   <- expression(rCauchy(kappa, omega, mu, win, ...))
           dont.complain.about(kappa, omega, mu)
         },
         VarGamma = {
           kappa  <- mp$kappa
           omega  <- mp$omega
           nu.ker <- object$covmodel$margs$nu.ker
           cmd    <- expression(rVarGamma(kappa, nu.ker, omega, mu, win, ...))
           dont.complain.about(kappa, nu.ker, omega, mu)
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
             # model will be simulated in 'win'
             cmd <- expression(rLGCP(model=model, mu=mu, param=param,
                               ..., win=win))
             #' check that RandomFields package recognises parameter format
             rfmod <- try(rLGCP(model, mu=mu, param=param, win=win,
                              ..., modelonly=TRUE))
           } else {
             # model will be simulated in as.owin(mu), then change window
             cmd <- expression(rLGCP(model=model, mu=mu, param=param,
                               ...)[win])
             #' check that RandomFields package recognises parameter format
             rfmod <- try(rLGCP(model, mu=mu, param=param, 
                              ..., modelonly=TRUE))
           }
           #' suppress warnings from code checker
           dont.complain.about(model, mu, param)
           #' check that model is recognised
           if(inherits(rfmod, "try-error"))
             stop(paste("Internal error in simulate.kppm:",
                        "unable to build Random Fields model",
                        "for log-Gaussian Cox process"))
         })
  
  # run
  if(verbose) {
    cat(paste("Generating", nsim, "simulations... "))
    state <- list()
  }
  for(i in 1:nsim) {
    out[[i]] <- try(eval(cmd))
    if(verbose) state <- progressreport(i, nsim, state=state)
  }
  # detect failures
  if(any(bad <- unlist(lapply(out, inherits, what="try-error")))) {
    nbad <- sum(bad)
    gripe <- paste(nbad,
                   ngettext(nbad, "simulation was", "simulations were"),
                   "unsuccessful")
    if(verbose) splat(gripe)
    if(retry <= 0) {
      fate <- "returned as NULL"
      out[bad] <- list(NULL)
    } else {
      if(verbose) cat("Retrying...")
      ntried <- 0
      while(ntried < retry) {
        ntried <- ntried + 1
        for(j in which(bad))
          out[[j]] <- try(eval(cmd))
        bad <- unlist(lapply(out, inherits, what="try-error"))
        nbad <- sum(bad)
        if(nbad == 0) break
      }
      if(verbose) cat("Done.\n")
      fate <- if(nbad == 0) "all recomputed" else
              paste(nbad, "simulations still unsuccessful")
      fate <- paste(fate, "after", ntried,
                    ngettext(ntried, "further try", "further tries"))
    }
    warning(paste(gripe, fate, sep=": "))
  }
  if(verbose)
    cat("Done.\n")
  #' pack up
  out <- simulationresult(out, nsim, drop)
  out <- timed(out, starttime=starttime)
  attr(out, "seed") <- RNGstate
  return(out)
}

condSimCox <- function(object, nsim=1,
                       ..., window=NULL,
                       n.cond=NULL, w.cond=NULL,
                       giveup=1000, maxchunk=100,
                       verbose=TRUE, drop=FALSE) {
  stopifnot(is.kppm(object))
  shortcut <- isFALSE(object$isPCP)

  w.sim <- as.owin(window)
  fullwindow <- is.null(w.cond)
  if(fullwindow) {
    w.cond <- w.sim
    w.free <- NULL
  } else {
    stopifnot(is.owin(w.cond))
    w.free <- setminus.owin(w.sim, w.cond)
  }
  
  nremaining <- nsim
  ntried <- 0
  accept <- FALSE
  nchunk <- 1
  phistory <- mhistory <- numeric(0)
  results <- list()
  while(nremaining > 0) {
    ## increase chunk length
    nchunk <- min(maxchunk, giveup - ntried, 2 * nchunk)
    ## bite off next chunk of simulations
    if(shortcut) {
      lamlist <- simulate(object, nsim=nchunk,
                          Lambdaonly=TRUE,
                          ..., drop=FALSE, verbose=FALSE)
    } else {
      Xlist <- simulate(object, nsim=nchunk,
                        saveLambda=TRUE,
                        ..., drop=FALSE, verbose=FALSE)
      lamlist <- lapply(unname(Xlist), attr, which="Lambda", exact=TRUE)
    }
    ## compute acceptance probabilities
    lamlist <- lapply(lamlist, "[", i=w.sim, drop=FALSE, tight=TRUE)
    if(fullwindow) {
      mu <- sapply(lamlist, integral)
    } else {
      mu <- sapply(lamlist, integral, domain=w.cond)
    }
    p <- exp(n.cond * log(mu/n.cond) + n.cond - mu)
    phistory <- c(phistory, p)
    mhistory <- c(mhistory, mu)
    ## accept/reject
    accept <- (runif(length(p)) < p)
    if(any(accept)) {
      jaccept <- which(accept)
      if(length(jaccept) > nremaining)
        jaccept <- jaccept[seq_len(nremaining)]
      naccepted <- length(jaccept)
      if(verbose)
        splat("Accepted the",
              commasep(ordinal(ntried + jaccept)),
              ngettext(naccepted, "proposal", "proposals"))
      nremaining <- nremaining - naccepted
      for(j in jaccept) {
        lamj <- lamlist[[j]]
        if(min(lamj) < 0)
          lamj <- eval.im(pmax(lamj, 0))
        if(fullwindow) {
          Y <- rpoint(n.cond, lamj, win=w.sim, forcewin=TRUE)
        } else {
          lamj.cond <- lamj[w.cond, drop=FALSE, tight=TRUE]
          lamj.free <- lamj[w.free, drop=FALSE, tight=TRUE]
          Ycond <- rpoint(n.cond, lamj.cond, win=w.cond)
          Yfree <- rpoispp(lamj.free)
          Y <- superimpose(Ycond, Yfree, W=w.sim)
        }
        results <- append(results, list(Y))
      }
    }
    ntried <- ntried + nchunk
    if(ntried >= giveup && nremaining > 0) {
      message(paste("Gave up after", ntried,
                    "proposals with", nsim - nremaining, "accepted"))
      message(paste("Mean acceptance probability =",
                    signif(mean(phistory), 3)))
      break
    }
  }
  nresults <- length(results)
  results <- simulationresult(results, nresults, drop)
  attr(results, "history") <- data.frame(mu=mhistory, p=phistory)
  if(verbose && nresults == nsim)
    splat("Mean acceptance probability", signif(mean(phistory), 3))
  return(results)
}
