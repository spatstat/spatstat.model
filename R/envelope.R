#
#   envelope.R
#
#   computes simulation envelopes 
#
#   $Revision: 2.134 $  $Date: 2026/02/28 12:22:18 $
#

## Code for envelope() and envelope.ppp() is moved to spatstat.explore


## //////////////// METHODS FOR "ppm" //////////////////////////

envelope.ppm <- 
  function(Y, fun=Kest, nsim=99, nrank=1, ..., 
           funargs=list(), funYargs=funargs,
           simulate=NULL, fix.n=FALSE, fix.marks=FALSE,
           verbose=TRUE, clipdata=TRUE, 
           start=NULL,
           control=update(default.rmhcontrol(Y), nrep=nrep), nrep=1e5, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL, 
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE, 
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL,
           maxnerr=nsim, rejectNA=FALSE, silent=FALSE,
           do.pwrong=FALSE,
           envir.simul=NULL)
{
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kest
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())

  # Extract data pattern X from fitted model Y
  X <- data.ppm(Y)
  
  if(is.null(simulate)) {
    #' Simulated realisations of the fitted model Y will be generated
    simrecipe <- make.simulrecipe(Y, envir.here,
                                  fix.n=fix.n,
                                  fix.marks=fix.marks,
                                  start=start,
                                  control=control,
                                  nrep=nrep)
  } else {
    #' Simulations are determined by 'simulate' argument
    #' Processing is deferred to envelopeEngine
    simrecipe <- simulate
  }
  
  envelopeEngine(X=X, fun=fun, simul=simrecipe, 
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=clipdata,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp, 
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname,
                 maxnerr=maxnerr, rejectNA=rejectNA, silent=silent,
                 cl=cl, envir.user=envir.user, do.pwrong=do.pwrong)
}

make.simulrecipe.ppm <- function(object, envir, ...,
                                 fix.n=FALSE, fix.marks=FALSE,
                                 start=NULL,
                                 control=update(default.rmhcontrol(object),
                                                nrep=nrep),
                                 nrep=1e5) {
  #' Simulated realisations of the fitted model will be generated
  model <- object
  pois <- is.poisson(model)
  csr <- is.stationary(model) && pois
  type <- if(csr) "csr" else "rmh"
  X <- response(model)
  #' Set up parameters for rmh
  rmodel <- rmhmodel(model, verbose=FALSE)
  if(is.null(start))
    start <- list(n.start=npoints(X))
  rstart <- rmhstart(start)
  rcontr <- rmhcontrol(control)
  if(fix.marks) {
    rcontr <- update(rcontr, fixall=TRUE, p=1, expand=1)
    nst <- if(is.multitype(X)) table(marks(X)) else npoints(X)
    rstart <- update(rstart, n.start=nst)
    constraints <- "with fixed number of points of each type"
  } else if(fix.n) {
    rcontr <- update(rcontr, p=1, expand=1)
    rstart <- update(rstart, n.start=X$n)
    constraints <- "with fixed number of points"
  } else constraints <- ""
  #' pre-digest arguments
  rmhinfolist <- rmh(rmodel, rstart, rcontr, preponly=TRUE, verbose=FALSE)
  assign("rmhinfolist", rmhinfolist, envir=envir)
  #' expression that will be evaluated
  simexpr <- expression(rmhEngine(rmhinfolist, verbose=FALSE))
  rlz <- paste("simulated realisations of",
               if(csr) "CSR" else paste("fitted",
                                        if(pois) "Poisson" else "Gibbs",
                                        "process"))
  simrecipe <- simulrecipe(type  = type,
                           expr  = simexpr,
                           envir = envir,
                           csr   = csr,
                           pois  = pois,
                           constraints = constraints,
                           value="result of rmh()",
                           realisations=rlz)
  return(simrecipe)
}

## //////////////// METHODS FOR "kppm" //////////////////////////

envelope.kppm <-
  function(Y, fun=Kest, nsim=99, nrank=1, ..., 
           funargs=list(), funYargs=funargs,
           simulate=NULL, verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2, Yname=NULL, 
           maxnerr=nsim, rejectNA=FALSE, silent=FALSE,
           do.pwrong=FALSE, envir.simul=NULL)
{
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kest
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())
  
  # Extract data pattern X from fitted model Y
  X <- Y$X
  
  if(is.null(simulate)) {
    #' Simulated realisations of the fitted model Y will be generated
    simrecipe <- make.simulrecipe(Y, envir.here)
  } else {
    #' Simulations are determined by 'simulate' argument
    #' Processing is deferred to envelopeEngine
    simrecipe <- simulate
  }
  envelopeEngine(X=X, fun=fun, simul=simrecipe, 
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=clipdata,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp,
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname, 
                 maxnerr=maxnerr, rejectNA=rejectNA, silent=silent,
                 cl=cl, envir.user=envir.user, do.pwrong=do.pwrong)

}

make.simulrecipe.kppm <- function(object, envir, ...) {
  #' Simulated realisations of the fitted model Y
  #' will be generated using simulate.kppm
  kmodel <- object
  assign("kmodel", kmodel, envir=envir)
  #' expression that will be evaluated
  simexpr <- expression(simulate(kmodel, drop=TRUE))
  rlz <- paste("simulated realisations of fitted",
               if(is.poissonclusterprocess(kmodel)) "cluster" else "Cox",
               "process")
  simrecipe <- simulrecipe(type = "kppm",
                           expr = simexpr,
                           envir = envir,
                           csr   = FALSE,
                           pois  = FALSE,
                           value="result of simulate.kppm()",
                           realisations=rlz)
  return(simrecipe)
}

## //////////////// METHODS FOR "dppm" //////////////////////////

envelope.dppm <-
  function(Y, fun=Kest, nsim=99, nrank=1, ..., 
           funargs=list(), funYargs=funargs,
           simulate=NULL, verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2, Yname=NULL, 
           maxnerr=nsim, rejectNA=FALSE, silent=FALSE,
           do.pwrong=FALSE, envir.simul=NULL)
{
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kest
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())
  
  # Extract data pattern X from fitted model Y
  X <- Y$X
  
  if(is.null(simulate)) {
    #' Simulated realisations of the fitted model Y will be generated
    simrecipe <- make.simulrecipe(Y, envir.here)
  } else {
    #' Simulations are determined by 'simulate' argument
    #' Processing is deferred to envelopeEngine
    simrecipe <- simulate
  }
  envelopeEngine(X=X, fun=fun, simul=simrecipe, 
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=clipdata,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp,
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname, 
                 maxnerr=maxnerr, rejectNA=rejectNA, silent=silent,
                 cl=cl, envir.user=envir.user, do.pwrong=do.pwrong)

}

make.simulrecipe.dppm <- function(object, envir, ...) {
  #' Simulated realisations of the fitted model Y
  #' will be generated using simulate.dppm
  dmodel <- object
  assign("dmodel", dmodel, envir=envir)
  #' expression that will be evaluated
  simexpr <- expression(simulate(dmodel, drop=TRUE))
  rlz <- "simulated realisations of fitted determinantal point process model"
  simrecipe <- simulrecipe(type = "dppm",
                           expr = simexpr,
                           envir = envir,
                           csr   = FALSE,
                           pois  = FALSE,
                           value="result of simulate.dppm()",
                           realisations=rlz)
  return(simrecipe)
}

make.simulrecipe.detpointprocfamily <- function(object, envir, ..., W) {
  #' Simulated realisations of the fitted model Y
  #' will be generated using simulate.detpointprocfamily
  assign("simprocess", object,     envir=envir)
  assign("simdomain",       W, envir=envir)
  #' expression that will be evaluated
  simexpr <- expression(simulate(simprocess, W=simdomain, drop=TRUE))
  rlz <- "simulated realisations of determinantal point process family"  
  simrecipe <- simulrecipe(type = "process",
                           expr = simexpr,
                           envir = envir,
                           csr   = FALSE,
                           pois  = FALSE,
                           value="result of simulate.detpointprocfamily()",
                           realisations=rlz)
  return(simrecipe)
}

## //////////////// METHODS FOR "slrm" //////////////////////////

envelope.slrm <-
  function(Y, fun=Kest, nsim=99, nrank=1, ..., 
           funargs=list(), funYargs=funargs,
           simulate=NULL, verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2, Yname=NULL, 
           maxnerr=nsim, rejectNA=FALSE, silent=FALSE,
           do.pwrong=FALSE, envir.simul=NULL)
{
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kest
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())
  
  # Extract data pattern X from fitted model Y
  X <- response(Y)
  
  if(is.null(simulate)) {
    #' Simulated realisations of the fitted model Y
    #' will be generated using simulate.slrm
    simrecipe <- make.simulrecipe(Y, envir.here)
  } else {
    # ...................................................
    # Simulations are determined by 'simulate' argument
    # Processing is deferred to envelopeEngine
    simrecipe <- simulate
  }
  envelopeEngine(X=X, fun=fun, simul=simrecipe, 
                 nsim=nsim, nrank=nrank, ...,
                 funargs=funargs, funYargs=funYargs,
                 verbose=verbose, clipdata=clipdata,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp,
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname, 
                 maxnerr=maxnerr, rejectNA=rejectNA, silent=silent,
                 cl=cl, envir.user=envir.user, do.pwrong=do.pwrong)

}

make.simulrecipe.slrm <- function(object, envir, ...) {
  #' Simulated realisations of the fitted model Y
  #' will be generated using simulate.slrm
  model <- object
  assign("model", model, envir=envir)
  #' expression that will be evaluated
  simexpr <- expression(simulate(model, drop=TRUE))
  rlz <- "simulated realisations of spatial logistic regression model"
  simrecipe <- simulrecipe(type = "slrm",
                           expr = simexpr,
                           envir = envir,
                           csr   = FALSE,
                           pois  = FALSE,
                           value = "result of simulate.slrm()",
                           realisations = rlz)

  return(simrecipe)
}



