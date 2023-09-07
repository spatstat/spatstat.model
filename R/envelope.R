#
#   envelope.R
#
#   computes simulation envelopes 
#
#   $Revision: 2.125 $  $Date: 2023/08/15 08:07:52 $
#

## Code for envelope() and envelope.ppp() is moved to spatstat.explore


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
           envir.simul=NULL) {
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kest
  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())

  # Extract data pattern X from fitted model Y
  X <- data.ppm(Y)
  
  if(is.null(simulate)) {
    # ...................................................
    # Simulated realisations of the fitted model Y
    # will be generated
    pois <- is.poisson(Y)
    csr <- is.stationary(Y) && pois
    type <- if(csr) "csr" else "rmh"
    # Set up parameters for rmh
    rmodel <- rmhmodel(Y, verbose=FALSE)
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
    # pre-digest arguments
    rmhinfolist <- rmh(rmodel, rstart, rcontr, preponly=TRUE, verbose=FALSE)
    # expression that will be evaluated
    simexpr <- expression(rmhEngine(rmhinfolist, verbose=FALSE))
    dont.complain.about(rmhinfolist)
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type  = type,
                             expr  = simexpr,
                             envir = envir.here,
                             csr   = csr,
                             pois  = pois,
                             constraints = constraints)
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
    # Simulated realisations of the fitted model Y
    # will be generated using simulate.kppm
    kmodel <- Y
    # expression that will be evaluated
    simexpr <- expression(simulate(kmodel)[[1L]])
    dont.complain.about(kmodel)
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "kppm",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = FALSE,
                             pois  = FALSE)
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
    # Simulated realisations of the fitted model Y
    # will be generated using simulate.slrm
    smodel <- Y
    # expression that will be evaluated
    simexpr <- expression(simulate(smodel)[[1L]])
    dont.complain.about(smodel)
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "slrm",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = FALSE,
                             pois  = FALSE)
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


