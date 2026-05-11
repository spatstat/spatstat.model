#'   makepspace.R
#'
#'   Create 'pspace' argument for kppm
#'
#'   Including default penalty for cluster scale
#'
#'   $Revision: 1.13 $ $Date: 2026/05/11 08:04:39 $
#' 
#'   Copyright (c) Tilman Davies, Martin Hazelton and Adrian Baddeley 2022-2026
#'  GNU Public Licence >= 2.0


make.pspace <- function(...,
                        fixedpar=NULL,
                        canonical=FALSE,
                        adjusted=FALSE,
                        trace=FALSE,
                        save=trajectory,
                        trajectory=FALSE,
                        nhalfgrid=NULL,
                        strict=TRUE,
                        penalised=NULL,
                        penalty=NULL,
                        penal.args=NULL,
                        tau=NULL,
                        clusters="Thomas",
                        fitmethod=c("mincon", "clik2", "waag", "palm"),
                        flatness=2,
                        C0factor=0.05,
                        xval=FALSE,
                        xval.args=list(),
                        debug=FALSE,
                        transfo=NULL) {
  ## validate cluster model
  info <- spatstatClusterModelInfo(clusters)
  ## validate fixed parameters
  if(length(fixedpar)) {
    nam <- names(fixedpar)
    if(!all(nzchar(nam))) stop("fixedpar must be a named list", call.=FALSE)
    pn.native    <- info$parnames
    pn.generic   <- info$pn.generic
    pn.canonical <- info$pn.canonical
    unknown <- setdiff(nam, c(pn.native, pn.generic, pn.canonical))
    if(length(unknown)) {
      if(length(unknown) == 1) {
        stop(paste("Name of fixed parameter", sQuote(unknown), "is not recognised"), call.=FALSE)
      } else {
        stop(paste("Names of fixed parameters", commasep(sQuote(unknown)),
                   "are not recognised"), call.=FALSE)
      }
    }
    can.be.native <- all(nam %in% pn.native)
    can.be.generic <- all(nam %in% info$pn.generic)
    can.be.canon <- all(nam %in% info$pn.canonical)
    if(!(can.be.native || can.be.generic || can.be.canon))
      stop(paste("The combination of parameters", sQuote(nam),
                 "does not match any of the standard parametrisations"),
           call.=FALSE)
    must.be.canon <- can.be.canon && !can.be.native && !can.be.generic
    if(missing(canonical)) {
      ## use canonical parametrisation only when it is implied 
      canonical <- must.be.canon
    } else if(isFALSE(canonical) && must.be.canon) {
      ## given canonical=FALSE but fixedpar implies canonical parameters
      warning("Re-setting canonical=TRUE (implied by fixedpar)",
              call.=FALSE)
      canonical <- TRUE
    } else if(isTRUE(canonical) && !can.be.canon) {
      ## given canonical=TRUE but fixedpar does not match canonical parameters
      warning("Re-setting canonical=FALSE (implied by fixedpar)",
              call.=FALSE)
      canonical <- FALSE
    }
  }
  ## assemble all recognised arguments
  p <- list(fixedpar  = fixedpar,
            canonical = isTRUE(canonical),
            adjusted  = isTRUE(adjusted),
            trace     = isTRUE(trace),
            save      = isTRUE(save),
            nhalfgrid = nhalfgrid,
            strict    = !isFALSE(strict),
            xval      = isTRUE(xval),
            xval.args = as.list(xval.args),
            debug     = debug,
            transfo   = transfo)
    
  ## penalise cluster scale?
  penalised <- isTRUE(penalised)
  if(is.function(penalty)) {
    ## user-specified penalty
    penalised <- TRUE
  } else if(penalised && is.null(penalty)) {
    ## default penalty function
    if(flatness <= 0 || flatness %% 2 != 0)
      stop("'flatness' of penalty must be even and positive", call.=FALSE)
    ## penalty is applied to generic 'scale' parameter
    native2generic <- info[["native2generic"]]
    if(!is.function(native2generic))
      stop(paste("Unable to determine generic scale parameter, for clusters=", sQuote(clusters)),
           call.=FALSE)
    HazeltonPenalty <- function(par, A, B, flatness) {
      s <- native2generic(par)[["scale"]] 
      u <- sqrt(s/A) - sqrt(B/s)
      v <- 1 - sqrt(B/A)
      (u/v)^flatness
    }
    penalty <- HazeltonPenalty
  }
  if(penalised) {
    ## compute arguments of penalty
    if(is.null(penal.args)) {
      penal.args <- function(X, rho=flatness) {
        nnd <- nndist(X)
        p <- list(A = median(nnd),
                  B = diameter(Window(X))/2,
                  flatness=rho)
        return(p)
      }
    }
    if(is.null(tau)) {
      fitmethod <- match.arg(fitmethod)
      tau <- switch(fitmethod,
                    mincon = function(X, poisval, f=C0factor) { f * poisval },
                    palm = 1,
                    waag = 1,
                    clik2 = 1)
    }
    ## add arguments of penalty to pspace
    p <- append(p,
                list(penalty    = penalty,
                     penal.args = penal.args,
                     tau        = tau))
  }
  return(p)
}


