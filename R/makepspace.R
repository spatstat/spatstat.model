#'   makepspace.R
#'
#'   Create 'pspace' argument for kppm
#'
#'   Including default penalty for cluster scale
#'
#'   $Revision: 1.11 $ $Date: 2026/05/04 06:49:57 $
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
    can.be.native <- all(nam %in% info$parnames)
    can.be.canon <- all(nam %in% info$pn.canonical)
    must.be.native <- can.be.native && !can.be.canon
    must.be.canon <- can.be.canon && !can.be.native
    if(missing(canonical)) {
      canonical <- must.be.canon
    } else if(isFALSE(canonical) && must.be.canon) {
      ## given canonical=FALSE but fixedpar implies canonical parameters
      warning("Re-setting canonical=TRUE (implied by fixedpar)",
              call.=FALSE)
      canonical <- TRUE
    } else if(isTRUE(canonical) && must.be.native) {
      ## given canonical=TRUE but fixedpar implies native parameters
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


