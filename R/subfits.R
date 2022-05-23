#
#
#  $Revision: 1.55 $   $Date: 2021/12/29 07:50:25 $
#
#

subfits.new <- local({

  subfits.new <- function(object, what="models", verbose=FALSE) {
    stopifnot(inherits(object, "mppm"))

    what <- match.arg(what, c("models", "interactions", "basicmodels"))

    if(what == "interactions")
      return(subfits.old(object, what, verbose))
  
    ## extract stuff
    announce <- if(verbose) Announce else Ignore
    
    announce("Extracting stuff...")
    fitter   <- object$Fit$fitter
    FIT      <- object$Fit$FIT
    trend    <- object$trend
    random   <- object$random
    info     <- object$Info
    npat     <- object$npat
    Inter    <- object$Inter
    interaction <- Inter$interaction
    itags    <- Inter$itags
    Vnamelist <- object$Fit$Vnamelist
    Isoffsetlist <- object$Fit$Isoffsetlist
    has.design <- info$has.design
#    has.random <- info$has.random
    announce("done.\n")

    ## fitted parameters
    coefs.full <- coef(object)
    if(is.null(dim(coefs.full))) {
      ## fixed effects model: replicate vector to matrix
      coefs.names <- names(coefs.full)
      coefs.full <- matrix(coefs.full, byrow=TRUE,
                           nrow=npat, ncol=length(coefs.full),
                           dimnames=list(NULL, coefs.names))
    } else {
      ## random/mixed effects model: coerce to matrix
      coefs.names <- colnames(coefs.full)
      coefs.full <- as.matrix(coefs.full)
    }
    
    ## determine which interaction(s) are active on each row
    announce("Determining active interactions...")
    active <- active.interactions(object)
    announce("done.\n")

    ## exceptions
    if(any(rowSums(active) > 1))
      stop(paste("subfits() is not implemented for models",
                 "in which several interpoint interactions",
                 "are active on the same point pattern"))
    if(!is.null(random) && any(variablesinformula(random) %in% itags))
      stop(paste("subfits() is not yet implemented for models",
                 "with random effects that involve",
                 "the interpoint interactions"))
  
    ## implied coefficients for each active interaction
    announce("Computing implied coefficients...")
    implcoef <- list()
    for(tag in itags) {
      announce(tag)
      implcoef[[tag]] <- impliedcoefficients(object, tag)
      announce(", ")
    }
    announce("done.\n")

    ## Fisher information and vcov
    fisher <- varcov <- NULL
    if(what == "models") {
      announce("Fisher information...")
      fisher   <- vcov(object, what="fisher", err="null")
      varcov   <- try(solve(fisher), silent=TRUE)
      if(inherits(varcov, "try-error"))
        varcov <- NULL
      announce("done.\n")
    } 
  
    ## Extract data frame 
    announce("Extracting data...")
    datadf   <- object$datadf
    rownames <- object$Info$rownames
    announce("done.\n")

    ## set up lists for results 
    models <- rep(list(NULL), npat)
    interactions <- rep(list(NULL), npat)
    
    ## interactions
    announce("Determining interactions...")
    pstate <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      ## Find relevant interaction
      acti <- active[i,]
      nactive <- sum(acti)
      interi <- if(nactive == 0) Poisson() else interaction[i, acti, drop=TRUE]
      tagi <- names(interaction)[acti]
      ## Find relevant coefficients
      coefs.avail  <- coefs.full[i,]
      names(coefs.avail) <- coefs.names
      if(nactive == 1) {
        ic <- implcoef[[tagi]]
        coefs.implied <- ic[i, ,drop=TRUE]
        names(coefs.implied) <- colnames(ic)
        ## overwrite any existing values of coefficients; add new ones.
        coefs.avail[names(coefs.implied)] <- coefs.implied
      }
      ## create fitted interaction with these coefficients
      vni <- if(nactive > 0) Vnamelist[[tagi]] else character(0)
      iso <- if(nactive > 0) Isoffsetlist[[tagi]] else logical(0)
      interactions[[i]] <- fii(interi, coefs.avail, vni, iso)
    }
    announce("Done!\n")
    names(interactions) <- rownames

    ##
    if(what=="interactions") 
      return(interactions)
  
    ## Extract data required to reconstruct complete model fits
    announce("Extracting more data...")
    data  <- object$data
    Y     <- object$Y
    Yname <- info$Yname
    moadf <- object$Fit$moadf
    fmla  <- object$Fit$fmla
    ## deal with older formats of mppm
    if(is.null(Yname)) Yname <- info$Xname
    if(is.null(Y)) Y <- data[ , Yname, drop=TRUE]
    ## 
    used.cov.names <- info$used.cov.names
    has.covar <- info$has.covar
    if(has.covar) {
      covariates.hf <- data[, used.cov.names, drop=FALSE]
      dfvar <- used.cov.names %in% names(datadf)
    }
    announce("done.\n")

    ## Construct template for fake ppm object
    spv <- package_version(versionstring.spatstat())
    fake.version <- list(major=spv$major,
                         minor=spv$minor,
                         release=spv$patchlevel,
                         date="$Date: 2021/12/29 07:50:25 $")
    fake.call <- call("cannot.update", Q=NULL, trend=trend,
                      interaction=NULL, covariates=NULL,
                      correction=object$Info$correction,
                      rbord     = object$Info$rbord)
    fakemodel <- list(
                    method       = "mpl",
                    fitter       = fitter,
                    coef         = coef(object),
                    trend        = object$trend,
                    interaction  = NULL,
                    fitin        = NULL,
                    Q            = NULL,
                    maxlogpl     = NA,
                    internal     = list(glmfit = FIT,
                                        glmdata  = NULL,
                                        Vnames   = NULL,
                                        IsOffset  = NULL,
                                        fmla     = fmla,
                                        computed = list()),
                    covariates   = NULL,
                    correction   = object$Info$correction,
                    rbord        = object$Info$rbord,
                    version      = fake.version,
                    problems     = list(),
                    fisher       = fisher,
                    varcov       = varcov,
                    call         = fake.call,
                    callstring   = "cannot.update()",
                    fake         = TRUE)
    class(fakemodel) <- "ppm"

    ## Loop through point patterns
    announce("Generating models for each row...")
    pstate <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      Yi <- Y[[i]]
      Wi <- if(is.ppp(Yi)) Yi$window else Yi$data$window
      ## assemble relevant covariate images
      covariates <-
        if(has.covar) covariates.hf[i, , drop=TRUE, strip=FALSE] else NULL
      if(has.covar && has.design) 
        ## Convert each data frame covariate value to an image
        covariates[dfvar] <- lapply(covariates[dfvar], as.im, W=Wi)

      ## Extract relevant interaction
      finte <- interactions[[i]]
      inte  <- finte$interaction
      if(is.poisson.interact(inte)) inte <- NULL
      Vnames <- finte$Vnames
      if(length(Vnames) == 0) Vnames <- NULL
      IsOffset <- finte$IsOffset
      if(length(IsOffset) == 0) IsOffset <- NULL
      
      ## Construct fake ppm object
      fakemodel$interaction <- inte
      fakemodel$fitin       <- finte
      fakemodel$Q           <- Yi
      fakemodel$covariates  <- covariates
      fakemodel$internal$glmdata <- moadf[moadf$id == i, ]
      fakemodel$internal$Vnames  <- Vnames
      fakemodel$internal$IsOffset <- IsOffset

      fake.call$Q <- Yi
      fake.call$covariates <- covariates
      fakemodel$call <- fake.call
      fakemodel$callstring <- short.deparse(fake.call)
      
      ## store in list
      models[[i]] <- fakemodel
    }
    announce("done.\n")
    names(models) <- rownames
    models <- as.anylist(models)
    return(models)
  }

  Announce <- function(...) cat(...)

  Ignore <- function(...) { NULL }

  subfits.new
})



## /////////////////////////////////////////////////////

subfits <-
subfits.old <- local({
    
  subfits.old <- function(object, what="models", verbose=FALSE, new.coef=NULL) {
    stopifnot(inherits(object, "mppm"))
    what <- match.arg(what, c("models","interactions", "basicmodels"))
    ## extract stuff
    announce <- if(verbose) Announce else Ignore
    
    announce("Extracting stuff...")
    trend    <- object$trend
    random   <- object$random
    use.gam  <- object$Fit$use.gam
    info     <- object$Info
    npat     <- object$npat
    Inter    <- object$Inter
    interaction <- Inter$interaction
    itags    <- Inter$itags
    Vnamelist <- object$Fit$Vnamelist
    Isoffsetlist <- object$Fit$Isoffsetlist
    has.design <- info$has.design
    has.random <- info$has.random
    moadf    <- object$Fit$moadf
    announce("done.\n")

    ## levels of any factors
    levelslist <- lapply(as.list(moadf), levelsAsFactor)
    isfactor <- !sapply(levelslist, is.null)
    
    ## fitted parameters
    coefs.full <- new.coef %orifnull% coef(object)
    if(is.null(dim(coefs.full))) {
      ## fixed effects model: replicate vector to matrix
      coefs.names <- names(coefs.full)
      coefs.full <- matrix(coefs.full, byrow=TRUE,
                           nrow=npat, ncol=length(coefs.full),
                           dimnames=list(NULL, coefs.names))
    } else {
      ## random/mixed effects model: coerce to matrix
      coefs.names <- colnames(coefs.full)
      coefs.full <- as.matrix(coefs.full)
    }
  
    ## determine which interaction(s) are active on each row
    announce("Determining active interactions...")
    active <- active.interactions(object)
    announce("done.\n")

    ## exceptions
    if(any(rowSums(active) > 1))
      stop(paste("subfits() is not implemented for models",
                 "in which several interpoint interactions",
                 "are active on the same point pattern"))
    if(!is.null(random) && any(variablesinformula(random) %in% itags))
      stop(paste("subfits() is not yet implemented for models",
                 "with random effects that involve",
                 "the interpoint interactions"))
  
    ## implied coefficients for each active interaction
    announce("Computing implied coefficients...")
    implcoef <- list()
    for(tag in itags) {
      announce(tag)
      implcoef[[tag]] <- impliedcoefficients(object, tag, new.coef=new.coef)
      announce(", ")
    }
    announce("done.\n")

    ## This code is currently not usable because the mapping is wrong
    reconcile <- FALSE
    if(reconcile) {
      ## determine which coefficients of main model are interaction terms
      announce("Identifying interaction coefficients...")
      md <- model.depends(object$Fit$FIT)
      usetags <- unlist(lapply(implcoef, colnames))
      isVname <- apply(md[, usetags, drop=FALSE], 1, any)
      mainVnames <- row.names(md)[isVname]
      announce("done.\n")
    }

    ## Fisher information and vcov
    fisher <- varcov <- NULL
    if(what == "models") {
      announce("Fisher information...")
      fisher   <- vcov(object, what="fisher", err="null", new.coef=new.coef)
      ## note vcov.mppm calls subfits(what="basicmodels") to avoid infinite loop
      varcov   <- try(solve(fisher), silent=TRUE)
      if(inherits(varcov, "try-error"))
        varcov <- NULL
      announce("done.\n")
    }
  
    ## Extract data frame 
    announce("Extracting data...")
    datadf   <- object$datadf
    rownames <- object$Info$rownames
    announce("done.\n")

    ## set up list for results 
    results <- rep(list(NULL), npat)
  
    if(what == "interactions") {
      announce("Determining interactions...")
      pstate <- list()
      for(i in 1:npat) {
        if(verbose) pstate <- progressreport(i, npat, state=pstate)
        ## Find relevant interaction
        acti <- active[i,]
        nactive <- sum(acti)
        interi <- if(nactive == 0) Poisson() else
                  interaction[i, acti, drop=TRUE]
        tagi <- names(interaction)[acti]
        ## Find relevant coefficients
        coefs.avail  <- coefs.full[i,]
        names(coefs.avail) <- coefs.names
        if(nactive == 1) {
          ic <- implcoef[[tagi]]
          coefs.implied <- ic[i, ,drop=TRUE]
          names(coefs.implied) <- colnames(ic)
          ## overwrite any existing values of coefficients; add new ones.
          coefs.avail[names(coefs.implied)] <- coefs.implied
        }
        ## create fitted interaction with these coefficients
        vni <- if(nactive > 0) Vnamelist[[tagi]] else character(0)
        iso <- if(nactive > 0) Isoffsetlist[[tagi]] else logical(0)
        results[[i]] <- fii(interi, coefs.avail, vni, iso)
      }
      announce("Done!\n")
      names(results) <- rownames
      return(results)
    }
  
    ## Extract data required to reconstruct complete model fits
    announce("Extracting more data...")
    data  <- object$data
    Y     <- object$Y
    Yname <- info$Yname
    ## deal with older formats of mppm
    if(is.null(Yname)) Yname <- info$Xname
    if(is.null(Y)) Y <- data[ , Yname, drop=TRUE]
    ##
    used.cov.names <- info$used.cov.names
    has.covar <- info$has.covar
    if(has.covar) {
      covariates.hf <- data[, used.cov.names, drop=FALSE]
      dfvar <- used.cov.names %in% names(datadf)
    }
    announce("done.\n")
  
    ## Loop through point patterns
    announce("Looping through rows...")
    pstate <- list()
    for(i in 1:npat) {
      if(verbose) pstate <- progressreport(i, npat, state=pstate)
      Yi <- Y[[i]]
      Wi <- if(is.ppp(Yi)) Yi$window else Yi$data$window
      ## assemble relevant covariate images
      scrambled <- FALSE
      if(!has.covar) { 
        covariates <- NULL
      } else {
        covariates <- covariates.hf[i, , drop=TRUE, strip=FALSE] 
        if(has.design) {
          ## Convert each data frame covariate value to an image
          imrowi <- lapply(covariates[dfvar], as.im, W=Wi)
          ## Problem: constant covariate leads to singular fit
          ## --------------- Hack: ---------------------------
          scrambled <- TRUE
          ##  Construct fake data by resampling from possible values
          covar.vals <- lapply(as.list(covariates[dfvar, drop=FALSE]), possible)
          fake.imrowi <- lapply(covar.vals, scramble, W=Wi, Y=Yi$data)
          ## insert fake data into covariates 
          covariates[dfvar] <- fake.imrowi
          ## ------------------ end hack ----------------------------
        }
        ## identify factor-valued spatial covariates
        spatialfactors <- !dfvar & isfactor[names(covariates)]
        if(any(spatialfactors)) {
          ## problem: factor levels may be dropped
          ## more fakery...
          scrambled <- TRUE
          spfnames <- names(spatialfactors)[spatialfactors]
          covariates[spatialfactors] <-
            lapply(levelslist[spfnames],
                   scramble, W=Wi, Y=Yi$data)
        }
      }
      ## Fit ppm to data for case i only
      ## using relevant interaction
      acti <- active[i,]
      nactive <- sum(acti)
      if(nactive == 1){
        interi <- interaction[i, acti, drop=TRUE] 
        tagi <- names(interaction)[acti]
        fiti <- PiPiM(Yi, trend, interi, covariates=covariates,
                      allcovar=has.random,
                      use.gam=use.gam,
                      vnamebase=tagi, vnameprefix=tagi)
      } else {
        fiti <- PiPiM(Yi, trend, Poisson(), covariates=covariates,
                      allcovar=has.random,
                      use.gam=use.gam)
      }
      fiti$scrambled <- scrambled
      ## fiti determines which coefficients are required 
      coefi.fitted <- fiti$coef
      coefnames.wanted <- coefnames.fitted <- names(coefi.fitted)
      ## reconcile interaction coefficient names
      if(reconcile) {
        coefnames.translated <- coefnames.wanted
        ma <- match(coefnames.fitted, fiti$internal$Vnames)
        hit <- !is.na(ma)
        if(any(hit)) 
          coefnames.translated[hit] <- mainVnames[ ma[hit] ]
      }
      ## take the required coefficients from the full mppm fit
      coefs.avail  <- coefs.full[i,]
      names(coefs.avail) <- coefs.names
      if(nactive == 1) {
        ic <- implcoef[[tagi]]
        coefs.implied <- ic[i, ,drop=TRUE]
        names(coefs.implied) <- colnames(ic)
        ## overwrite any existing values of coefficients; add new ones.
        coefs.avail[names(coefs.implied)] <- coefs.implied
      }
      ## check
      if(!all(coefnames.wanted %in% names(coefs.avail))) 
        stop("Internal error: some fitted coefficients not accessible")

      ## hack entries in ppm object 
      fiti$method <- "mppm"
      ## reset coefficients
      coefi.new <- coefs.avail[coefnames.wanted]
      fiti$coef.orig <- coefi.fitted ## (detected by summary.ppm, predict.ppm)
      fiti$theta <- fiti$coef <- coefi.new
      ## reset interaction coefficients in 'fii' object
      coef(fiti$fitin)[] <- coefi.new[names(coef(fiti$fitin))]
      ## ... and replace fake data by true data
      if(has.design) {
        fiti$internal$glmdata.scrambled <- gd <- fiti$internal$glmdata
        fixnames <- intersect(names(imrowi), colnames(gd))
        for(nam in fixnames) {
          fiti$covariates[[nam]] <- imrowi[[nam]]
          fiti$internal$glmdata[[nam]] <- data[i, nam, drop=TRUE]
        }
      }
      ## Adjust rank of glm fit object
#      fiti$internal$glmfit$rank <- FIT$rank 
      fiti$internal$glmfit$rank <- sum(is.finite(fiti$coef))
      ## Fisher information and variance-covariance if known
      ## Extract submatrices for relevant parameters
      if(reconcile) {
        #' currently disabled because mapping is wrong
        if(!is.null(fisher)) {
          if(!reconcile) {
            fiti$fisher <-
              fisher[coefnames.wanted, coefnames.wanted, drop=FALSE]
          } else {
            fush <-
              fisher[coefnames.translated, coefnames.translated, drop=FALSE]
            dimnames(fush) <- list(coefnames.wanted, coefnames.wanted)
            fiti$fisher <- fush
          }
        }
        if(!is.null(varcov)) {
          if(!reconcile) {
            fiti$varcov <-
              varcov[coefnames.wanted, coefnames.wanted, drop=FALSE]
          } else {
            vc <- varcov[coefnames.translated, coefnames.translated, drop=FALSE]
            dimnames(vc) <- list(coefnames.wanted, coefnames.wanted)
            fiti$varcov <- vc
          }
        }
      }
      ## store in list
      results[[i]] <- fiti
    }
    announce("done.\n")
    names(results) <- rownames
    results <- as.anylist(results)
    return(results)
  }

  PiPiM <- function(Y, trend, inter, covariates, ...,
                    allcovar=FALSE, use.gam=FALSE,
                    vnamebase=c("Interaction", "Interact."),
                    vnameprefix=NULL) {
    # This ensures that the model is fitted in a unique environment
    # so that it can be updated later.
    force(Y)
    force(trend)
    force(inter)
    force(covariates)
    force(allcovar)
    force(use.gam)
    force(vnamebase)
    force(vnameprefix)
    feet <- ppm(Y, trend, inter, covariates=covariates,
                allcovar=allcovar, use.gam=use.gam,
                forcefit=TRUE, vnamebase=vnamebase, vnameprefix=vnameprefix)
    return(feet)
  }
  
  possible <- function(z) {
    if(!is.factor(z)) unique(z) else factor(levels(z), levels=levels(z))
  }
  
  scramble <- function(vals, W, Y) {
    W <- as.mask(W)
    npixels <- prod(W$dim)
    nvalues <- length(vals)
    npts <- npoints(Y)
    ## sample the possible values randomly at the non-data pixels
    sampled <- sample(vals, npixels, replace=TRUE)
    Z <- im(sampled, xcol=W$xcol, yrow=W$yrow)
    ## repeat the possible values cyclically at the data points
    if(npts >= 1)
      Z[Y] <- vals[1 + ((1:npts) %% nvalues)]
    return(Z)
  }

  Announce <- function(...) cat(...)

  Ignore <- function(...) { NULL }

  subfits.old
})

cannot.update <- function(...) {
  stop("This model cannot be updated")
}

mapInterVars <- function(object, subs=subfits(object), mom=model.matrix(object)) {
  #' Map the canonical variables of each 'subs[[i]]' to those of 'object'
  #' This is only needed for interaction variables.
  #'
  #' (1) Information about the full model
  #'     Names of interaction variables
  Vnamelist    <- object$Fit$Vnamelist
  Isoffsetlist <- object$Fit$Isoffsetlist
  #'     Dependence map of canonical variables on the original variables/interactions
  md <- model.depends(object$Fit$FIT)
  cnames <- rownames(md)
  #' (2) Information about the fit on each row
  #'     Identify the (unique) active interaction in each row
  activeinter <- active.interactions(object)
  #'     Determine which canonical variables of full model are active in each row
  mats <- split.data.frame(mom, object$Fit$moadf$id)
  activevars <- sapply(mats, function(df) { apply(df != 0, 2, any) })
  activevars <- if(ncol(mom) > 1) t(activevars) else matrix(activevars, ncol=1)
  if(ncol(activevars) != ncol(mom)) warning("Internal error: activevars columns do not match canonical variables")
  if(nrow(activevars) != length(mats)) warning("Internal error: activevars rows do not match hyperframe rows")
  #' (3) Process each submodel
  n <- length(subs)
  result <- rep(list(list()), n)
  for(i in seq_len(n)) {
    #' the submodel in this row
    subi <- subs[[i]]
    if(!is.poisson(subi)) {
      cnames.i <- names(coef(subi))
      #' the (unique) tag name of the interaction in this model
      tagi <- colnames(activeinter)[activeinter[i,]]
      #' the corresponding variable name(s) in glmdata and coef(subi)
      vni <- Vnamelist[[tagi]]
      iso <- Isoffsetlist[[tagi]]
      #' ignore offset variables
      vni <- vni[!iso]
      if(length(vni)) {
        #' retain only interaction variables
        e <- cnames.i %in% vni
        cnames.ie <- cnames.i[e]
        #' which coefficients of the full model are active in this row
        acti <- activevars[i,]
        #' for each interaction variable name in the submodel,
        #' find the coefficient(s) in the main model to which it contributes
        nie <- length(cnames.ie)
        cmap <- vector(mode="list", length=nie)
        names(cmap) <- cnames.ie
        for(j in seq_len(nie)) {
          cj <- cnames.ie[j]
          cmap[[j]] <- cnames[ md[,cj] & acti ]
        }
        #' The result 'cmap' is a named list of character vectors
        #' where each name is an interaction variable in subs[[i]]
        #' and the corresponding value is the vector of names of
        #' corresponding canonical variables in the full model
        result[[i]] <- cmap
      }
    }
  }
  return(result)
}
