#  Variance-covariance matrix for mppm objects
#
# $Revision: 1.23 $ $Date: 2021/12/29 07:50:32 $
#
#

vcov.mppm <- local({

  errhandler <- function(whinge, err) {
    switch(err,
           fatal=stop(whinge),
           warn={
             warning(whinge)
             return(NA)
           },
           null= return(NULL),
           stop(paste("Unrecognised option: err=", dQuote(err))))
  }
    
  vcov.mppm <- function(object, ..., what="vcov", err="fatal") {

    what <- match.arg(what,
                      c("vcov", "corr", "fisher", "Fisher", "internals", "all"))
    if(what == "Fisher") what <- "fisher"
    
    if(is.poisson.mppm(object) && object$Fit$fitter == "glm")
      return(vcmPois(object, ..., what=what, err=err))

    return(vcmGibbs(object, ..., what=what, err=err))
  }

  vcmPois <- function(object, ..., what, err,
                      nacoef.action=c("warn", "fatal", "silent"),
                      new.coef=NULL
                      ) {
    #' legacy algorithm for Poisson case

    #' detect NA coefficients
    if(missing(nacoef.action) && !missing(err) && !is.null(err)) {
      nacoef.action <- err
    } else {
      nacoef.action <- match.arg(nacoef.action)
    }
    if(!all(is.finite(coef(object)))) {
      gripe <- "Cannot compute variance; some coefficients are NA, NaN or Inf"
      switch(nacoef.action,
             fatal = stop(gripe, call.=FALSE),
             warn = warning(gripe, call.=FALSE),
             silent = {})
      return(NULL)
    }

    #' get to work
    gf <- object$Fit$FIT
    gd <- object$Fit$moadf
    wt <- gd$.mpl.W

    fo <- object$trend
    if(is.null(fo)) fo <- (~1)

    mof <- model.frame(fo, gd)
    mom <- model.matrix(fo, mof)
    momnames <- dimnames(mom)[[2]]

    ## fitted intensity
    if(!is.null(new.coef) && inherits(gf, c("gam", "lme"))) {
      warn.once("vcovGAMnew", "'new.coef' is not supported by vcov.mppm for GAM or LME models; ignored")
      new.coef <- NULL
    }
    fi <- if(is.null(new.coef)) fitted(gf) else GLMpredict(gf, gd, new.coef, changecoef=TRUE, type="response")

    fisher <- sumouter(mom, fi * wt)
    dimnames(fisher) <- list(momnames, momnames)

    switch(what,
           fisher = { return(fisher) },
           vcov   = {
             vc <- try(solve(fisher), silent=(err == "null"))
             if(inherits(vc, "try-error"))
               return(errhandler("Fisher information is singular", err))
             else
               return(vc)
           },
           corr={
             co <- try(solve(fisher), silent=(err == "null"))
             if(inherits(co, "try-error"))
               return(errhandler("Fisher information is singular", err))
             sd <- sqrt(diag(co))
             return(co / outer(sd, sd, "*"))
           })
  }

  vcmGibbs <- function(object, ..., what, err,
                       matrix.action=c("warn", "fatal", "silent"),
                       gam.action=c("warn", "fatal", "silent"),
                       logi.action=c("warn", "fatal", "silent"),
                       nacoef.action=c("warn", "fatal", "silent"),
                       new.coef=NULL
                       ) {
    if(!missing(err)) {
      if(err == "null") err <- "silent" 
      matrix.action <-
        if(missing(matrix.action)) err else match.arg(matrix.action)
      gam.action <- if(missing(gam.action)) err else match.arg(gam.action)
      logi.action <- if(missing(logi.action)) err else match.arg(logi.action)
      nacoef.action <- if(missing(nacoef.action)) err else match.arg(nacoef.action)
    } else {
      matrix.action <- match.arg(matrix.action)
      gam.action <- match.arg(gam.action)
      logi.action <- match.arg(logi.action)
      nacoef.action <- match.arg(nacoef.action)
    }
    #' detect NA coefficients
    if(!all(is.finite(as.matrix(coef(object))))) {
      gripe <- "Cannot compute variance; some coefficients are NA, NaN or Inf"
      switch(nacoef.action,
             fatal = stop(gripe, call.=FALSE),
             warn = warning(gripe, call.=FALSE),
             silent = {})
      return(NULL)
    }
    #' extract stuff from fitted model
    Inter        <- object$Inter
    interaction  <- Inter$interaction
    itags        <- Inter$itags
    Vnamelist    <- object$Fit$Vnamelist
    Isoffsetlist <- object$Fit$Isoffsetlist
    glmdata      <- object$Fit$moadf
    fitter       <- object$Fit$fitter
    fitobj       <- object$Fit$FIT
    #' compute fitted intensity
    if(is.null(new.coef)) {
      fi <- fitted(fitobj)
    } else if(fitter != "glm") {
      warn.once("vcovMppmGAMnew", "'new.coef' is not supported by vcov.mppm for GAM or LME models; ignored")
      new.coef <- NULL
      fi <- fitted(fitobj)
    } else {
      fi <- GLMpredict(fitobj, glmdata, new.coef, changecoef=TRUE, type="response")
    }
    #' initialise
    cnames <- names(fixed.effects(object))
    nc <- length(cnames)
    A2 <- A3 <- matrix(0, nc, nc, dimnames=list(cnames, cnames))    
    #' (1) Compute matrix A1 directly
    glmsub  <- glmdata$.mpl.SUBSET
    wt      <- glmdata$.mpl.W
    mom <- model.matrix(object)
    lam <- unlist(fitted(object, new.coef=new.coef))
    A1 <- sumouter(mom, lam * wt * glmsub)
    #' (2) compute matrices A2 and A3 for submodels
    #' compute submodels 
    subs <- subfits(object, what="basicmodels", new.coef=new.coef)
    n <- length(subs)
    #' identify the (unique) active interaction in each row
    activeinter <- active.interactions(object)
    ## compute A2 and A3 for each submodel
    guts <- lapply(subs,
                   vcov,
                   what="internals",
                   matrix.action=matrix.action,
                   gam.action=gam.action,
                   logi.action=logi.action,
                   dropcoef=TRUE,
                   ...)
    a2   <- lapply(guts, getElement, name="A2")
    a3   <- lapply(guts, getElement, name="A3")
    #' (3) Determine map from interaction variables of subfits
    #'     to canonical variables of 'object'
    maps <- mapInterVars(object, subs, mom)
    #' (4) Process each row, summing A2 and A3
    for(i in seq_len(n)) {
      subi <- subs[[i]]
      cmap <- maps[[i]]
      #' contributes to second order terms only if non-Poisson
      if(!is.poisson(subi)) {
        cnames.i <- names(coef(subi))
        a2i <- a2[[i]]
        a3i <- a3[[i]]
        #' the (unique) tag name of the interaction in this model
        tagi <- colnames(activeinter)[activeinter[i,]]
        #' the corresponding canonical variable name(s) for this interaction
        vni <- Vnamelist[[tagi]]
        #' ignore offset variables
        iso <- Isoffsetlist[[tagi]]
        vni <- vni[!iso]
        if(length(vni)) {
          #' retain only interaction rows & columns (the rest are zero anyway)
          e <- cnames.i %in% vni
          a2ie <- a2i[e, e, drop=FALSE]
          a3ie <- a3i[e, e, drop=FALSE]
          #' all possible mappings 
          mappings <- do.call(expand.grid,
                              append(cmap, list(stringsAsFactors=FALSE)))
          nmappings <- nrow(mappings)
          if(nmappings == 0) {
            warning("Internal error: Unable to map submodel to full model")
          } else {
            for(irow in 1:nmappings) {
              for(jcol in 1:nmappings) {
                cmi <- as.character(mappings[irow,])
                cmj <- as.character(mappings[jcol,])
                if(anyDuplicated(cmi) || anyDuplicated(cmj)) {
                  warning("Internal error: duplicated labels in submodel map")
                } else if(!is.null(a2ie)) {
                  A2[cmi,cmj] <- A2[cmi,cmj] + a2ie
                  A3[cmi,cmj] <- A3[cmi,cmj] + a2ie
                }
              }
            }
          }
        }
      }
    }
    #' (5) pack up
    internals <- list(A1=A1, A2=A2, A3=A3)
    if(what %in% c("internals", "all"))
      internals <- c(internals, list(suff=mom))
    if(what %in% c("vcov", "corr", "all")) {
      #' variance-covariance matrix required
      U <- checksolve(A1, matrix.action, , "variance")
      vc <- if(is.null(U)) NULL else (U %*% (A1 + A2 + A3) %*% U)
    }
    out <- switch(what,
                  fisher = A1 + A2 + A3,
                  vcov   = vc,
                  corr   = {
                    if(is.null(vc)) return(NULL)
                    sd <- sqrt(diag(vc))
                    vc / outer(sd, sd, "*")
                  },
                  internals = internals,
                  all = list(internals=internals,
                             fisher=A1+A2+A3,
                             varcov=vc,
                             invgrad=A1)
                  )
    return(out)
  }

  addsubmatrix <- function(A, B, guessnames) {
    if(is.null(B)) return(A)
    if(is.null(colnames(B)) && !missing(guessnames)) {
      if(is.character(guessnames))
        guessnames <- list(guessnames, guessnames)
      if(all(lengths(guessnames) == dim(B)))
        colnames(B) <- guessnames
    }
    if(is.null(colnames(B))) {
      #' unusual
      if(!all(dim(A) == dim(B))) 
        stop("Internal error: no column names, and matrices non-conformable")
      A <- A + B
      return(A)
    }
    j <- match(colnames(B), colnames(A))
    if(anyNA(j)) 
      stop("Internal error: unmatched column name(s)")
    A[j,j] <- A[j,j] + B
    return(A)
  }

  bindsubmatrix <- function(A, B) {
    if(is.null(B)) return(A)
    if(is.null(colnames(B))) {
      if(ncol(A) != ncol(B))
        stop("Internal error: no column names, and matrices non-conformable")
      A <- rbind(A, B)
      return(A)
    }
    j <- match(colnames(B), colnames(A))
    if(anyNA(j))
      stop("Internal error: unmatched column name(s)")
    BB <- matrix(0, nrow(B), ncol(A))
    BB[,j] <- B
    A <- rbind(A, BB)
    return(A)
  }

  mergeAlternatives <- function(A, B) {
    okA <- !sapply(A, is.null)
    okB <- !sapply(B, is.null)
    if(any(override <- !okA & okB))
      A[override] <- B[override]
    return(A)
  }

##  notallzero <- function(df) { apply(df != 0, 2, any) }
  
  vcov.mppm
  
})
