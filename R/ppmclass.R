#
#	ppmclass.R
#
#	Class 'ppm' representing fitted point process models.
#
#
#	$Revision: 2.159 $	$Date: 2024/10/06 10:32:04 $
#
#       An object of class 'ppm' contains the following:
#
#            $method           model-fitting method (currently "mpl")
#
#            $coef             vector of fitted regular parameters
#                              as given by coef(glm(....))
#
#            $trend            the trend formula
#                              or NULL 
#
#            $interaction      the interaction family 
#                              (an object of class 'interact') or NULL
#
#            $Q                the quadrature scheme used
#
#            $maxlogpl         the maximised value of log pseudolikelihood
#
#            $internal         list of internal calculation results
#
#            $correction       name of edge correction method used
#            $rbord            erosion distance for border correction (or NULL)
#
#            $the.call         the originating call to ppm()
#
#            $the.version      version of mpl() which yielded the fit
#
#
#------------------------------------------------------------------------

is.ppm <- function(x) { inherits(x, "ppm") }

print.ppm <-
function(x, ...,
         what=c("all", "model", "trend", "interaction", "se", "errors")) {

  verifyclass(x, "ppm")

  misswhat <- missing(what) 

  opts <- c("model", "trend", "interaction", "se", "errors")
  what <- match.arg(what, c("all", opts), several.ok=TRUE)
  if("all" %in% what) what <- opts

  np <- length(coef(x))
  terselevel <- spatstat.options("terse")
  digits <- getOption('digits')
  showextras <- waxlyrical('extras', terselevel)
  
  ## secret option
  showname <- !isFALSE(list(...)$showname)
  
  ## Determine whether SE is required 
  want.SE <- force.SE <- force.no.SE <- FALSE
  if(!misswhat && ("se" %in% what)) {
    ## SE was explicitly requested
    force.SE <- TRUE
  } else {
    ## Default rule: compute SE only if the model is Poisson
    switch(spatstat.options("print.ppm.SE"),
           always = { force.SE <- TRUE }, 
           never  = { force.no.SE <- TRUE },
           poisson = {
             want.SE <- is.poisson(x) && showextras
           })
  }
  do.SE <- (want.SE || force.SE) && !force.no.SE
  if(do.SE) {
    ## Check whether able to compute SE
    unable.SE <- (np == 0) || any(x$fitter %in% "gam") ||
      !is.null(x$internal$VB) || 
      (any(x$method %in% "mppm") && is.null(x$varcov))
    ## resolve
    if(force.SE && unable.SE) 
      warning("Unable to compute variances for this model", call.=FALSE)
    do.SE <- do.SE && !unable.SE
  }

  s <- summary.ppm(x, quick=if(do.SE) FALSE else "no variances")
        
  notrend <-    s$no.trend
#  stationary <- s$stationary
  poisson <-    s$poisson
  markeddata <- s$marked
  multitype  <- s$multitype
  dataname   <- s$dataname
  
#  markedpoisson <- poisson && markeddata
  csr <- poisson && notrend && !markeddata

  special <- csr && all(c("model", "trend") %in% what)
  if(special) {
    ## ---------- Trivial/special cases -----------------------
    splat("Stationary Poisson process")
    if(showname && showextras)
      splat("Fitted to point pattern dataset", sQuote(dataname))
    cat("Intensity:", signif(s$trend$value, digits), fill=TRUE)
  } else {
    ## ----------- Print model type -------------------
    if("model" %in% what) {
      splat(s$name)
      if(showname && showextras)
        splat("Fitted to point pattern dataset", sQuote(dataname))
      parbreak(terselevel)
        
      if(markeddata) mrk <- s$entries$marks
      if(multitype) {
        splat(paste("Possible marks:",
                    commasep(sQuote(levels(mrk)))))
        parbreak(terselevel)
      }
    }
    ## ----- trend --------------------------
    if("trend" %in% what) {
      if(!notrend) {
        splat("Log",
              if(poisson) "intensity: " else "trend: ",
              pasteFormula(s$trend$formula))
        parbreak(terselevel)
      }

      if(waxlyrical('space', terselevel) || !do.SE) {
        ## print trend coefficients, unless redundant and space is tight
        tv <- s$trend$value
      
        if(length(tv) == 0) 
          splat("[No trend coefficients]")
        else {
          thead <- paste0(s$trend$label, ":")
          if(is.list(tv)) {
            splat(thead)
            for(i in seq_along(tv))
              print(tv[[i]])
          } else if(is.numeric(tv) && length(tv) == 1) {
            ## single number: append to end of current line
            tvn <- names(tv)
            tveq <- if(is.null(tvn)) "\t" else paste(" ", tvn, "= ")
            splat(paste0(thead, tveq, signif(tv, digits)))
          } else {
            ## some other format 
            splat(thead)
            print(tv)
          }
        }
        parbreak(terselevel)
      }
    }

    if(waxlyrical("space", terselevel) &&
       !is.null(cfa <- s$covfunargs) && length(cfa) > 0) {
      cfafitter <- s$cfafitter
      if(is.null(cfafitter)) {
        cat("Covariate", "function", "arguments", "(covfunargs)",
            "provided:", fill=TRUE)
      } else {
        cat("Irregular", "parameters", "(covfunargs)",
            "fitted", "by", paste0(sQuote(cfafitter), ":"),
            fill=TRUE)
      }
      for(i in seq_along(cfa)) {
        cat(paste(names(cfa)[i], "= "))
        cfai <- cfa[[i]]
        if(is.numeric(cfai) && length(cfai) == 1) {
          cfai <- signif(cfai, digits)
          cat(paste(cfai, "\n"))
        } else print(cfai)
      }
    }
  }
  
  # ---- Interaction ----------------------------

  if("interaction" %in% what) {
    if(!poisson) {
      print(s$interaction, family=FALSE, banner=FALSE, brief=!showextras)
      parbreak(terselevel)
    }
  }
  
  # ----- parameter estimates with SE and 95% CI --------------------
  if(showextras && ("se" %in% what) && (np > 0)) {
    if(!is.null(cose <- s$coefs.SE.CI)) {
      print(cose, digits=digits)
    } else if(do.SE) {
      # standard error calculation failed
      splat("Standard errors unavailable; Fisher information matrix is singular")
    } else if(!force.no.SE) {
      # standard error was voluntarily omitted
      if(waxlyrical('space', terselevel))
        splat("For standard errors, type coef(summary(x))\n")
    }
  }
  
  # ---- Warnings issued in mpl.prepare  ---------------------

  if(waxlyrical("errors", terselevel) && "errors" %in% what) {
    probs <- s$problems
    if(!is.null(probs) && is.list(probs) && (length(probs) > 0)) 
      lapply(probs,
             function(x) {
               if(is.list(x) && !is.null(p <- x$print))
                 splat(paste("Problem:\n", p, "\n\n"))
             })
    
    if(s$old)
      warning(paste("Model fitted by old spatstat version", s$version))
        
  # ---- Algorithm status ----------------------------

    fitter <- s$fitter
    converged <- s$converged
    if(!is.null(fitter) && fitter %in% c("glm", "gam") && !converged)
      splat("*** Fitting algorithm for", sQuote(fitter),
            "did not converge ***")
  }

  if(showextras && s$projected) {
    parbreak()
    splat("Fit was emended to obtain a valid point process model")
  }

  if(identical(s$valid, FALSE) && waxlyrical("errors", terselevel)) {
    parbreak()
    splat("*** Model is not valid ***")
    if(!all(is.finite(s$entries$coef))) {
      splat("*** Some coefficients are NA or Inf ***")
    } else {
      splat("*** Interaction parameters are outside valid range ***")
    }
  } else if(showextras && is.na(s$valid)) {
    parbreak()
    splat("[Validity of model could not be checked]")
  }
  
  return(invisible(NULL))
}

# Extract version string from ppm object

versionstring.ppm <- function(object) {
  verifyclass(object, "ppm")
  v <- object$version
  if(is.null(v) || !is.list(v))
    v <- list(major=1, minor=3, release=4)
  vs <- paste(v$major, ".", v$minor, "-", v$release, sep="")
  return(vs)
}

# Extract quadrature scheme

quad.ppm <- function(object, drop=FALSE, clip=FALSE) {
  if(!is.ppm(object)) {
    if(is.kppm(object)) object <- object$po else
    if(is.lppm(object)) object <- object$fit else
    stop("object is not of class ppm, kppm or lppm")
  }
  Q <- object$Q
  if(is.null(Q))
    return(Q)
  if(drop || clip) {
    ok <- getglmsubset(object)
    if(!is.null(ok))
      Q <- Q[ok]
  }
  if(clip && object$correction == "border") {
    Wminus <- erosion(as.owin(object), object$rbord)
    Q <- Q[Wminus]
  }
  return(Q)
}

data.ppm <- function(object) { 
  verifyclass(object, "ppm")
  object$Q$data
}

dummy.ppm <- function(object, drop=FALSE) { 
  return(quad.ppm(object, drop=drop)$dummy)
}
  
# method for 'coef'
coef.ppm <- function(object, ...) {
  verifyclass(object, "ppm")
  object$coef
}

# extract internal fit

getglmdata.ppm <- function(object, ..., drop=FALSE) {
  verifyclass(object, "ppm")
  gd <- object$internal$glmdata
  if(!drop || is.null(gd)) return(gd)
  return(gd[getglmsubset(object), , drop=FALSE])
}

getglmfit.ppm <- function(object, ...) {
  verifyclass(object, "ppm")
  glmfit <- object$internal$glmfit
  if(is.null(glmfit))
      return(NULL)
  if(object$method != "mpl")
    glmfit$coefficients <- object$coef
  return(glmfit)
}

getglmsubset.ppm <- function(object, ...) {
  verifyclass(object, "ppm")
  gd <- object$internal$glmdata
  if(is.null(gd)) return(NULL)
  if(object$method=="logi") gd$.logi.ok else gd$.mpl.SUBSET
}

hasglmfit.ppm <- function(object) {
  verifyclass(object, "ppm")
  return(!is.null(object$internal$glmfit))
}

## internal, not exposed to user

getppmdatasubset <- function(object) {
  ## Equivalent to getglmsubset(object)[is.data(quad.ppm(object))]
  ## but also works for models fitted exactly, etc
  ##
  sub <- getglmsubset(object)
  if(!is.null(sub)) {
    Z <- is.data(quad.ppm(object))
    subZ <- sub[Z]
  } else {
    X <- data.ppm(object)
    subZ <- if(object$correction == "border") {
              (bdist.points(X) >= object$rbord)
            } else rep(TRUE, npoints(X))
  }
  return(subZ)
}


getppmOriginalCovariates <- function(object) {
  df <- as.data.frame(as.ppp(quad.ppm(object)))
  cova <- object$covariates
  if(length(cova) > 0) {
    df2 <- mpl.get.covariates(object$covariates,
                              union.quad(quad.ppm(object)),
                              "quadrature points",
                              object$covfunargs)
    df <- cbind(df, df2)
  } 
  return(df)
}
  
# ??? method for 'effects' ???

valid <- function(object, ...) {
  UseMethod("valid")
}

valid.ppm <- function(object, warn=TRUE, ...) {
  verifyclass(object, "ppm")
  coeffs <- coef(object)
  # ensure all coefficients are fitted, and finite
  if(!all(is.finite(coeffs)))
    return(FALSE)
  # inspect interaction
  inte <- object$interaction
  if(is.poisson(object))
    return(TRUE) # Poisson process
  # extract fitted interaction coefficients
  Vnames <- object$internal$Vnames
  IsOffset <- object$internal$IsOffset  
  Icoeffs <- coeffs[Vnames[!IsOffset]]
  # check interaction
  checker <- inte$valid
  if(is.null(checker) || !newstyle.coeff.handling(inte)) {
    if(warn) warning("Internal error: unable to check validity of model")
    return(NA)
  }
  #' remove prefix to obtain coefficient names expected by interaction
  if(npre <- sum(nchar(object$internal$vnameprefix)))
    names(Icoeffs) <- substring(names(Icoeffs), npre+1L)
  answer <- checker(Icoeffs, inte)
  return(answer)
}

emend <- function(object, ...) {
  UseMethod("emend")
}

emend.ppm <- project.ppm <- local({
  tracemessage <- function(depth, ...) {
    if(depth == 0) return(NULL)
    spacer <- paste(rep.int("  ", depth), collapse="")
    marker <- ngettext(depth, "trace", paste("trace", depth))
    marker <- paren(marker, "[")
    splat(paste0(spacer, marker, " ", paste(...)))
  }
  leaving <- function(depth) {
    tracemessage(depth, ngettext(depth, "Returning.", "Exiting level."))
  }
  emend.ppm <- function(object, ..., fatal=FALSE, trace=FALSE) {
    verifyclass(object, "ppm")
    fast <- spatstat.options("project.fast")
    # user specifies 'trace' as logical
    # but 'trace' can also be integer representing trace depth
    td <- as.integer(trace)
    trace <- (td > 0)
    tdnext <- if(trace) td+1 else 0
    if(valid.ppm(object)) {
      tracemessage(td, "Model is valid.")
      leaving(td)
      return(object)
    }
    # First ensure trend coefficients are all finite
    coeffs <- coef(object)
    # Which coefficients are trend coefficients
    coefnames  <- names(coeffs)
    internames <- object$internal$Vnames
    trendnames <- coefnames[!(coefnames %in% internames)]
    # Trend terms in trend formula
    trendterms <- attr(terms(object), "term.labels")
    # Mapping from coefficients to terms of GLM
    coef2term  <- attr(model.matrix(object), "assign")
    istrend <- (coef2term > 0) & (coefnames %in% trendnames)
    # Identify non-finite trend coefficients
    bad <- istrend & !is.finite(coeffs)
    if(!any(bad)) {
      tracemessage(td, "Trend terms are valid.")
    } else {
      nbad <- sum(bad)
      tracemessage(td,
                   "Non-finite ",
                   ngettext(nbad,
                            "coefficient for term ",
                            "coefficients for terms "),
                   commasep(sQuote(trendterms[coef2term[bad]])))
      if(fast) {
        # remove first illegal term
        firstbad <- min(which(bad))
        badterm <- trendterms[coef2term[firstbad]]
        # remove this term from model
        tracemessage(td, "Removing term ", sQuote(badterm))
        removebad <- as.formula(paste("~ . - ", badterm), env=object$callframe)
        newobject <- update(object, removebad)
        if(trace) {
          tracemessage(td, "Updated model:")
          print(newobject)
        }
        # recurse
        newobject <- emend.ppm(newobject, fatal=fatal, trace=tdnext)
        # return
        leaving(td)
        return(newobject)
      } else {
        # consider all illegal terms
        bestobject <- NULL
        for(i in which(bad)) {
          badterm <- trendterms[coef2term[i]]
          # remove this term from model
          tracemessage(td, "Considering removing term ", sQuote(badterm))
          removebad <- as.formula(paste("~ . - ", badterm),
                                  env=object$callframe)
          object.i <- update(object, removebad)
          if(trace) {
            tracemessage(td, "Considering updated model:")
            print(object.i)
          }
          # recurse
          object.i <- emend.ppm(object.i, fatal=fatal, trace=tdnext)
          # evaluate logPL
          logPL.i   <- logLik(object.i, warn=FALSE)
          tracemessage(td, "max log pseudolikelihood = ", logPL.i)
          # optimise
          if(is.null(bestobject) || (logLik(bestobject, warn=FALSE) < logPL.i))
            bestobject <- object.i
        }
        if(trace) {
          tracemessage(td, "Best submodel:")
          print(bestobject)
        }
        # return
        leaving(td)
        return(bestobject)
      }
    } 
    # Now handle interaction
    inte <- object$interaction
    if(is.null(inte)) {
      tracemessage(td, "No interaction to check.")
      leaving(td)
      return(object)
    }
    tracemessage(td, "Inspecting interaction terms.")
    proj <- inte$project
    if(is.null(proj)) {
      whinge <- "Internal error: interaction has no projection operator"
      if(fatal) stop(whinge) 
      warning(whinge)
      leaving(td)
      return(object)
    }
    # ensure the same edge correction is used!
    correction <- object$correction
    rbord      <- object$rbord
    # apply projection 
    coef.orig <- coeffs <- coef(object)
    Vnames   <- object$internal$Vnames
    Icoeffs  <- coeffs[Vnames]
    change <- proj(Icoeffs, inte)
    if(is.null(change)) {
      tracemessage(td, "Interaction does not need updating.")
      leaving(td)
      return(object)
    }
    tracemessage(td, "Interaction is not valid.")
    if(is.numeric(change)) {
      tracemessage(td, "Interaction coefficients updated without re-fitting.")
      # old style: 'project' returned a vector of updated coefficients
      Icoeffs <- change
      # tweak interaction coefficients
      object$coef[Vnames] <- Icoeffs
      # recompute fitted interaction
      object$fitin <- NULL
      object$fitin <- fitin(object)
    } else if(is.interact(change)) {
      # new style: 'project' returns an interaction
      if(trace) {
        tracemessage(td, "Interaction changed to:")
        print(change)
      }
      # refit the whole model 
      #      (using the same edge correction)
      #      (and the same quadrature scheme)
      newobject <- update(object, interaction=change,
                          correction=correction, rbord=rbord,
                          forcefit=TRUE,
                          envir=object$callframe)
      if(trace) {
        tracemessage(td, "Updated model:")
        print(newobject)
      }
      # recurse
      newobject <- emend.ppm(newobject, fatal=fatal, trace=tdnext)
      object <- newobject
    } else if(is.list(change) && all(unlist(lapply(change, is.interact)))) {
      # new style: 'project' returns a list of candidate interactions
      nchange <- length(change)
      tracemessage(td, "Considering", nchange,
                   ngettext(nchange, "submodel", "submodels"))
      bestobject <- NULL
      for(i in seq_len(nchange)) {
        change.i <- change[[i]]
        if(trace) {
          tracemessage(td,
                       "Considering", ordinal(i), 
                       "candidate submodel, with interaction:")
          print(change.i)
        }
        # refit the whole model
        object.i <- update(object, interaction=change.i,
                           correction=correction, rbord=rbord,
                           forcefit=TRUE,
                           envir=object$callframe)
        if(trace) {
          tracemessage(td, "Considering", ordinal(i),
                       "candidate updated model:")
          print(object.i)
        }
        # recurse
        object.i <- emend.ppm(object.i, fatal=fatal, trace=tdnext)
        # evaluate logPL
        logPL.i   <- logLik(object.i, warn=FALSE)
        tracemessage(td, "max log pseudolikelihood = ", logPL.i)
        # optimise
        if(is.null(bestobject) || (logLik(bestobject, warn=FALSE) < logPL.i))
          bestobject <- object.i
      }
      # end loop through submodels
      if(trace) {
        tracemessage(td, "Best submodel:")
        print(bestobject)
      }
      object <- bestobject
    } else stop("Internal error: unrecognised format of update")
    object$projected <- TRUE
    object$coef.orig  <- coef.orig
    leaving(td)
    return(object)
  }
  emend.ppm
})

# more methods

deviance.ppm <- function(object, ...) {
  force(object)
  satlogpl <- object$satlogpl
  if(is.null(satlogpl)) {
    object <- update(object, forcefit=TRUE)
    satlogpl <- object$satlogpl
  }
  if(is.null(satlogpl) || !is.finite(satlogpl))
    return(NA)
  ll <- do.call(logLik,
                resolve.defaults(list(object=quote(object), absolute=FALSE),
                                 list(...)))
  ll <- as.numeric(ll)
  2 * (satlogpl - ll)
}

logLik.ppm <- function(object, ..., new.coef=NULL, warn=TRUE, absolute=FALSE) {
  if(!is.poisson.ppm(object) && warn) 
    warn.once("ppmLogLik",
              "log likelihood is not available for non-Poisson model;",
              "log pseudolikelihood returned")
  ## degrees of freedom
  nip <- if(!inherits(object, "ippm")) 0 else
           length(attr(object$covfunargs, "free"))
  df <- length(coef(object)) + nip
  ## compute adjustment constant
  if(absolute && object$method %in% c("exact", "mpl", "ho")) {
    X <- data.ppm(object)
    W <- Window(X)
    areaW <-
      if(object$correction == "border" && object$rbord > 0) 
      eroded.areas(W, object$rbord) else area(W)
    constant <- areaW * markspace.integral(X)
  } else constant <- 0
  ##
  if(is.null(new.coef)) {
    ## extract from object
    ll <- object$maxlogpl + constant
    attr(ll, "df") <- df
    class(ll) <- "logLik"
    return(ll)
  } 
  ## recompute for new parameter values
  method <- object$method
  if(method == "exact")
    method <- update(method, forcefit=TRUE)
  Q <- quad.ppm(object, drop=TRUE)
  Z <- is.data(Q)
  cif <- fitted(object, type="cif", new.coef=new.coef, drop=TRUE)
  cifdata <- cif[Z]
  switch(method,
         mpl=,
         exact=,
         ho = {
           w <- w.quad(Q)
           ll <- sum(log(cifdata[cifdata > 0])) - sum(w * cif)
         },
         logi=,
         VBlogi={
           B <- getglmdata(object, drop=TRUE)$.logi.B
           p <- cif/(B+cif)
           ll <- sum(log(p/(1-p))[Z]) + sum(log(1-p)) + sum(log(B[Z]))
         },
         stop(paste("Internal error: unrecognised ppm method:",
                    dQuote(method)))
         )
  ll <- ll + constant
  attr(ll, "df") <- df
  class(ll) <- "logLik"
  return(ll)
}

pseudoR2 <- function(object, ...) {
  UseMethod("pseudoR2")
}

pseudoR2.slrm <- pseudoR2.ppm <- function(object, ..., keepoffset=TRUE) {
  dres <- deviance(object, ..., warn=FALSE)
  nullfmla <- . ~ 1
  if(keepoffset && has.offset.term(object)) {
    off <- attr(model.depends(object), "offset")
    offterms <- row.names(off)[apply(off, 1, any)]
    if(length(offterms)) {
      nullrhs <- paste(offterms, collapse=" + ") 
      nullfmla <- as.formula(paste(". ~ ", nullrhs))
    }
  } 
  nullmod <- update(object, nullfmla, forcefit=TRUE)
  dnul <- deviance(nullmod, warn=FALSE)
  return(1 - dres/dnul)
}

formula.ppm <- function(x, ...) {
  return(x$trend)
}

terms.ppm <- function(x, ...) {
  terms(x$terms, ...)
}

labels.ppm <- function(object, ...) {
  # extract fitted trend coefficients
  co <- coef(object)
  Vnames <- object$internal$Vnames
  is.trend <- !(names(co) %in% Vnames)
  # model terms
  tt <- terms(object)
  lab <- attr(tt, "term.labels")
  if(length(lab) == 0)
    return(character(0))
  # model matrix
  mm <- model.matrix(object)
  ass <- attr(mm, "assign")
  # 'ass' associates coefficients with model terms
  # except ass == 0 for the Intercept
  coef.ok <- is.finite(co)
  relevant <- (ass > 0) & is.trend
  okterms <- unique(ass[coef.ok & relevant])
  return(lab[okterms])
}

AIC.ppm <- function(object, ..., k=2, takeuchi=TRUE) {
  ll <- logLik(object, warn=FALSE)
  pen <- attr(ll, "df")
  if(takeuchi && !is.poisson(object)) {
    vv <- vcov(object, what="internals")
    logi <- (object$method == "logi")
    J  <- with(vv, if(!logi) Sigma else (Sigma1log+Sigma2log))
    H  <- with(vv, if(!logi) A1 else Slog)
    ## Takeuchi penalty = trace of J H^{-1} = trace of H^{-1} J
    JiH <- try(solve(H, J), silent=TRUE)
    if(!inherits(JiH, "try-error")) 
      pen <- sum(diag(JiH))
  } 
  return(- 2 * as.numeric(ll) + k * pen)
}

extractAIC.ppm <- function (fit, scale = 0, k = 2, ..., takeuchi=TRUE)
{
  edf <- length(coef(fit))
  aic <- AIC(fit, k=k, takeuchi=takeuchi)
  c(edf, aic)
}

#
# method for model.frame

model.frame.ppm <- function(formula, ...) {
  object <- formula
  gf <- getglmfit(object)
  if(is.null(gf)) {
    warning("Model re-fitted with forcefit=TRUE")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
  }
  argh <- resolve.defaults(list(formula=quote(gf)),
                           list(...),
                           list(data = getglmdata(object),
                                subset = TRUE))
  result <- switch(object$fitter,
                   gam = do.call(modelFrameGam, argh),
                   do.call(model.frame, argh))
  return(result)
}

#' a hacked version of model.frame.glm that works for gam objects (mgcv)
modelFrameGam <- function(formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
                      0L)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
#    fcall$method <- "model.frame"
    fcall[[1L]] <- quote(mgcv::gam)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) 
      env <- parent.frame()
    refut <- eval(fcall, env)
    refut$model
  } else formula$model
}

#
# method for model.matrix

model.matrix.ppm <- function(object,
                             data=model.frame(object, na.action=NULL),
                             ..., Q=NULL, keepNA=TRUE) {
  if(missing(data)) data <- NULL
  PPMmodelmatrix(object, data=data, ..., Q=Q, keepNA=keepNA)
}

model.matrix.ippm <- function(object,
                              data=model.frame(object, na.action=NULL),
                              ..., Q=NULL, keepNA=TRUE, irregular=FALSE) {
  if(missing(data)) data <- NULL 
  PPMmodelmatrix(object, data=data, ...,
                 Q=Q, keepNA=keepNA, irregular=irregular)
}

PPMmodelmatrix <- function(object,
                           data = NULL, 
                           ...,
                           subset, Q=NULL, keepNA=TRUE, irregular=FALSE,
                           splitInf=FALSE,
                           eps=NULL, dimyx=NULL, xy=NULL, rule.eps=NULL) {
  ## trap pixel resolution arguments - not used 
  dont.complain.about(eps, dimyx, xy, rule.eps)
  # handles ppm and ippm			      
  data.given <- !is.null(data)
  irregular <- irregular && inherits(object, "ippm") && !is.null(object$iScore)
  if(splitInf && !data.given && is.null(Q)) {
    #' force re-computation
    Q <- quad.ppm(object)
  }
  if(!is.null(Q)) {
    if(data.given) stop("Arguments Q and data are incompatible")
    if(!inherits(Q, c("ppp", "quad")))
      stop("Q should be a point pattern or quadrature scheme")
    if(is.ppp(Q)) Q <- quad(Q, Q[FALSE])
    ## construct Berman-Turner frame
    needed <- c("trend", "interaction", "covariates", "covfunargs",
                "correction", "rbord")
    bt <- do.call(bt.frame,
                  c(list(Q), object[needed], list(splitInf=splitInf)))
    forbid <- bt$forbid
    ## compute model matrix (including NA's if present)
    mf <- model.frame(bt$fmla, bt$glmdata, ..., na.action=NULL)
    mm <- model.matrix(bt$fmla, mf, ...)
    ass <- attr(mm, "assign")
    if(irregular) {
      ## add irregular score components
      U <- union.quad(Q)
      mi <- sapply(object$iScore, do.call,
                   args=append(list(x=U$x, y=U$y), object$covfunargs),
                   envir=environment(terms(object)))
      if(nrow(mi) != nrow(mm))
        stop("Internal error: incorrect number of rows in iScore")
      mm <- cbind(mm, mi)
    }
    ## subset
    if(!missing(subset)) {
      ok <- eval(substitute(subset), envir=bt$glmdata)
      mm <- mm[ok, , drop=FALSE]
      if(!is.null(forbid)) forbid <- forbid[ok]
    }
    ## remove NA's ?
    if(!keepNA) {
      ok <- complete.cases(mm)
      mm <- mm[ok, , drop=FALSE]
      if(!is.null(forbid)) forbid <- forbid[ok]
    }
    attr(mm, "assign") <- ass 
    attr(mm, "-Inf") <- forbid
    return(mm)
  }

  #' extract GLM fit 
  gf <- getglmfit(object)
  if(is.null(gf)) {
    warning("Model re-fitted with forcefit=TRUE")
    object <- update(object, forcefit=TRUE)
    gf <- getglmfit(object)
    if(is.null(gf))
      stop("internal error: unable to extract a glm fit")
  }

  if(data.given) {
    #' new data. Must contain the Berman-Turner variables as well.
    bt <- list(.mpl.Y=1, .mpl.W=1, .mpl.SUBSET=TRUE)
    if(any(forgot <- !(names(bt) %in% names(data)))) 
      data <- do.call(cbind, append(list(data), bt[forgot]))
    mm <- model.matrix(gf, data=data, ..., subset=NULL)
    ass <- attr(mm, "assign")
    if(irregular) {
      ## add irregular score components 
      mi <- sapply(object$iScore, do.call,
                   args=append(list(x=data$x, y=data$y), object$covfunargs),
                   envir=environment(terms(object)))
      if(nrow(mi) != nrow(mm))
        stop("Internal error: incorrect number of rows in iScore")
      mm <- cbind(mm, mi)
      attr(mm, "assign") <- ass
    }
    if(inherits(gf, "gam")) 
      attr(mm, "assign") <- gf$assign
    return(mm)
  }

  scrambled <- object$scrambled %orifnull% FALSE
  ## if TRUE, this object was produced by 'subfits' using jittered covariate

  if(!keepNA && !irregular && !scrambled) {
    # extract model matrix of glm fit object
    # restricting to its 'subset' 
    mm <- model.matrix(gf, ...)
    if(inherits(gf, "gam")) 
      attr(mm, "assign") <- gf$assign
    return(mm)
  }
  
  ## extract model matrix for all cases
  gd <- getglmdata(object, drop=FALSE) 
  if(!scrambled) {
    ## 'gf' was fitted to correct data. Use internals.
    mm <- model.matrix(gf, ..., subset=NULL, na.action=NULL)
    ass <- attr(mm, "assign")
  } else {
    ## 'gf' was originally fitted using jittered data:
    ## Use correct data given by 'gd'
    ## Temporarily add scrambled data to avoid singular matrices etc
    gds <- object$internal$glmdata.scrambled
    gdplus <- rbind(gd, gds)
    mm <- model.matrix(gf, ..., data=gdplus, subset=NULL, na.action=NULL)
    ass <- attr(mm, "assign")
    ## Now remove rows corresponding to scrambled data
    mm <- mm[seq_len(nrow(gd)), , drop=FALSE]
    attr(mm, "assign") <- ass
  } 
  cn <- colnames(mm)
  if(nrow(mm) != nrow(gd)) {
    # can occur if covariates include NA's or interaction is -Inf
    insubset <- getglmsubset(object)
    isna <- is.na(insubset) | !insubset
    if(sum(isna) + nrow(mm) == nrow(gd)) {
      # insert rows of NA's
      mmplus <- matrix( , nrow(gd), ncol(mm))
      mmplus[isna, ] <- NA
      mmplus[!isna, ] <- mm
      mm <- mmplus
      attr(mm, "assign") <- ass
    } else 
    stop("internal error: model matrix does not match glm data frame")
  }
  if(irregular) {
    ## add irregular score components 
    U <- union.quad(quad.ppm(object, drop=FALSE))
    mi <- sapply(object$iScore, do.call,
                 args=append(list(x=U$x, y=U$y), object$covfunargs),
	  envir=environment(terms(object)))
    if(nrow(mi) != nrow(mm))
      stop("Internal error: incorrect number of rows in iScore")
    mm <- cbind(mm, mi)
    attr(mm, "assign") <- ass
    cn <- c(cn, colnames(mi))
  }
  ## subset
  if(!missing(subset)) {
    ok <- eval(substitute(subset), envir=gd)
    mm <- mm[ok, , drop=FALSE]
    attr(mm, "assign") <- ass
  }
  ## remove NA's
  if(!keepNA) {
    mm <- mm[complete.cases(mm), , drop=FALSE]
    attr(mm, "assign") <- ass
  }
  if(inherits(gf, "gam")) 
    attr(mm, "assign") <- gf$assign
  colnames(mm) <- cn
  return(mm)
}

model.images <- function(object, ...) {
  UseMethod("model.images")
}

model.images.ppm <- function(object, W=as.owin(object), ...) {
  X <- data.ppm(object)
#  irregular <- resolve.1.default(list(irregular=FALSE), list(...))
  ## make a quadscheme with a dummy point at every pixel
  Q <- pixelquad(X, W, ...)  # recognises dimyx etc
  ## compute model matrix
  mm <- model.matrix(object, Q=Q, ...) # ignores dimyx etc
  ## retain only the entries for dummy points (pixels)
  mm <- mm[!is.data(Q), , drop=FALSE]
  mm <- as.data.frame(mm)
  ## create template image
  Z <- as.im(attr(Q, "M"))
  ok <- !is.na(Z$v)
  ## make images
  imagenames <- colnames(mm)
  if(!is.multitype(object)) {
    result <- lapply(as.list(mm), replace, list=ok, x=Z)
    result <- as.solist(result)
    names(result) <- imagenames
  } else {
    marx <- marks(Q$dummy)
    mmsplit <- split(mm, marx)
    result <- vector(mode="list", length=length(mmsplit))
    for(i in seq_along(mmsplit))
      result[[i]] <- as.solist(lapply(as.list(mmsplit[[i]]),
                                      replace, list=ok, x=Z))
    names(result) <- names(mmsplit)
    result <- do.call(hyperframe, result)
    row.names(result) <- imagenames
  }
  return(result)
}

unitname.ppm <- function(x) {
  return(unitname(x$Q))
}

"unitname<-.ppm" <- function(x, value) {
  unitname(x$Q) <- value
  return(x)
}

nobs.ppm <- function(object, ...) { npoints(data.ppm(object)) }

as.interact.ppm <- function(object) {
 verifyclass(object, "ppm")
 inte <- object$interaction
 if(is.null(inte))
   inte <- Poisson()
 return(inte)
}

as.ppm <- function(object) {
  UseMethod("as.ppm")
}

as.ppm.ppm <- function(object) {
  object
}

## method for as.owin

as.owin.ppm <- function(W, ..., from=c("points", "covariates"), fatal=TRUE) {
  if(!verifyclass(W, "ppm", fatal=fatal))
    return(NULL)
  from <- match.arg(from)
  datawin <- as.owin(data.ppm(W))
  if(from == "points")
    return(datawin)
  covs <- W$covariates
  isim <- unlist(lapply(covs, is.im))
  if(!any(isim))
    return(datawin)
  cwins <- lapply(covs[isim], as.owin)
  covwin <- do.call(intersect.owin, unname(cwins))
  result <- intersect.owin(covwin, datawin)
  return(result)
}

domain.ppm <- Window.ppm <- function(X, ..., from=c("points", "covariates")) {
  from <- match.arg(from)
  as.owin(X, ..., from=from)
}

hardcoredist.ppm <- function(x, ..., epsilon=0) {
  hardcoredist.fii(fitin(x), ..., epsilon=epsilon)
}

## change the coefficients in a ppm or other model

tweak.coefs <- function(model, new.coef) {
  if(is.null(new.coef)) return(model)
  co <- coef(model)
  check.nvector(new.coef, length(co), things="coefficients", vname="new.coef")
  model$coef.orig <- co
  model$coef <- new.coef
  return(model)
}

spatialCovariateUnderModel <- function(model, covariate, ...) {
  UseMethod("spatialCovariateUnderModel")
}

spatialCovariateUnderModel.ppm <-
spatialCovariateUnderModel.kppm <-
spatialCovariateUnderModel.dppm <- function(model, covariate, ...) {
  Q <- quad.ppm(as.ppm(model))
  loc <- as.ppp(Q)
  df <- mpl.get.covariates(list(Z=covariate), loc, covfunargs=list(...))
  df$wt <- fitted(model) * w.quad(Q)
  return(df)
}

