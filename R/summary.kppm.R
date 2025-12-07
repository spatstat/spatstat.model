#'
#'       summary.kppm.R
#'
#'   $Revision: 1.43 $  $Date: 2025/12/07 02:25:28 $
#' 

summary.kppm <- function(object, ..., quick=FALSE) {
  nama <- names(object)
  result <- unclass(object)[!(nama %in% c("X", "po", "call", "callframe"))]
  ## handle old format
  if(is.null(result$isPCP)) result$isPCP <- TRUE
  ## extract 'optim' object
  Fit <- object$Fit
  opt <- switch(Fit$method,
                mincon = Fit$mcfit$opt,
                clik  =,
                clik2 = Fit$clfit,
                palm = Fit$clfit,
                waag = Fit$clfit,
                adapcl = Fit$cladapfit,
                warning(paste("Unrecognised fitting method",
                              sQuote(Fit$method)))
                )
  if(Fit$method != "adapcl") {
    result$optim.converged <- optimConverged(opt)
    result$optim.status    <- optimStatus(opt)
    result$optim.nsteps    <- optimNsteps(opt)
  }
  ## summarise trend component
  result$trend <- summary(as.ppm(object), ..., quick=quick)
  if(isFALSE(quick)) {
    theta <- coef(object)
    if(length(theta) > 0) {
      vc <- vcov(object, matrix.action="warn")
      if(!is.null(vc)) {
        se <- if(is.matrix(vc)) sqrt(diag(vc)) else
        if(length(vc) == 1) sqrt(vc) else NULL
      }
      if(!is.null(se)) {
        two <- qnorm(0.975)
        lo <- theta - two * se
        hi <- theta + two * se
        zval <- theta/se
        pval <- 2 * pnorm(abs(zval), lower.tail=FALSE)
        psig <- cut(pval, c(0,0.001, 0.01, 0.05, 1),
                    labels=c("***", "**", "*", "  "),
                    include.lowest=TRUE)
        ## table of coefficient estimates with SE and 95% CI
        result$coefs.SE.CI <- data.frame(Estimate=theta, S.E.=se,
                                         CI95.lo=lo, CI95.hi=hi,
                                         Ztest=psig,
                                         Zval=zval)
      }
    }
  }

  if(object$isPCP) {
    #' ---------- information from cluster parameters paper -------------
    #' sibling probability
    result$psib <- mean(psib(object))
    #' overdispersion index
    win <- as.owin(object, from="points")
    vac <- varcount(object, B=win)
    Lam <- integral(predict(object, window=win))
    result$odi <- vac/Lam
    #' detect penalised fit
    result$penalised <- !is.null(Fit$pspace.used$penalty)
    #' optimization trace
    result$trace <- attr(object, "h")
    #' spatial persistence (over window)
    g <- pcfmodel(object)
    d <- diameter(win)
    result$persist <- (g(d)-1)/(g(0)-1)
    result$ispo <- poisson.fits.better(object)
    #' bounds on distance from Poisson and mixed Poisson
    aW <- area(win)
    if(is.stationary(object)) {
      lambda <- object$lambda
      mu <- object$mu
      EN <- lambda * aW
    } else {
      EN <- Lam # integral of intensity
      mu <- max(object$mu)
    }
    #' first bound
    tvbound1 <- 2 * EN * (1-exp(-mu))
    rules <- spatstatClusterModelInfo(object$clusters)
    newpar <- object$clustpar
    oldpar <- rules$checkpar(newpar, native=TRUE, strict=FALSE)
    if(all(oldpar > 0)) {
      scal <- newpar[["scale"]]
      kappa <- newpar[["kappa"]]
      result$phi <- phinew <- g(0) - 1
      A10 <- phinew * kappa * scal^2
      h10 <- rules$kernel(oldpar, 0)
      #' second and third bounds
      tvbound2 <- phinew * EN^2 * (1 + scal^2/(aW * A10))
      tvbound3 <- phinew * EN^2 * (1 + h10/A10)
      #' fourth (new) bound
      tvbound4 <- EN * sqrt(phinew)
      #' save bounds
      result$tvbound1 <- tvbound1
      result$tvbound2 <- tvbound2
      result$tvbound3 <- tvbound3
      result$tvbound4 <- tvbound4
      result$tvbound <- min(1, tvbound1, tvbound2, tvbound3, tvbound4)
      if(is.stationary(object)) {
        #' characteristics of nonempty clusters
        em <- exp(-mu)
        result$kappa1 <- kappa * (1-em)
        result$mu1 <- mu/(1-em)
        result$kappa2 <- kappa * (1 - em - mu * em)
        result$eta <- kappa/(lambda + kappa)
        result$panysib <- 1-em
        ## distance to mixed Poisson
        A1d <- (g(d)-1) * kappa * scal^2
        tvbmix <- EN * (phinew/A10) * sqrt(2 * (A10 - A1d))
        result$tvbmix <- min(1, tvbmix)
      }
    }
    ## -----------------------------------------------------
  }
  
  class(result) <- "summary.kppm"
  return(result)
}

coef.summary.kppm <- function(object, ...) {
  return(object$coefs.SE.CI)
}

print.summary.kppm <- function(x, ...) {
  terselevel <- spatstat.options('terse')
  digits <- getOption('digits')
  isPCP <- x$isPCP
  splat(if(x$stationary) "Stationary" else "Inhomogeneous",
        if(isPCP) "cluster" else "Cox",
        "point process model")

  if(waxlyrical('extras', terselevel) && nchar(x$Xname) < 20)
    splat("Fitted to point pattern dataset", sQuote(x$Xname))

  Fit <- x$Fit
  
  if(waxlyrical('gory', terselevel)) {
    fittedby <- "Fitted by"
    #' detect whether fit used a penalty
    if(isTRUE(x$penalised))
      fittedby <- "Fitted by penalised"
    switch(Fit$method,
           mincon = {
             splat(fittedby, "minimum contrast")
             splat("\tSummary statistic:", Fit$StatName)
             print(Fit$mcfit)
           },
           clik  =,
           clik2 = {
             splat(fittedby, "maximum Guan's second order composite likelihood")
             splat("\trmax =", Fit$rmax)
             if(!is.null(wtf <- Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
             printStatus(x$optim.status)
           },
           palm = {
             splat(fittedby, "maximum Palm likelihood")
             splat("\trmax =", Fit$rmax)
             if(!is.null(wtf <- Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
             printStatus(x$optim.status)
           },
           waag = {
             splat(fittedby, "maximum Waagepetersen's second order composite likelihood")
             splat("\trmax =", Fit$rmax)
             if(!is.null(wtf <- Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
             printStatus(x$optim.status)
           },
           adapcl = {
             splat("Fitted by adaptive second order composite likelihood")
             splat("\tepsilon =", x$Fit$epsilon)
             if(!is.null(wtf <- x$Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
           },
           warning(paste("Unrecognised fitting method", sQuote(Fit$method)))
           )
    #' optimization trace
    if(!is.null(x$trace)) {
      parbreak()
      splat("[Includes history of evaluations of objective function]")
    }
  }

  # ............... trend .........................

  parbreak()
  splat("----------- TREND  -----")
  print(x$trend, ...)

  # ..................... clusters ................

  tableentry <- spatstatClusterModelInfo(x$clusters)
  
  parbreak()
  splat("-----------", 
        if(isPCP) "CLUSTER" else "COX",
        "",
        "-----------")
  splat("Model:", tableentry$printmodelname(x))
  parbreak()
  
  cm <- x$covmodel
  if(!isPCP) {
    # Covariance model - LGCP only
    splat("\tCovariance model:", cm$model)
    margs <- cm$margs
    if(!is.null(margs)) {
      nama <- names(margs)
      tags <- ifelse(nzchar(nama), paste(nama, "="), "")
      tagvalue <- paste(tags, margs)
      splat("\tCovariance parameters:",
            paste(tagvalue, collapse=", "))
    }
  }
  pa <- x$clustpar
  if (!is.null(pa)) {
    splat("Fitted",
          if(isPCP) "cluster" else "covariance",
          "parameters:")
    print(pa, digits=digits)
  }

  if(!is.null(mu <- x$mu)) {
    moo <- signif(mean(mu), digits)
    if(is.im(moo)) moo <- paste0("[pixel image] [mean = ", moo, "]")
    if(isPCP) {
      splat("Mean cluster size: ", moo, if(!is.im(mu)) "points" else NULL)
    } else {
      splat("Fitted mean of log of random intensity:", moo)
    }
  }
  #' table of coefficient estimates with SE and 95% CI
  if(!is.null(cose <- x$coefs.SE.CI)) {
    parbreak()
    splat("Final standard error and CI")
    splat("(allowing for correlation of",
          if(isPCP) "cluster" else "Cox",
          "process):")
    print(cose)
  }

  #' Cluster strength indices
  psi <- x$psib
  odi <- x$odi
  if(!is.null(psi) || !is.null(odi)) {
    parbreak()
    splat("----------- cluster strength indices ---------- ")
    if(!is.null(psi)) {
      psi <- signif(psi, digits)
      if(isTRUE(x$stationary)) {
        splat("Sibling probability", psi)
      } else splat("Mean sibling probability", psi)
    }
    if(!is.null(odi))
      splat("Count overdispersion index (on original window):",
            signif(odi, digits))
  }
  
  #' cluster strength 
  if(!is.null(phi <- x$phi))
      splat("Cluster strength:", signif(phi, digits))
  #' spatial persistence (over window)
  if(!is.null(x$persist)) {
    parbreak()
    splat("Spatial persistence index (over window):", signif(x$persist, digits))
  }
  #' bound on total-variation distance from Poisson
  if(!is.null(x$tvbound)) {
    parbreak()
    splat("Bound on distance from Poisson process (over window):",
          signif(x$tvbound, digits))
    if(!is.null(x$tvbound1) && !is.null(x$tvbound2)
       && !is.null(x$tvbound3) && !is.null(x$tvbound4)) {
      tvb <- as.numeric(x[c("tvbound1", "tvbound2", "tvbound3", "tvbound4")])
      splat("\t = min", paren(paste(c(1, signif(tvb, digits)), collapse=", ")))
    }
  }
  if(!is.null(x$tvbmix)) {
    parbreak()
    splat("Bound on distance from MIXED Poisson process (over window):",
          signif(x$tvbmix, digits))
  }
  #'
  if(!is.null(x$kappa1)) {
    parbreak()
    splat("Intensity of parents of nonempty clusters:",
          signif(x$kappa1, digits))
    splat("Mean number of offspring in a nonempty cluster:",
          signif(x$mu1, digits))
    splat("Intensity of parents of clusters of more than one offspring point:",
          signif(x$kappa2, digits))
    splat("Ratio of parents to parents-plus-offspring:", signif(x$eta, digits),
          "(where 1 = Poisson process)")
    splat("Probability that a typical point belongs to a nontrivial cluster:",
          signif(x$panysib, digits))
  }
  if(isTRUE(x$ispo)) {
    parbreak()
    splat(">>> The Poisson process is a better fit <<< ")
  }
  #'
  invisible(NULL)
}
