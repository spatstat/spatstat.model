#'
#'       summary.kppm.R
#'
#'   $Revision: 1.38 $  $Date: 2022/10/23 07:37:37 $
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

  #' sibling probability
  if(object$isPCP) result$psib <- mean(psib(object))
  #' overdispersion index
  win <- as.owin(object, from="points")
  vac <- varcount(object, B=win)
  Lam <- integral(predict(object, window=win))
  result$odi <- vac/Lam
  #'
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
    switch(Fit$method,
           mincon = {
             splat(fittedby, "minimum contrast")
             splat("\tSummary statistic:", Fit$StatName)
             print(Fit$mcfit)
           },
           clik  =,
           clik2 = {
             splat(fittedby, "maximum second order composite likelihood")
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
    if(isPCP) {
      splat("Mean cluster size: ",
            if(!is.im(mu)) paste(signif(mu, digits), "points") else "[pixel image]")
    } else {
      splat("Fitted mean of log of random intensity:",
            if(!is.im(mu)) signif(mu, digits) else "[pixel image]")
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
  

  #'
  invisible(NULL)
}
