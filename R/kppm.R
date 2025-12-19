#'
#'   kppm.R
#'
#'  kluster/kox point process models
#'
#'  $Revision: 1.240 $ $Date: 2025/12/19 02:52:53 $
#'
#'  Copyright (c) 2001-2025 Adrian Baddeley, Rolf Turner, Ege Rubak,
#'                Abdollah Jalilian and Rasmus Plenge Waagepetersen
#'  GNU Public Licence (>= 2.0)

kppm <- function(X, ...) {
  UseMethod("kppm")
}

kppm.formula <-
  function(X, clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
           ..., data=NULL) {
  ## remember call
  callstring <- short.deparse(sys.call())
  cl <- match.call()

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(X, "formula"))
    stop(paste("Argument 'X' should be a formula"))
  formula <- X
  
  if(spatstat.options("expand.polynom"))
    formula <- expand.polynom(formula)

  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"))
  Yexpr <- formula[[2L]]
  trend <- formula[c(1L,3L)]
  
  ## FIT #######################################
  thecall <- call("kppm", X=Yexpr, trend=trend,
                  data=data, clusters=clusters)
  ncall <- length(thecall)
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    thecall[ncall + 1:nargh] <- argh
    names(thecall)[ncall + 1:nargh] <- names(argh)
  }
##  result <- eval(thecall, 
##                 envir=if(!is.null(data)) data else parent.frame(),
##                 enclos=if(!is.null(data)) parent.frame() else baseenv())
  callenv <- list2env(as.list(data), parent=parent.frame())
  result <- eval(thecall, envir=callenv, enclos=baseenv())

  result$call <- cl
  result$callframe <- parent.frame()
  if(!("callstring" %in% names(list(...))))
    result$callstring <- callstring
  
  return(result)
}

kppm.ppp <- kppm.quad <-
  function(X, trend = ~1,
           clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
           data=NULL,
           ...,
           covariates = data,
           subset, 
           method = c("mincon", "clik2", "palm", "waag", "adapcl"),
           penalised = FALSE,
           improve.type = c("none", "clik1", "wclik1", "quasi"),
           improve.args = list(),
           weightfun=NULL,
           control=list(),
           stabilize=TRUE,
           algorithm,
           trajectory=FALSE,
           statistic="K",
           statargs=list(),
           rmax = NULL,
           epsilon=0.01,
           covfunargs=NULL,
           use.gam=FALSE,
           nd=NULL, eps=NULL,
           ppm.improve.type = c("none", "ho", "enet"),
           ppm.improve.args=list()) {
  cl <- match.call()
  callstring <- paste(short.deparse(sys.call()), collapse="")
  Xname <- short.deparse(substitute(X))
  clusters <- match.arg(clusters)
  improve.type <- match.arg(improve.type)
  method <- match.arg(method)
  if(method == "mincon")
    statistic <- pickoption("summary statistic", statistic,
                            c(K="K", g="pcf", pcf="pcf"))
  
  if(missing(algorithm)) {
    algorithm <- if(method == "adapcl") "Broyden" else "Nelder-Mead"
  } else check.1.string(algorithm)
  
  if(isTRUE(trajectory) && method == "adapcl")
    warning("trajectory=TRUE is not supported for method 'adapcl'", call.=FALSE)

  ClusterArgs <- list(method = method,
                      improve.type = improve.type,
                      improve.args = improve.args,
                      weightfun=weightfun,
                      control=control,
                      stabilize=stabilize,
                      algorithm=algorithm,
                      statistic=statistic,
                      statargs=statargs,
                      rmax = rmax)
  Xenv <- list2env(as.list(covariates), parent=parent.frame())
  X <- eval(substitute(X), envir=Xenv, enclos=parent.frame())
  if(is.NAobject(X)) {
    ## response pattern is missing; result is an NA object
    DPP <- list(...)$DPP
    cla <- if(!is.null(DPP)) "dppm" else "kppm"
    result <- NAobject(cla)
    return(result)
  }
  isquad <- is.quad(X)
  if(!is.ppp(X) && !isquad)
    stop("X should be a point pattern (ppp) or quadrature scheme (quad)")
  if(is.marked(X))
    stop("Sorry, cannot handle marked point patterns")
  if(!missing(subset)) {
    W <- eval(subset, covariates, parent.frame())
    if(!is.null(W)) {
      if(is.im(W)) {
        W <- solutionset(W)
      } else if(!is.owin(W)) {
        stop("Argument 'subset' should yield a window or logical image",
             call.=FALSE)
      }
      X <- X[W]
    }
  }

  po <- ppm(Q=X, trend=trend, covariates=covariates,
            forcefit=TRUE, rename.intercept=FALSE,
            covfunargs=covfunargs, use.gam=use.gam, nd=nd, eps=eps,
            improve.type=ppm.improve.type, improve.args=ppm.improve.args)
  XX <- if(isquad) X$data else X
  
  if(is.character(weightfun)) {
    RmaxW <- (rmax %orifnull% rmax.rule("K", Window(XX), intensity(XX))) / 2
    switch(weightfun,
           threshold = {
             weightfun <- function(d) { as.integer(d <= RmaxW) }
             attr(weightfun, "selfprint") <-
               paste0("Indicator(distance <= ", RmaxW, ")")
           },
           taper = {
             weightfun <- function(d) { pmin(1, RmaxW/d)^2 }
             attr(weightfun, "selfprint") <-
               paste0("min(1, ", RmaxW, "/d)^2")
           },
           stop(paste("Unrecognised option", sQuote(weightfun), "for weightfun"))
           )
  }
           
  ## default weight function
  if(is.null(weightfun))
    switch(method,
           adapcl = {
             weightfun <- function(d) { as.integer(abs(d) <= 1)*exp(1/(d^2-1)) }
             attr(weightfun, "selfprint") <-
               "Indicator(-1 <= distance <= 1) * exp(1/(distance^2-1))"
           },
           mincon = { },
           {
             RmaxW <- (rmax %orifnull% rmax.rule("K", Window(XX),
                                                 intensity(XX))) / 2
             weightfun <- function(d) { as.integer(d <= RmaxW) }
             attr(weightfun, "selfprint") <-
               paste0("Indicator(distance <= ", RmaxW, ")")
           })
  
  ## fit
  out <- switch(method,
         mincon = kppmMinCon(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, stabilize=stabilize,
                             statistic=statistic,
                             statargs=statargs, rmax=rmax,
                             algorithm=algorithm,
                             penalised = penalised,
                             trajectory = trajectory,
                             ...),
         clik2   = kppmComLik(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, stabilize=stabilize,
                             weightfun=weightfun, 
                             rmax=rmax, algorithm=algorithm,
                             penalised = penalised,
                             trajectory = trajectory,
                             ...),
         palm   = kppmPalmLik(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, stabilize=stabilize,
                             weightfun=weightfun, 
                             rmax=rmax, algorithm=algorithm,
                             penalised = penalised,
                             trajectory = trajectory,
                             ...),
         waag   = kppmWaagLik(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, stabilize=stabilize,
                             weightfun=weightfun, 
                             rmax=rmax, algorithm=algorithm,
                             penalised = penalised,
                             trajectory = trajectory,
                             ...),
         adapcl   = kppmCLadap(X=XX, Xname=Xname, po=po, clusters=clusters,
                             control=control, epsilon=epsilon, 
                             weightfun=weightfun, rmax=rmax,
                             algorithm=algorithm,
                             ...))
  ##
  h <- attr(out, "h")
  out <- append(out, list(ClusterArgs=ClusterArgs,
                          call=cl,
                          callframe=parent.frame(),
                          callstring=callstring))
  # Detect DPPs
  DPP <- list(...)$DPP
  class(out) <- unique(c(ifelse(is.null(DPP), "kppm", "dppm"), class(out)))

  # Update intensity estimate with improve.kppm if necessary:
  if(improve.type != "none")
    out <- do.call(improve.kppm,
                   append(list(object = out, type = improve.type),
                          improve.args))
  
  if(!is.null(h)) class(h) <- unique(c("traj", class(h)))

  attr(out, "h") <- h
  return(out)
}


improve.kppm <- local({

  fnc <- function(r, eps, g){ (g(r) - 1)/(g(0) - 1) - eps}

  improve.kppm <- function(object, type=c("quasi", "wclik1", "clik1"),
                           rmax = NULL, eps.rmax = 0.01,
                           dimyx = 50, maxIter = 100, tolerance = 1e-06,
                           fast = TRUE, vcov = FALSE, fast.vcov = FALSE,
                           verbose = FALSE,
                           save.internals = FALSE) {
    verifyclass(object, "kppm")
    type <- match.arg(type)
    gfun <- pcfmodel(object)
    X <- object$X
    win <- as.owin(X)
    ## simple (rectangular) grid quadrature scheme
    ## (using pixels with centers inside owin only)
    mask <- as.mask(win, dimyx = dimyx)
    wt <- pixellate(win, W = mask)
    wt <- wt[mask]
    Uxy <- rasterxy.mask(mask)
    U <- ppp(Uxy$x, Uxy$y, window = win, check=FALSE)
    U <- U[mask]
#    nU <- npoints(U)
    Yu <- pixellate(X, W = mask)
    Yu <- Yu[mask]
    
    ## covariates at quadrature points
    po <- object$po
    Z <- model.images(po, mask)
    Zu <- sapply(Z, "[", i=U, drop=FALSE)

    ## remove pixels which have NA values in covariates
    ok <- complete.cases(Zu)
    if(!all(ok)) {
      Zu <- if(is.null(dim(Zu))) Zu[ok] else Zu[ok, , drop=FALSE] 
      Yu <- Yu[ok]
      wt <- wt[ok]
      M <- mask
      M$m[M$m] <- ok
      U <- U[M]
    }

    ##obtain initial beta estimate using composite likelihood
    beta0 <- coef(po)
    
    ## determining the dependence range
    if (type != "clik1" && is.null(rmax))
      {
        diamwin <- diameter(win)
        rmax <- if(fnc(diamwin, eps.rmax, gfun) >= 0) diamwin else
                uniroot(fnc, lower = 0, upper = diameter(win),
                        eps=eps.rmax, g=gfun)$root
        if(verbose)
          splat(paste0("type: ", type, ", ",
                       "dependence range: ", rmax, ", ",
                       "dimyx: ", dimyx, ", g(0) - 1:", gfun(0) -1))
      }
    ## preparing the WCL case
    if (type == "wclik1")
      Kmax <- 2*pi * integrate(function(r){r * (gfun(r) - 1)},
                               lower=0, upper=rmax)$value * exp(c(Zu %*% beta0))
    ## the g()-1 matrix without tapering
    if (!fast || (vcov && !fast.vcov)){
      if (verbose)
        cat("computing the g(u_i,u_j)-1 matrix ...")
      gminus1 <- matrix(gfun(c(pairdist(U))) - 1, U$n, U$n)
      if (verbose)
        cat("..Done.\n")
    }
    
    if ( (fast && type == "quasi") | fast.vcov ){
      if (verbose)
        cat("computing the sparse G-1 matrix ...\n")
      ## Non-zero gminus1 entries (when using tapering)
      cp <- crosspairs(U,U,rmax,what="ijd")
      if (verbose)
        cat("crosspairs done\n")
      Gtap <- (gfun(cp$d) - 1)
      if(vcov){
        if(fast.vcov){
          gminus1 <- Matrix::sparseMatrix(i=cp$i, j=cp$j,
                                          x=Gtap, dims=c(U$n, U$n))
        } else{
          if(fast)
            gminus1 <- matrix(gfun(c(pairdist(U))) - 1, U$n, U$n)
        }
      }
      if (verbose & type!="quasi")
        cat("..Done.\n")
    }
       
    if (type == "quasi" && fast){
      mu0 <- exp(c(Zu %*% beta0)) * wt
      mu0root <- sqrt(mu0)
      sparseG <- Matrix::sparseMatrix(i=cp$i, j=cp$j,
                                      x=mu0root[cp$i] * mu0root[cp$j] * Gtap,
                                      dims=c(U$n, U$n))
      Rroot <- Matrix::Cholesky(sparseG, perm = TRUE, Imult = 1)
      ##Imult=1 means that we add 1*I
      if (verbose)
        cat("..Done.\n")
    }
    
    ## iterative weighted least squares/Fisher scoring
    bt <- beta0
    noItr <- 1
    repeat {
      mu <- exp(c(Zu %*% bt)) * wt
      mu.root <- sqrt(mu)
      ## the core of estimating equation: ff=phi
      ## in case of ql, \phi=V^{-1}D=V_\mu^{-1/2}x where (G+I)x=V_\mu^{1/2} Z
      ff <- switch(type,
                   clik1 = Zu,
                   wclik1= Zu/(1 + Kmax),
                   quasi = if(fast){
                     Matrix::solve(Rroot, mu.root * Zu)/mu.root
                   } else{
                     solve(diag(U$n) + t(gminus1 * mu), Zu)
                   }
                   )
      ##alternative
      ##R=chol(sparseG+sparseMatrix(i=c(1:U$n),j=c(1:U$n),
      ##                            x=rep(1,U$n),dims=c(U$n,U$n)))
      ##ff2 <- switch(type,
      ##              clik1 = Zu,
      ##              wclik1= Zu/(1 + Kmax),
      ##              quasi = if (fast)
      ##                         solve(R,solve(t(R), mu.root * Zu))/mu.root
      ##                      else solve(diag(U$n) + t(gminus1 * mu), Zu))
      ## print(summary(as.numeric(ff)-as.numeric(ff2)))
      ## the estimating equation: u_f(\beta)
      uf <- (Yu - mu) %*% ff
      ## inverse of minus expectation of Jacobian matrix: I_f
      Jinv <- solve(t(Zu * mu) %*% ff)
      if(maxIter==0){
        ## This is a built-in early exit for vcov internal calculations
        break
      }
      deltabt <- as.numeric(uf %*% Jinv)
      if (any(!is.finite(deltabt))) {
        warning(paste("Infinite value, NA or NaN appeared",
                      "in the iterative weighted least squares algorithm.",
                      "Returning the initial intensity estimate unchanged."),
                call.=FALSE)
        return(object)
      }
      ## updating the present estimate of \beta
      bt <- bt + deltabt
      if (verbose)
        splat(paste0("itr: ", noItr, ",\nu_f: ", as.numeric(uf),
                     "\nbeta:", bt, "\ndeltabeta:", deltabt))
      if (max(abs(deltabt/bt)) <= tolerance || max(abs(uf)) <= tolerance)
        break
      if (noItr > maxIter)
        stop("Maximum number of iterations reached without convergence.")
      noItr <- noItr + 1
    }
    out <- object
    out$po$coef.orig <- beta0
    out$po$coef <- bt
    loc <- if(is.sob(out$lambda)) as.mask(out$lambda) else mask
    out$lambda <- predict(out$po, locations = loc)
    out$improve <- list(type = type,
                        rmax = rmax,
                        dimyx = dimyx,
                        fast = fast,
                        fast.vcov = fast.vcov)
    if(save.internals){
      out$improve <- append(out$improve, list(ff=ff, uf=uf, J.inv=Jinv))
    }
    if(vcov){
      if (verbose)
        cat("computing the asymptotic variance ...\n")
      ## variance of the estimation equation: Sigma_f = Var(u_f(bt))
      trans <- if(fast) Matrix::t else t
      Sig <- trans(ff) %*% (ff * mu) + trans(ff * mu) %*% gminus1 %*% (ff * mu)
      ## note Abdollah's G does not have mu.root inside...
      ## the asymptotic variance of \beta:
      ##         inverse of the Godambe information matrix
      out$vcov <- as.matrix(Jinv %*% Sig %*% Jinv)
    }
    return(out)
  }
  improve.kppm
})


is.kppm <- function(x) { inherits(x, "kppm")}

print.kppm <- print.dppm <- function(x, ...) {

  isPCP <- x$isPCP
  # detect DPP
  isDPP <- inherits(x, "dppm")
  # handle outdated objects - which were all cluster processes
  if(!isDPP && is.null(isPCP)) isPCP <- TRUE

  terselevel <- spatstat.options('terse')
  digits <- getOption('digits')
  
  splat(if(x$stationary) "Stationary" else "Inhomogeneous",
        if(isDPP) "determinantal" else if(isPCP) "cluster" else "Cox",
        "point process model")

  Xname <- x$Xname
  if(waxlyrical('extras', terselevel) && nchar(Xname) < 20) {
    has.subset <- ("subset" %in% names(x$call))
    splat("Fitted to",
          if(has.subset) "(a subset of)" else NULL,
          "point pattern dataset", sQuote(Xname))
  }

  if(waxlyrical('gory', terselevel)) {
    fittedby <- "Fitted by"
    #' detect whether fit used a penalty
    if(!is.null(x$Fit$pspace$penalty))
      fittedby <- "Fitted by penalised"
    switch(x$Fit$method,
           mincon = {
             splat(fittedby, "minimum contrast")
             splat("\tSummary statistic:", x$Fit$StatName)
           },
           clik  =,
           clik2 = {
             splat(fittedby, "maximum Guan second order composite likelihood")
             splat("\trmax =", x$Fit$rmax)
             if(!is.null(wtf <- x$Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
           },
           palm = {
             splat(fittedby, "maximum Palm likelihood")
             splat("\trmax =", x$Fit$rmax)
             if(!is.null(wtf <- x$Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
           },
           waag = {
             splat(fittedby, "maximum Waagepetersen second order composite likelihood")
             splat("\trmax =", x$Fit$rmax)
             if(!is.null(wtf <- x$Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
           },
           adapcl = {
             splat("Fitted by adaptive second order composite likelihood")
             splat("\tepsilon =", x$Fit$epsilon)
             if(!is.null(wtf <- x$Fit$weightfun)) {
               a <- attr(wtf, "selfprint") %orifnull% pasteFormula(wtf)
               splat("\tweight function:", a)
             }
           },
           warning(paste("Unrecognised fitting method", sQuote(x$Fit$method)))
           )
  }

  #' optimization trace
  if(!is.null(attr(x, "h")))
    splat("[Includes history of evaluations of objective function]")
  
  parbreak(terselevel)
  
  # ............... trend .........................

  if(!(isDPP && is.null(x$fitted$intensity)))
    print(x$po, what="trend")

  # ..................... clusters ................

  # DPP case
  if(isDPP){
      splat("Fitted DPP model:")
      print(x$fitted)
      return(invisible(NULL))
  }

  tableentry <- spatstatClusterModelInfo(x$clusters)
  
  splat(if(isPCP) "Cluster" else "Cox",
        "model:", tableentry$printmodelname(x))
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
  pc <- x$par.canon
  if(!is.null(pc)) {
    splat("Fitted canonical parameters:")
    print(pc, digits=digits)
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
  if(isDPP) {
    rx <- repul(x)
    splat(if(is.im(rx)) "(Average) strength" else "Strength",
          "of repulsion:", signif(mean(rx), 4))
  }
  if(isPCP) {
    parbreak(terselevel)
    g <- pcfmodel(x)
    phi <- g(0) - 1
    splat("Cluster strength: phi = ", signif(phi, 4))
    psib <- phi/(1+phi)
    splat("Sibling probability: psib = ", signif(psib, 4))
    
  }
  invisible(NULL)
}

plot.kppm <- local({

  plotem <- function(x, ..., main=dmain, dmain) { plot(x, ..., main=main) }
  
  plot.kppm <- function(x, ...,
                        what=c("intensity", "statistic", "cluster"),
                        pause=interactive(),
                        xname) {
    ## catch objectname from dots if present otherwise deparse x:
    if(missing(xname)) xname <- short.deparse(substitute(x))
    nochoice <- missing(what)
    what <- pickoption("plot type", what,
                       c(statistic="statistic",
                         intensity="intensity",
                         cluster="cluster"),
                       multi=TRUE)
    ## handle older objects
    Fit <- x$Fit
    if(is.null(Fit)) {
      warning("kppm object is in outdated format")
      Fit <- x
      Fit$method <- "mincon"
    }
    ## Catch locations for clusters if given
    loc <- list(...)$locations
    inappropriate <- (nochoice & ((what == "intensity") & (x$stationary))) |
             ((what == "statistic") & (Fit$method != "mincon")) |
             ((what == "cluster") & (identical(x$isPCP, FALSE))) | 
             ((what == "cluster") & (!x$stationary) & is.null(loc))

    if(!nochoice && !x$stationary && "cluster" %in% what && is.null(loc))
      stop("Please specify additional argument ", sQuote("locations"),
           " which will be passed to the function ",
           sQuote("clusterfield"), ".")

    if(any(inappropriate)) {
      what <- what[!inappropriate]
      if(length(what) == 0){
        message("Nothing meaningful to plot. Exiting...")
        return(invisible(NULL))
      }
    }
    pause <- pause && (length(what) > 1)
    if(pause) {
       opa <- par(ask=TRUE)
       on.exit(par(opa))
    }
    for(style in what)
      switch(style,
             intensity={
               plotem(x$po, ...,
                      dmain=c(xname, "Intensity"),
                      how="image", se=FALSE)
             },
             statistic={
               plotem(Fit$mcfit, ...,
                      dmain=c(xname, Fit$StatName))
             },
             cluster={
               plotem(clusterfield(x, locations = loc, verbose=FALSE), ...,
                      dmain=c(xname, "Fitted cluster"))
             })
    return(invisible(NULL))
  }

  plot.kppm
})


predict.kppm <- predict.dppm <- function(object, ...) {
  se <- resolve.1.default(list(se=FALSE), list(...))
  interval <- resolve.1.default(list(interval="none"), list(...))
  if(se) warning("Standard error calculation assumes a Poisson process")
  if(interval != "none")
    warning(paste(interval, "interval calculation assumes a Poisson process"))
  predict(as.ppm(object), ...)
}

fitted.kppm <- fitted.dppm <- function(object, ...) {
  fitted(as.ppm(object), ...)
}

residuals.kppm <- residuals.dppm <- function(object, ...) {
  type <- resolve.1.default(list(type="raw"), list(...))
  if(type != "raw")
    warning(paste("calculation of", type, "residuals",
                  "assumes a Poisson process"))
  residuals(as.ppm(object), ...)
}


formula.kppm <- formula.dppm <- function(x, ...) {
  formula(x$po, ...)
}

terms.kppm <- terms.dppm <- function(x, ...) {
  terms(x$po, ...)
}

labels.kppm <- labels.dppm <- function(object, ...) {
  labels(object$po, ...)
}

update.kppm <- update.dppm <- function(object, ..., evaluate=TRUE, envir=environment(terms(object))) {
  argh <- list(...)
  nama <- names(argh)
  callframe <- object$callframe
  #' look for a formula argument
  fmla <- formula(object)
  jf <- integer(0)
  if(!is.null(trend <- argh$trend)) {
    if(!can.be.formula(trend))
      stop("Argument \"trend\" should be a formula")
    fmla <- newformula(formula(object), trend, callframe, envir)
    jf <- which(nama == "trend")
  } else if(any(isfo <- sapply(argh, can.be.formula))) {
    if(sum(isfo) > 1) {
      if(!is.null(nama)) isfo <- isfo & nzchar(nama)
      if(sum(isfo) > 1)
        stop(paste("Arguments not understood:",
                   "there are two unnamed formula arguments"))
    }
    jf <- which(isfo)
    fmla <- argh[[jf]]
    fmla <- newformula(formula(object), fmla, callframe, envir)
  }

  #' look for a point pattern or quadscheme
  if(!is.null(X <- argh$X)) {
    if(!inherits(X, c("ppp", "quad")))
      stop(paste("Argument X should be a formula,",
                 "a point pattern or a quadrature scheme"))
    jX <- which(nama == "X")
  } else if(any(ispp <- sapply(argh, inherits, what=c("ppp", "quad")))) {
    if(sum(ispp) > 1) {
      if(!is.null(nama)) ispp <- ispp & nzchar(nama)
      if(sum(ispp) > 1)
        stop(paste("Arguments not understood:",
                   "there are two unnamed point pattern/quadscheme arguments"))
    }
    jX <- which(ispp)
    X <- argh[[jX]]
  } else {
    X <- object$X
    jX <- integer(0)
  }
  Xexpr <- if(length(jX) > 0) sys.call()[[2L + jX]] else NULL
  #' remove arguments just recognised, if any
  jused <- c(jf, jX)
  if(length(jused) > 0) {
    argh <- argh[-jused]
    nama <- names(argh)
  }
  #' update the matched call
  thecall <- getCall(object)
  methodname <- as.character(thecall[[1L]])
  switch(methodname,
         kppm.formula = {
	   ## original call has X = [formula with lhs]
	   if(!is.null(Xexpr)) {
	     lhs.of.formula(fmla) <- Xexpr
	   } else if(is.null(lhs.of.formula(fmla))) {
	     lhs.of.formula(fmla) <- as.name('.')
	   }
           oldformula <- getCall(object)$X
           oldformula <- eval(oldformula, callframe)
           thecall$X <- newformula(oldformula, fmla, callframe, envir)
         },
         {
	   ## original call has X = ppp and trend = [formula without lhs]
           oldformula <- getCall(object)$trend %orifnull% (~1)
           oldformula <- eval(oldformula, callframe)
	   fom <-  newformula(oldformula, fmla, callframe, envir)
	   if(!is.null(Xexpr))
	      lhs.of.formula(fom) <- Xexpr
	   if(is.null(lhs.of.formula(fom))) {
	      # new call has same format
	      thecall$trend <- fom
  	      if(length(jX) > 0)
  	        thecall$X <- X
	   } else {
	      # new call has formula with lhs
	      thecall$trend <- NULL
	      thecall$X <- fom
	   }
         })
  knownnames <- unique(c(names(formals(kppm.ppp)),
                         names(formals(mincontrast)),
                         names(formals(optim))))
  knownnames <- setdiff(knownnames,
                        c("X", "trend",
                          "observed", "theoretical",
                          "fn", "gr", "..."))
  ok <- nama %in% knownnames
  thecall <- replace(thecall, nama[ok], argh[ok])
  thecall$formula <- NULL # artefact of 'step', etc
  thecall[[1L]] <- as.name("kppm")
  if(!evaluate)
    return(thecall)
  out <- eval(thecall, envir=parent.frame(), enclos=envir)
  #' update name of data
  if(length(jX) == 1) {
    mc <- match.call()
    Xlang <- mc[[2L+jX]]
    out$Xname <- short.deparse(Xlang)
  }
  #'
  return(out)
}

updateData.kppm <- function(model, X, ...) {
  update(model, X=X)
}

unitname.kppm <- unitname.dppm <- function(x) {
  return(unitname(x$X))
}

"unitname<-.kppm" <- "unitname<-.dppm" <- function(x, value) {
  unitname(x$X) <- value
  if(!is.null(x$Fit$mcfit)) {
    unitname(x$Fit$mcfit) <- value
  } else if(is.null(x$Fit)) {
    warning("kppm object in outdated format")
    if(!is.null(x$mcfit))
      unitname(x$mcfit) <- value
  }
  return(x)
}

as.fv.kppm <- as.fv.dppm <- function(x) {
  if(x$Fit$method == "mincon") 
    return(as.fv(x$Fit$mcfit))
  gobs <- if(is.stationary(x)) {
            pcf(x$X,
                correction="good", divisor="a", zerocor="J")
          } else {
            pcfinhom(x$X,
                     lambda=x, update=FALSE,
                     correction="good", divisor="a", zerocor="J")
          }
  gname <- attr(gobs, "fname")
  gfit <- (pcfmodel(x))(gobs$r)
  g <- bind.fv(gobs,
               data.frame(fit=gfit),
               labl = makefvlabel(fname=gname, sub="fit"),
               desc = "predicted %s for fitted model")
  return(g)
}

coef.kppm <- coef.dppm <- function(object, ...) {
  return(coef(object$po))
}


Kmodel.kppm <- function(model, ...) {
  Kpcf.kppm(model, what="K")
}

pcfmodel.kppm <- function(model, ...) {
  Kpcf.kppm(model, what="pcf")
}

Kpcf.kppm <- function(model, what=c("K", "pcf", "kernel")) {
  what <- match.arg(what)
  # Extract function definition from internal table
  clusters <- model$clusters
  tableentry <- spatstatClusterModelInfo(clusters)
  if(is.null(tableentry))
    stop("No information available for", sQuote(clusters), "cluster model")
  fun <- tableentry[[what]]
  if(is.null(fun))
    stop("No expression available for", what, "for", sQuote(clusters),
         "cluster model")
  # Extract model parameters
  par <- model$par
  # Extract covariance model (if applicable)
  cm <- model$covmodel
  model <- cm$model
  margs <- cm$margs
  #
  f <- function(r) as.numeric(fun(par=par, rvals=r,
                                  model=model, margs=margs))
  return(f)
}

is.stationary.kppm <- is.stationary.dppm <- function(x) {
  return(x$stationary)
}

is.poisson.kppm <- function(x) {
  switch(x$clusters,
         Cauchy=,
         VarGamma=,
         Thomas=,
         MatClust={
           # Poisson cluster process
           mu <- x$mu
           return(!is.null(mu) && (max(mu) == 0))
         },
         LGCP = {
           # log-Gaussian Cox process
           sigma2 <- x$par[["sigma2"]]
           return(sigma2 == 0)
         },
         return(FALSE))
}

# extract ppm component

as.ppm.kppm <- as.ppm.dppm <- function(object) {
  object$po
}

# other methods that pass through to 'ppm'

as.owin.kppm <- as.owin.dppm <- function(W, ..., from=c("points", "covariates"), fatal=TRUE) {
  from <- match.arg(from)
  as.owin(as.ppm(W), ..., from=from, fatal=fatal)
}

domain.kppm <- Window.kppm <- domain.dppm <-
  Window.dppm <- function(X, ..., from=c("points", "covariates")) {
  from <- match.arg(from)
  as.owin(X, from=from)
}

model.images.kppm <-
  model.images.dppm <- function(object, W=as.owin(object), ...) {
  model.images(as.ppm(object), W=W, ...)
}

model.matrix.kppm <-
  model.matrix.dppm <- function(object,
                                data=model.frame(object, na.action=NULL), ...,
                                Q=NULL, 
                                keepNA=TRUE) {
  if(missing(data)) data <- NULL
  model.matrix(as.ppm(object), data=data, ..., Q=Q, keepNA=keepNA)
}

model.frame.kppm <- model.frame.dppm <- function(formula, ...) {
  model.frame(as.ppm(formula), ...)
}

logLik.kppm <- logLik.dppm <- function(object, ...) {
  cl <- object$Fit$maxlogcl
  if(is.null(cl))
    stop(paste("logLik is only available for kppm objects fitted with",
               "method='palm' or method='clik2'"),
         call.=FALSE)
  ll <- logLik(as.ppm(object)) # to inherit class and d.f.
  ll[] <- cl
  return(ll)
}
 
AIC.kppm <- AIC.dppm <- function(object, ..., k=2) {
  cl <- logLik(object)
  df <- attr(cl, "df")
  return(- 2 * as.numeric(cl) + k * df)
}

extractAIC.kppm <- extractAIC.dppm <- function (fit, scale = 0, k = 2, ...) {
  cl <- logLik(fit)
  edf <- attr(cl, "df")
  aic <- - 2 * as.numeric(cl) + k * edf
  return(c(edf, aic))
}

nobs.kppm <- nobs.dppm <- function(object, ...) { nobs(as.ppm(object)) }

psib <- function(object) UseMethod("psib")

psib.kppm <- function(object) {
  clus <- object$clusters
  info <- spatstatClusterModelInfo(clus)
  if(!info$isPCP) {
    warning("The model is not a cluster process")
    return(NA)
  }
  g <- pcfmodel(object)
  p <- 1 - 1/g(0)
  return(p)
}

reach.kppm <- function(x, ..., epsilon) {
  thresh <- if(missing(epsilon)) NULL else epsilon
  if(x$isPCP) 
    return(2 * clusterradius.kppm(x, ..., thresh=thresh))
  ## use pair correlation
  g <- pcfmodel(x)
  ## find upper bound
  if(is.null(thresh)) thresh <- 0.01
  f <- function(r) { g(r) - 1 - thresh }
  scal <- as.list(x$par)$scale %orifnull% 1
  for(a in scal * 2^(0:10)) { if(f(a) < 0) break; }
  if(f(a) > 0) return(Inf)
  ## solve g(r) = 1 + epsilon
  b <- uniroot(f, c(0, a))$root
  return(b)
}

is.poissonclusterprocess <- function(model) {
  UseMethod("is.poissonclusterprocess")
}

is.poissonclusterprocess.default <- function(model) { FALSE }
  
is.poissonclusterprocess.kppm <- function(model) { isTRUE(model$isPCP) }

#' cluster strength

clusterstrength <- function(object) {
  verifyclass(object, "kppm")
  g <- pcfmodel(object)
  phi <- g(0) - 1
  return(phi)
}

#' spatial persistence index

persist <- function(object, W=Window(object)) {
  verifyclass(object, "kppm")
  stopifnot(is.owin(W))
  g <- pcfmodel(object)
  d <- diameter(W)
  v <- (g(d)-1)/(g(0)-1)
  return(v)
}
