#
#	Gcom.R
#
#	Model compensator of G 
#
#	$Revision: 1.12 $	$Date: 2023/02/02 02:39:48 $
#
################################################################################
#


Gcom <- function(object, r=NULL, breaks=NULL, ...,
                 correction=c("border", "Hanisch"),
                 conditional=!is.poisson(object),
                 restrict=FALSE,
                 model=NULL,
                 trend=~1, interaction=Poisson(),
                 rbord=reach(interaction),
                 ppmcorrection="border",
                 truecoef=NULL, hi.res=NULL) {
  if(is.ppm(object)) {
    fit <- object
  } else if(is.ppp(object) || is.quad(object)) {
    if(is.ppp(object)) object <- quadscheme(object, ...)
    if(!is.null(model)) {
      stopifnot(is.ppm(model))
      ## update 'model' using new point pattern 'object' 
      e <- list2env(list(object=object), parent=model$callframe)
      fit <- update(model, Q=object, forcefit=TRUE, envir=e)
    } else {
      fit <- ppm(object, trend=trend, interaction=interaction, rbord=rbord,
                 forcefit=TRUE)
    }
  } else 
    stop("object should be a fitted point process model or a point pattern")

  if(missing(conditional) || is.null(conditional))
    conditional <- !is.poisson(fit)
  
  restrict <- isTRUE(restrict)
  if(restrict && !conditional) {
    warning("restrict=TRUE ignored because conditional=FALSE", call.=FALSE)
    restrict <- FALSE
  }
  
#  rfixed <- !is.null(r) || !is.null(breaks)
  
  # selection of edge corrections
#  correction.given <- !missing(correction) && !is.null(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             Hanisch="Hanisch",
                             hanisch="Hanisch",
                             best="Hanisch"),
                           multi=TRUE)

  # Extract data and quadrature points
  Q <- quad.ppm(fit, drop=FALSE)
  X <- data.ppm(fit)
  Win <- X$window

  # edge correction algorithm 
  algo <- if(!conditional) "classical" else
          if(restrict) "restricted" else "reweighted"

  # conditioning on border region?
  if(!conditional) {
    Wfree <- Win
  } else {
    rbord <- fit$rbord
    Wfree <- erosion(Win, rbord)
    if(restrict) {
      retain <- inside.owin(union.quad(Q), , Wfree)
      Q <- Q[Wfree]
      X <- X[Wfree]
      Win <- Wfree
    } 
  }

  # Extract quadrature info
  U <- union.quad(Q)
  Z <- is.data(Q) # indicator data/dummy
#  E <- equalsfun.quad(Q)
  WQ <- w.quad(Q)  # quadrature weights

  # basic statistics
  npts <- npoints(X)
  areaW <- area(Win)
  lambda <- npts/areaW
  
  # quadrature points used
  USED <- if(algo == "reweighted") (bdist.points(U) > rbord) else rep.int(TRUE, U$n)
  
  # adjustments to account for restricted domain 
  if(conditional && spatstat.options("eroded.intensity")) {
    npts.used <- sum(Z & USED)
    area.used <- sum(WQ[USED])
    lambda.used <- npts.used/area.used
  } else {
    npts.used <- npts
    area.used <- areaW
    lambda.used <- lambda
  }
  
  #  determine breakpoints for r values
  rmaxdefault <- rmax.rule("G", if(restrict) Wfree else Win, lambda)
  breaks <- handle.r.b.args(r, breaks, Wfree, rmaxdefault=rmaxdefault)
  rvals <- breaks$r
  rmax  <- breaks$max
  
  # residuals
  resid <- residuals(fit, type="raw",drop=FALSE,
                    new.coef=truecoef, quad=hi.res)
  rescts  <- with(resid, "continuous")
  if(restrict) {
    # keep only data inside Wfree
    rescts  <- rescts[retain]
  }
  # absolute weight for continuous integrals
#  wc   <- -rescts

  # nearest neighbours (quadrature point to data point)
  nn <- nncross(U, X, seq(U$n), seq(X$n))
  dIJ <- nn$dist
  I <- seq(U$n)
#  J <- nn$which
  DD <- Z <- (I <= X$n)  # TRUE for data points
  wcIJ <- -rescts

  # determine whether a quadrature point will be used in integral
  okI <- USED[I]

   # initialise fv object
  r <- breaks$r
  df <- data.frame(r=r, pois=1 - exp(-pi * lambda.used * r^2))
  G <- fv(df, "r", substitute(G(r), NULL), "pois", . ~ r,
          alim=c(0, rmax),
          labl=c("r","%s[pois](r)"),
          desc=c("distance argument r", "theoretical Poisson %s"),
          fname="G")

  #  distance to boundary
  b <- bI <- bdist.points(U)

  dotnames <- character(0)

  # Border method
  if("border" %in% correction) {
    # reduced sample for G(r) of data only
    ZUSED <- Z & USED
    RSX <- Kount(dIJ[DD & okI], bI[DD & okI], b[ZUSED], breaks)
    Gb <- RSX$numerator/RSX$denom.count
    G <- bind.fv(G, data.frame(border=Gb), "hat(%s)[bord](r)",
                 "border-corrected nonparametric estimate of %s",
                 "border")
    # reduced sample for adjustment integral
    RSD <- Kwtsum(dIJ[okI], bI[okI], wcIJ[okI], b[ZUSED],
                  rep.int(1, sum(ZUSED)), breaks, fatal=FALSE)
    Gbcom <- RSD$numerator/(1 + RSD$denominator)
    
    G <- bind.fv(G, data.frame(bcom=Gbcom), "bold(C)~hat(%s)[bord](r)",
                 "model compensator of border-corrected %s",
                 "bcom")

    dotnames <- c("border", "bcom", "pois")
  }

  # Hanisch correction for data
  if("Hanisch" %in% correction) {
    nnd <- dIJ[DD & okI]
    bdry <- bI[DD & okI]
    # weights
    ea <- eroded.areas(Win, rvals)
    if(algo == "reweighted") {
      # replace weight(r) by weight(max(rbord,r))
      ea[rvals < rbord] <- eroded.areas(Win, rbord)
    }
    # compute
    x <- nnd[nnd <= bdry]
    h <- whist(x[x <= rmax], breaks=breaks$val)
    H <- (1/lambda.used) * cumsum(h/ea)
    # glue on 
    G <- bind.fv(G, data.frame(han=H), "hat(%s)[han](r)",
                 "Hanisch correction estimate of %s",
                 "han")
    # Hanisch correction for adjustment integral
    nnd <- dIJ[okI]
    bdry <- bI[okI]
    wt   <- wcIJ[okI]
    x <- nnd[nnd <= bdry]
    wt <- wt[nnd <= bdry]
    h <- whist(x[x <= rmax], breaks=breaks$val, weights=wt[x <= rmax])
    lambdaplus <- (npts.used + 1)/area.used
    Hint <- (1/lambdaplus) * cumsum(h/ea)
    # glue on 
    G <- bind.fv(G, data.frame(hcom=Hint), "bold(C)~hat(%s)[han](r)",
                 "model compensator of Hanisch-corrected %s",
                 "hcom")
    # pseudovariance for Hanisch residual
    Hvar <- (1/lambdaplus^2) * cumsum(h/ea^2)
    G <- bind.fv(G, data.frame(hvar=Hvar), "bold(C)^2~hat(%s)[han](r)",
                 "Poincare variance for Hanisch corrected %s",
                 "hcom")
    # default plot does not show all components
    dotnames <- c("han", "hcom", dotnames)
  }
  # compute sensible 'alim'
  endpoint <- function(y, r, f) { min(r[y >= f * max(y)]) }
  amax <- endpoint(G$pois, G$r, 0.99)
  if(length(dotnames) > 0) 
    amax <- max(amax,
                unlist(lapply(as.data.frame(G)[,dotnames,drop=FALSE],
                              endpoint,
                              r=r, f=0.9)))
  attr(G, "alim") <- c(0, amax)
  #  
  fvnames(G, ".") <- dotnames
  unitname(G) <- unitname(X)
  # secret tag used by 'Gres'
  attr(G, "maker") <- "Gcom"
  return(G)
}

