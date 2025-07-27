#'
#'   intensity.ppm.R
#'
#'    Intensity and intensity approximations for fitted point process models
#'
#'    $Revision: 1.7 $ $Date: 2025/07/27 06:50:57 $
#'
#'    Adrian Baddeley with contributions from Frederic Lavancier

intensity.ppm <- function(X, ..., approx=c("Poisson", "DPP")) {
  approx <- match.arg(approx)
  if(!isTRUE(valid.ppm(X))) {
    warning("Model is invalid - projecting it")
    X <- project.ppm(X)
  }
  if(is.poisson(X)) {
    if(is.stationary(X)) {
      # stationary univariate/multivariate Poisson
      sX <- summary(X, quick="no variances")
      lam <- sX$trend$value
      if(sX$multitype && sX$no.trend) {
        ## trend is ~1; lam should be replicated for each mark
        lev <- levels(marks(data.ppm(X)))
        lam <- rep(lam, length(lev))
        names(lam) <- lev
      }
      return(lam)
    }
    # Nonstationary Poisson
    return(predict(X, ...))
  }
  # Gibbs process
  if(is.multitype(X))
    stop("Not yet implemented for multitype Gibbs processes")
  # Compute first order term
  if(is.stationary(X)) {
    ## activity parameter
    sX <- summary(X, quick="no variances")
    beta <- sX$trend$value
  } else {
    ## activity function (or values of it, depending on '...')
    beta <- predict(X, ...)
  }
  ## apply approximation
  lambda <- SaddleApprox(beta, fitin(X), approx=approx)
  return(lambda)
}

SaddleApprox <- function(beta, fi, approx=c("Poisson", "DPP")) {
  approx <- match.arg(approx)
  z <- switch(approx,
              Poisson = PoisSaddle(beta, fi),
              DPP     = DPPSaddle(beta, fi))
  return(z)
}

PoisSaddle <- function(beta, fi) {
  ## apply Poisson-Saddlepoint approximation
  ## given first order term and fitted interaction
  stopifnot(inherits(fi, "fii"))
  if(interactionorder(fi) == 2)
    return(PoisSaddlePairwise(beta, fi))
  modelname <- as.interact(fi)$name
  if(identical(modelname, "Geyer saturation process"))
    return(PoisSaddleGeyer(beta, fi))
  if(identical(modelname, "Area-interaction process"))
    return(PoisSaddleArea(beta, fi))
  stop(paste("Poisson-saddlepoint intensity approximation",
             "is not yet available for",
             modelname), call.=FALSE)
}

PoisSaddlePairwise <- function(beta, fi) {
  # compute second Mayer cluster integral
  G <- Mayer(fi)
  if(is.null(G) || !is.finite(G)) 
    stop("Internal error in computing Mayer cluster integral")
  if(G < 0)
    stop(paste("Unable to apply Poisson-saddlepoint approximation:",
               "Mayer cluster integral is negative"))
  ## solve
  if(is.im(beta)) {
    lambda <- if(G == 0) beta else eval.im(LambertW(G * beta)/G)
  } else {
    lambda <- if(G == 0) beta else (LambertW(G * beta)/G)
    if(length(lambda) == 1) lambda <- unname(lambda)
  }
  return(lambda)
}

# Lambert's W-function

LambertW <- local({

  yexpyminusx <- function(y,x){y*exp(y)-x}

  W <- function(x) {
    result <- rep.int(NA_real_, length(x))
    ok <- is.finite(x) & (x >= 0)
    if(requireNamespace("gsl", quietly=TRUE)) {
      result[ok] <- gsl::lambert_W0(x[ok])
    } else {
      for(i in which(ok))
        result[i] <- uniroot(yexpyminusx, c(0, x[i]), x=x[i])$root
    }
    return(result)
  }

  W
})

PoisSaddleGeyer <- local({

  PoisSaddleGeyer <- function(beta, fi) {
    gamma <- summary(fi)$sensible$param$gamma
    if(gamma == 1) return(beta)
    inte <- as.interact(fi)
    sat <- inte$par$sat
    R   <- inte$par$r
    if(sat == 0 || R == 0) return(beta)
    #' get probability distribution of Geyer statistic under reference model
    z <- Spatstat.Geyer.Nulldist # from sysdata
    if(is.na(m <- match(sat, z$sat)))
      stop(paste("Sorry, the Poisson-saddlepoint approximation",
                 "is not implemented for Geyer models with sat =",
                 paste0(sat, ";"),
                 "supported values are sat =", commasep(z$sat)),
           call.=FALSE)
    ## extract frequency table of values of Geyer statistic T(0, Z_n)
    ## from huge simulation, where Z_n is runifdisc(n, radius = 2*R).
    freqtable <- z$freq[[m]]
    possT <- attr(freqtable, "possvals")
    possN <- 0:(nrow(freqtable)-1)
    nsim <- attr(freqtable, "nsim") %orifnull% sum(freqtable[1,])
    ## normalise to obtain probability distributions
    probmatrix <- freqtable/nsim
    #' Compute conditional expectation of gamma^T given N=n for each n
    gammaT <- gamma^possT
    EgamTgivenN <- probmatrix %*% gammaT
    #' precompute some constants
    gammasat <- gamma^sat
    gammaTrange <- range(gammaT)
    fpiR2 <- 4 * pi * R^2
    #' apply approximation
    betavalues <- beta[]
    nvalues <- length(betavalues)
    lambdavalues <- numeric(nvalues)
    for(i in seq_len(nvalues)) {
      beta.i <- betavalues[i]
      lambdavalues[i] <- uniroot(diffapproxGeyer,
                                 interval    = beta.i * gammaTrange,
                                 beta        = beta.i,
                                 gammasat    = gammasat,
                                 fpiR2       = fpiR2,
                                 EgamTgivenN = EgamTgivenN,
                                 possN       = possN
                                 )$root
    }
    #' return result in same format as 'beta'
    lambda <- beta
    lambda[] <- lambdavalues
    if(length(lambda) == 1) lambda <- unname(lambda)
    return(lambda)
  }

  diffapproxGeyer <- function(lambda, beta,
                              gammasat,
                              fpiR2, EgamTgivenN, possN) {
    #' Compute approximation to E_Pois(lambda) Lambda(0,X) for Geyer
    pN <- dpois(possN, lambda * fpiR2)
    EgamT <- sum(pN * EgamTgivenN) + gammasat * (1-sum(pN))
    #' Last term on previous line assumes that, for n > max(possN),
    #' distribution of T is concentrated on T=sat
    lambdaApprox <- beta * EgamT
    return(lambdaApprox - lambda)
  }

  PoisSaddleGeyer
})

PoisSaddleArea <- local({

  PoisSaddleArea <- function(beta, fi) {
    eta <- summary(fi)$sensible$param$eta
    if(eta == 1) return(beta)
    etarange <- range(c(eta^2, 1.1, 0.9))
    inte <- as.interact(fi)
    R   <- inte$par$r
    #' reference distribution of canonical sufficient statistic
    zeroprob <- Spatstat.Area.Zeroprob
    areaquant <- Spatstat.Area.Quantiles
    # expectation of eta^A_n for each n = 0, 1, ....
    EetaAn <- c(1/eta,
                zeroprob + (1-zeroprob) * colMeans((eta^(-areaquant))))
    #' compute approximation
    betavalues <- beta[]
    nvalues <- length(betavalues)
    lambdavalues <- numeric(nvalues)
    for(i in seq_len(nvalues)) {
      beta.i <- betavalues[i]
      ra <- beta.i * etarange
      lambdavalues[i] <- uniroot(diffapproxArea, ra, beta=beta.i,
                                 eta=eta, r=R,
                                 EetaAn=EetaAn)$root
    }
    #' return result in same format as 'beta'
    lambda <- beta
    lambda[] <- lambdavalues
    if(length(lambda) == 1) lambda <- unname(lambda)
    return(lambda)
  }

  diffapproxArea <- function(lambda, beta, eta, r, EetaAn) {
    lambda - approxEpoisArea(lambda, beta, eta, r, EetaAn)
  }

  approxEpoisArea <- function(lambda, beta=1, eta=1, r=1, EetaAn) {
    #' Compute approximation to E_Pois(lambda) Lambda(0,X) for AreaInter
    mu <- lambda * pi * (2*r)^2
    zeta <- pi^2/2 - 1
    theta <-  -log(eta)
    zetatheta <- zeta * theta

    #' contribution from tabulated values
    Nmax <- length(EetaAn) - 1L
    possN <- 0:Nmax
    qN <- dpois(possN, mu)
    # expectation of eta^A when N ~ poisson (truncated)
    EetaA <- sum(qN * EetaAn)

    #' asymptotics for quite large n
    Nbig <- qpois(0.999, mu)
    qn <- 0
    if(Nbig > Nmax) {
      n <- (Nmax+1):Nbig
      #' asymptotic mean uncovered area conditional on this being positive
      mstarn <- (16/((n+3)^2)) * exp(n * (1/4 - log(4/3)))
      ztm <- zetatheta * mstarn
      ok <- (ztm < 1)
      if(!any(ok)) {
        Nbig <- Nmax
        qn <- 0
      } else {
        if(!all(ok)) {
          Nbig <- max(which(!ok)) - 1
          n <- (Nmax+1):Nbig
          ztm <- ztm[1:((Nbig-Nmax)+1)]
        }
        qn <- dpois(n, mu)
        #' asymptotic  probability of complete coverage
        pstarn <- 1 - pmin(1, 3 * (1 + n^2/(16*pi)) * exp(-n/4))
        Estarn <- (1 - ztm)^(-1/zeta)
        EetaA <- EetaA + sum(qn * (pstarn + (1-pstarn) * Estarn))
      }
    }
    #' for very large n, assume complete coverage, so A = 0
    EetaA <- EetaA + 1 - sum(qN) - sum(qn)
    return(beta * eta * EetaA)
  }

  PoisSaddleArea

})

## The following was contributed by Frederic Lavancier
## hacked by Adrian

DPPSaddle <- function(beta, fi) {
  if(interactionorder(fi) != 2)
    stop("DPP approximation is only available for pairwise interactions")
  if(!is.numeric(beta))
    stop("DPP approximation is only available for stationary models")
  if(length(beta) > 1)
    stop("DPP approximation is not available for multitype models")
  DPPSaddlePairwise(beta, fi)
}

DPPSaddlePairwise<-function (beta, fi) {
  stopifnot(inherits(fi, "fii"))
  ## second Mayer cluster integral
  G <- Mayer(fi)
  if (is.null(G) || !is.finite(G)) 
    stop("Internal error in computing Mayer cluster integral")
  if(G < (-exp(-1)/beta)){
    warning(paste("The second Mayer cluster integral",
                  "is less than exp(-1)/beta,",
                  "which may lead to an unreliable solution."))
  }
  ## integral of (1-g)^2
  G2 <- Mayer(fi, exponent=2)
  if (is.null(G2) || !is.finite(G2)) 
    stop("Internal error in computing integral of (1-g)^2")

  ## hard core distance
  hc <- hardcoredist(fi)

  ## interaction range
  R <- reach(fi, epsilon=0.001)
  if(!is.finite(R))
    stop("Unable to determine a finite range of interaction")

  ## solve
  Bhc <- pi * hc^2
  BR <- pi*R^2
  kappa<-(G2-Bhc)/(BR-Bhc)
  
  fun <-function(x) {
    log(beta) + (1+x*Bhc)*log(1-x*Bhc/(1+x*Bhc)) + 
      (1+x*(G-Bhc)/kappa)*log(1-x*(G-Bhc)/(1+x*(G-Bhc)/kappa)) - log(x)
  }
  lambda <- uniroot(fun,c(0,2*beta))$root

  return(lambda)
}

Mayer <- function(fi, exponent=1){
  stopifnot(inherits(fi, "fii"))
  ## compute second Mayer cluster integral for a PAIRWISE interaction
  if(exponent == 1) {
    ## check if there is an analytic expression
    inte <- as.interact(fi)
    MayerCode <- inte$Mayer
    if(is.function(MayerCode)) {
      ## interaction coefficients
      co <- with(fi, coefs[Vnames[!IsOffset]])
      z <- MayerCode(co, inte) # sic
      return(z)
    }
  }
  ## No specialised code provided.
  if(interactionorder(fi) != 2)
    stop("Mayer() is only defined for pairwise interactions")
  ## Compute by numerical integration
  f <- function(x) {
    log.g <- evalPairwiseTerm(fi, x)
    z <- 2 * pi * x * ifelse(is.finite(log.g),
                             (1 - exp(log.g))^exponent,
                             1)
    return(z)
  }
  R <- reach(fi)
  M <- integrate(f,0,R)$value
  return(M)
}

