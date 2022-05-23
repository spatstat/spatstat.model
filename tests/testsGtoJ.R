#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.model
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.model)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#       
#        tests/hobjects.R
#
#   Validity of methods for ppm(... method="ho")
#


if(FULLTEST) {
local({
  set.seed(42)
  fit  <- ppm(cells ~1,         Strauss(0.1), method="ho", nsim=10)
  fitx <- ppm(cells ~offset(x), Strauss(0.1), method="ho", nsim=10)

  a  <- AIC(fit)
  ax <- AIC(fitx)

  f  <- fitted(fit)
  fx <- fitted(fitx)

  p  <- predict(fit)
  px <- predict(fitx)
})
}


#'     tests/hypotests.R
#'     Hypothesis tests
#' 
#'  $Revision: 1.9 $ $Date: 2020/11/02 06:39:23 $

if(FULLTEST) {
local({

  #' scan test with baseline
  fit <- ppm(cells ~ x)
  lam <- predict(fit)
  rr <- c(0.05, 1)
  scan.test(cells, rr, nsim=5,
            method="poisson", baseline=fit, alternative="less")
  scan.test(cells, rr, nsim=5,
            method="poisson", baseline=lam, alternative="less")
})
}
#'
#'  tests/interact.R
#'
#'  Support for interaction objects
#'
#'  $Revision: 1.2 $ $Date: 2020/04/28 12:58:26 $

if(FULLTEST) {
local({
  #' print.intermaker
  Strauss
  Geyer
  Ord
  #' intermaker
  BS <- get("BlankStrauss", envir=environment(Strauss))
  BD <- function(r) { instantiate.interact(BS, list(r=r)) }
  BlueDanube <- intermaker(BD, BS) 
})
}

#'   tests/ippm.R
#'   Tests of 'ippm' class
#'   $Revision: 1.6 $ $Date: 2020/04/28 12:58:26 $

if(FULLTEST) {
local({
  # .......... set up example from help file .................
  nd <- 10
  gamma0 <- 3
  delta0 <- 5
  POW <- 3
  # Terms in intensity
  Z <- function(x,y) { -2*y }
  f <- function(x,y,gamma,delta) { 1 + exp(gamma - delta * x^POW) }
  # True intensity
  lamb <- function(x,y,gamma,delta) { 200 * exp(Z(x,y)) * f(x,y,gamma,delta) }
  # Simulate realisation
  lmax <- max(lamb(0,0,gamma0,delta0), lamb(1,1,gamma0,delta0))
  set.seed(42)
  X <- rpoispp(lamb, lmax=lmax, win=owin(), gamma=gamma0, delta=delta0)
  # Partial derivatives of log f
  DlogfDgamma <- function(x,y, gamma, delta) {
    topbit <- exp(gamma - delta * x^POW)
    topbit/(1 + topbit)
  }
  DlogfDdelta <- function(x,y, gamma, delta) {
    topbit <- exp(gamma - delta * x^POW)
    - (x^POW) * topbit/(1 + topbit)
  }
  # irregular score
  Dlogf <- list(gamma=DlogfDgamma, delta=DlogfDdelta)
  # fit model
  fit <- ippm(X ~Z + offset(log(f)),
              covariates=list(Z=Z, f=f),
              iScore=Dlogf,
              start=list(gamma=1, delta=1),
              nd=nd)
  # fit model with logistic likelihood but without iScore
  fitlo <- ippm(X ~Z + offset(log(f)),
                method="logi",
                covariates=list(Z=Z, f=f),
                start=list(gamma=1, delta=1),
                nd=nd)

  ## ............. test ippm class support ......................
  Ar <- model.matrix(fit)
  Ai <- model.matrix(fit, irregular=TRUE)
  An <- model.matrix(fit, irregular=TRUE, keepNA=FALSE)
  AS <- model.matrix(fit, irregular=TRUE, subset=(abs(Z) < 0.5))

  Zr <- model.images(fit)
  Zi <- model.images(fit, irregular=TRUE)
  ## update.ippm
  fit2 <- update(fit, . ~ . + I(Z^2))
  fit0 <- update(fit,
                 . ~ . - Z,
                 start=list(gamma=2, delta=4))
  oldfit <- ippm(X,
              ~Z + offset(log(f)),
              covariates=list(Z=Z, f=f),
              iScore=Dlogf,
              start=list(gamma=1, delta=1),
              nd=nd)
  oldfit2 <- update(oldfit, . ~ . + I(Z^2))
  oldfit0 <- update(oldfit,
                    . ~ . - Z,
                    start=list(gamma=2, delta=4))
  ## again with logistic
  fitlo2 <- update(fitlo, . ~ . + I(Z^2))
  fitlo0 <- update(fitlo,
                   . ~ . - Z,
                   start=list(gamma=2, delta=4))
  oldfitlo <- ippm(X,
                   ~Z + offset(log(f)),
                   method="logi",
                   covariates=list(Z=Z, f=f),
                   start=list(gamma=1, delta=1),
                   nd=nd)
  oldfitlo2 <- update(oldfitlo, . ~ . + I(Z^2))
  oldfitlo0 <- update(oldfitlo,
                      . ~ . - Z,
                      start=list(gamma=2, delta=4))
  ## anova.ppm including ippm objects
  fit0 <- update(fit, . ~ Z)
  fit0lo <- update(fitlo, . ~ Z)
  A <- anova(fit0, fit)
  Alo <- anova(fit0lo, fitlo)
})
}
