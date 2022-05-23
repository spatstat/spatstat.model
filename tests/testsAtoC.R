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
#'  tests/aucroc.R
#'
#'  AUC and ROC code
#'
#'  $Revision: 1.6 $ $Date: 2020/11/02 06:26:45 $

local({
  if(FULLTEST) {
    fit <- kppm(redwood ~ I(y-x))
    a <- roc(fit)
    b <- auc(fit)
    fet <- ppm(amacrine~x+y+marks)
    d <- roc(fet)
    e <- auc(fet)
  }
})
## tests/cdf.test.R


local({
  NSIM <- 9
  op <- spatstat.options(ndummy.min=16, npixel=32)
  if(FULLTEST) {
    ## Monte Carlo test for Gibbs model
    fit <- ppm(cells ~ 1, Strauss(0.07))
    cdf.test(fit, "x", nsim=NSIM)

    ## cdf.test.slrm
    fut <- slrm(japanesepines ~ x + y)
    Z <- distmap(japanesepines)
    cdf.test(fut, Z)
  }
  reset.spatstat.options()
})


#'
#'   tests/contrib.R
#'
#'   Tests for user-contributed code in spatstat
#'
#'   $Revision: 1.4 $  $Date: 2021/04/17 02:32:24 $

local({
  #' Jinhom
  #' Marie-Colette van Lieshout and Ottmar Cronie
  X <- redwood3
  if(FULLTEST) {
    fit <- ppm(X ~ polynom(x,y,2))
  } else {
    X <- X[c(TRUE,FALSE)]
    spatstat.options(npixel=32, ndummy.min=16)
    fit <- ppm(X ~ x)
  }
  lam <- predict(fit)
  lamX <- fitted(fit, dataonly=TRUE)
  lmin <- 0.9 * min(lam)
  g1 <- Ginhom(X, lambda=fit, update=TRUE)
  if(FULLTEST) {
    g2 <- Ginhom(X, lambda=fit, update=FALSE, lmin = lmin)
    g3 <- Ginhom(X, lambda=lam,  lmin=lmin)
    g4 <- Ginhom(X, lambda=lamX, lmin=lmin)
  }
  if(ALWAYS) {
    f2 <- Finhom(X, lambda=fit, update=FALSE)
  }
  if(FULLTEST) {
    f1 <- Finhom(X, lambda=fit, update=TRUE)
    f3 <- Finhom(X, lambda=lam,  lmin=lmin)
  }
  if(!FULLTEST) reset.spatstat.options()
})
