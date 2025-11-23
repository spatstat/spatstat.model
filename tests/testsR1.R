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
#'  tests/resid.R
#'
#'  Stuff related to residuals and residual diagnostics
#'           including residual summary functions
#'
#'   $Revision: 1.8 $  $Date: 2025/11/23 09:41:53 $
#'

local({
  fit <- ppm(cells ~x, Strauss(r=0.15))
  rr <- residuals(fit, quad=quadscheme(cells, nd=128))
  diagnose.ppm(fit, cumulative=FALSE, type="pearson")

  if(FULLTEST) {
    diagnose.ppm(fit, cumulative=FALSE)

    fitoff <- ppm(cells ~ sin(x) + offset(y))
    plot(a <- parres(fitoff, "x"))
    plot(b <- parres(fitoff, "y"))
    print(a)
    print(b)
  
    d <- diagnose.ppm(fit, which="marks")
    plot(d, plot.neg="discrete")
    plot(d, plot.neg="imagecontour")

    d <- diagnose.ppm(fit, type="pearson", which="smooth")
    plot(d, plot.smooth="image")
    plot(d, plot.smooth="contour")
    plot(d, plot.smooth="imagecontour")
  
    d <- diagnose.ppm(fit, type="pearson", which="x")
    plot(d)
    d <- diagnose.ppm(fit, type="pearson", which="y")
    plot(d)
  
    diagnose.ppm(fit, type="pearson", which="x", cumulative=FALSE)
    diagnose.ppm(fit, type="pearson", which="x", cumulative=FALSE)
    diagnose.ppm(fit, type="raw", plot.neg="discrete", plot.smooth="image")
    diagnose.ppm(fit, type="pearson", plot.neg="contour", plot.smooth="contour")

    diagnose.ppm(fitoff, type="raw", which="smooth", plot.smooth="persp")
    diagnose.ppm(fitoff, type="pearson", plot.neg="imagecontour")

    plot(Frame(letterR), main="")
    ploterodewin(letterR, erosion(letterR, 0.05), main="jeans")
    W <- as.mask(letterR)
    plot(Frame(W), main="")
    ploterodewin(W, erosion(W, 0.05), main="JeAnS")

    #' entangled terms in model
    U <- as.im(1, owin())
    Z <- as.im(function(x,y) x, owin())
    X <- runifpoint(40)
    fut <- ppm(X ~ Z:U)
    a <- parres(fut, "Z")
    futoff <- ppm(X ~ offset(Z*U))
    a <- parres(futoff, "Z")

    #' residual summary functions
    pt <- psst(cells, interaction=Strauss(0.1), fun=nndcumfun)

    #' code blocks in residualMeasure
    A <- residualMeasure(cells, -42)
    B <- residualMeasure(cells, function(x,y) { rep(-42, length(x)) })
  }
})



##
## tests/rhohat.R
##
## Test all combinations of options for rhohatCalc
##
## $Revision: 1.6 $ $Date: 2022/05/22 08:03:48 $

local({
  if(FULLTEST) {
    X <-  rpoispp(function(x,y){exp(3+3*x)})
    Z <- as.im(function(x,y) { x }, Window(X))
    f <- funxy(function(x,y) { y + 1 }, Window(X))
    
    
    ## rhohat.ppm
    fit <- ppm(X ~x)
    rhofitA <- rhohat(fit, "x")
    rhofitB <- rhohat(fit, "x", method="reweight")
    rhofitC <- rhohat(fit, "x", method="transform")
    rhofitD <- rhohat(fit, Z)
    rhofitD <- rhohat(fit, Z, positiveCI=TRUE)
    lam <- predict(fit)


    ## Horvitz-Thompson
    rhofitAH <- rhohat(fit, "x", horvitz=TRUE)
    rhofitBH <- rhohat(fit, "x", method="reweight", horvitz=TRUE)
    rhofitCH <- rhohat(fit, "x", method="transform", horvitz=TRUE)

    r2myx <- rho2hat(fit, "y", "x")
    r2myxw <- rho2hat(fit, "y", "x", method="reweight")
    plot(r2myx)
    plot(r2myxw)
    print(r2myxw)
    predict(r2myxw)
    predict(r2myxw, relative=TRUE)
  }
})
