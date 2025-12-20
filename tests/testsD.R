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
#'
#'    tests/deltasuffstat.R
#' 
#'    Explicit tests of 'deltasuffstat'
#' 
#' $Revision: 1.4 $ $Date: 2021/01/22 08:08:48 $

if(!FULLTEST)
  spatstat.options(npixel=32, ndummy.min=16)

if(ALWAYS) {  # depends on C code
local({
  
  disagree <- function(x, y, tol=1e-7) { !is.null(x) && !is.null(y) && max(abs(x-y)) > tol }

  flydelta <- function(model, modelname="") {
    ## Check execution of different algorithms for 'deltasuffstat' 
    dSS <- deltasuffstat(model, sparseOK=TRUE)
    dBS <- deltasuffstat(model, sparseOK=TRUE,  use.special=FALSE, force=TRUE)
    dBF <- deltasuffstat(model, sparseOK=FALSE, use.special=FALSE, force=TRUE)
    ## Compare results
    if(disagree(dBS, dSS))
      stop(paste(modelname, "model: Brute force algorithm disagrees with special algorithm"))
    if(disagree(dBF, dBS))
      stop(paste(modelname, "model: Sparse and full versions of brute force algorithm disagree"))
    return(invisible(NULL))
  }

  modelS <- ppm(cells ~ x, Strauss(0.13), nd=10)
  flydelta(modelS, "Strauss")

  antsub <- ants[c(FALSE,TRUE,FALSE)]
  rmat <- matrix(c(130, 90, 90, 60), 2, 2)
  
  modelM <- ppm(antsub ~ 1, MultiStrauss(rmat), nd=16)
  flydelta(modelM, "MultiStrauss")
                
  modelA <- ppm(antsub ~ 1, HierStrauss(rmat, archy=c(2,1)), nd=16)
  flydelta(modelA, "HierStrauss")
})

}

reset.spatstat.options()
#'
#'  tests/density.R
#'
#'  Test behaviour of density() methods,
#'                    relrisk(), Smooth()
#'                    and inhomogeneous summary functions
#'                    and idw, adaptive.density, intensity
#'                    and SpatialMedian, SpatialQuantile
#'
#'  $Revision: 1.70 $  $Date: 2025/07/27 07:21:08 $
#'

if(!FULLTEST)
  spatstat.options(npixel=32, ndummy.min=16)


local({
  ## likewise 'relrisk.ppm'
  fit <- ppm(ants ~ x)
  rants <- function(..., model=fit) {
    a <- relrisk(model, sigma=100, se=TRUE, ...)
    return(TRUE)
  }
  if(ALWAYS) {
    rants()
    rants(diggle=TRUE)
    rants(edge=FALSE)
    rants(at="points")
    rants(casecontrol=FALSE)
    rants(relative=TRUE)
  }
  if(FULLTEST) {
    rants(diggle=TRUE, at="points")
    rants(edge=FALSE, at="points")
    rants(casecontrol=FALSE, relative=TRUE)
    rants(casecontrol=FALSE,at="points")
    rants(relative=TRUE,at="points")
    rants(casecontrol=FALSE, relative=TRUE,at="points")
    rants(relative=TRUE, control="Cataglyphis", case="Messor")
    rants(relative=TRUE, control="Cataglyphis", case="Messor", at="points")
  }
  ## more than 2 types
  fut <- ppm(sporophores ~ x)
  if(ALWAYS) {
    rants(model=fut)
  }
  if(FULLTEST) {
    rants(model=fut, at="points")
    rants(model=fut, relative=TRUE, at="points")
  }
  if(FULLTEST) {
    ## intensity - multitype Poisson
    fitM <- ppm(amacrine ~ 1)
    lamM <- intensity(fitM)
    fitMx <- ppm(amacrine ~ x)
    lamMx <- intensity(fitMx, dimyx=32)
  }
})

reset.spatstat.options()

#'
#'      tests/diagnostique.R
#'
#'  Diagnostic tools such as diagnose.ppm, qqplot.ppm
#'
#'  $Revision: 1.7 $  $Date: 2025/12/19 01:33:20 $
#'

if(FULLTEST) {
local({
  cat("Testing diagnose.ppm on Poisson model...")
  fit <- ppm(cells ~ x)
  diagE <- diagnose.ppm(fit, type="eem")
  diagI <- diagnose.ppm(fit, type="inverse")
  diagP <- diagnose.ppm(fit, type="Pearson")
  plot(diagE, which="all")
  plot(diagI, which="smooth")
  plot(diagP, which="x")
  plot(diagP, which="marks", plot.neg="discrete")
  plot(diagP, which="marks", plot.neg="contour")
  plot(diagP, which="smooth", srange=c(-5,5))
  plot(diagP, which="smooth", plot.smooth="contour")
  plot(diagP, which="smooth", plot.smooth="image")

  cat("Testing diagnose.ppm on Strauss model...")
  fitS <- ppm(cells ~ x, Strauss(0.08))
  diagES <- diagnose.ppm(fitS, type="eem", clip=FALSE)
  diagIS <- diagnose.ppm(fitS, type="inverse", clip=FALSE)
  diagPS <- diagnose.ppm(fitS, type="Pearson", clip=FALSE)
  plot(diagES, which="marks", plot.neg="imagecontour")
  plot(diagPS, which="marks", plot.neg="discrete")
  plot(diagPS, which="marks", plot.neg="contour")
  plot(diagPS, which="smooth", plot.smooth="image")
  plot(diagPS, which="smooth", plot.smooth="contour")
  plot(diagPS, which="smooth", plot.smooth="persp")
  
  #' infinite reach, not border-corrected
  cat("Testing diagnose.ppm on Softcore model...")
  fut <- ppm(cells ~ x, Softcore(0.5), correction="isotropic")
  diagnose.ppm(fut)

  #' 
  cat("Testing diagnose.ppm with cumulative=FALSE...")
  diagPX <- diagnose.ppm(fit, type="Pearson", cumulative=FALSE)
  plot(diagPX, which="y")

  #' simulation based
  cat("Testing qqplot.ppm...")
  e <- envelope(cells, nsim=4, savepatterns=TRUE, savefuns=TRUE)
  Plist <- rpoispp(40, nsim=5)

  qf <- qqplot.ppm(fit, nsim=4, expr=e, plot.it=FALSE)
  print(qf)
  qp <- qqplot.ppm(fit, nsim=5, expr=Plist, fast=FALSE)
  print(qp)
  qp <- qqplot.ppm(fit, nsim=5, expr=expression(rpoispp(40)), plot.it=FALSE)
  print(qp)
  qg <- qqplot.ppm(fit, nsim=5, style="classical", plot.it=FALSE)
  print(qg)
  
  #' lurking.ppm
  cat("Testing lurking.ppm...")
  fitx <- ppm(japanesepines ~ x)
  #' covariate is name of coordinate (and has a unit name)
  lurking(fitx, "y")
  #' covariate is numeric vector
  yvals <- coords(as.ppp(quad.ppm(fitx)))[,"y"]
  lurking(fitx, yvals)
  #' covariate is function(x,y) and has a unit name
  D <- distfun(runifpoint(ex=japanesepines))
  lurking(fitx, D)
  #' covariate is stored but is not used in model
  Z <- as.im(function(x,y){ x+y }, Window(japanesepines))
  fitxx <- ppm(japanesepines ~ x, data=solist(Zed=Z), allcovar=TRUE)
  lurking(fitxx, expression(Zed))
  #' envelope is a ppplist; length < nsim default; glmdata=NULL
  fit <- ppm(japanesepines ~ 1)
  stuff <- lurking(fit, expression(x), envelope=Plist, plot.sd=FALSE)
  #' plot.lurk
  plot(stuff, shade=NULL)
  #' problem raised by Gabriela Calana Somoza
  A <- cells
  attr(Window(A), "vvv") <- 42
  fut <- ppm(A ~ x)
  lurking(fut, expression(x), type="raw", envelope=TRUE, nsim=3)
})
}

#'
#'    tests/deepeepee.R
#'
#'    Tests for determinantal point process models
#' 
#'    $Revision: 1.9 $ $Date: 2022/04/24 09:14:46 $

local({
  if(ALWAYS) {
    #' simulate.dppm
    jpines <- residualspaper$Fig1
    fit <- dppm(jpines ~ 1, dppGauss)
    set.seed(10981)
    simulate(fit, W=square(5))
  }
  if(FULLTEST) {
    #' simulate.detpointprocfamily - code blocks
    model <- dppGauss(lambda=100, alpha=.05, d=2)
    simulate(model, seed=1999, correction="border")
    u <- is.stationary(model)
    #' other methods for dppm
    kay <- Kmodel(fit)
    gee <- pcfmodel(fit)
    lam <- intensity(fit)
    arr <- reach(fit)
    pah <- parameters(fit)
    #' a user bug report - matrix dimension error
    set.seed(256)
    dat <- simulate( dppGauss(lambda = 8.5, alpha = 0.1, d = 2), nsim = 1)
  }
  if(FULLTEST) {
    ## cover print.summary.dppm
    jpines <- japanesepines[c(TRUE,FALSE,FALSE,FALSE)]
    print(summary(dppm(jpines ~ 1, dppGauss)))
    print(summary(dppm(jpines ~ 1, dppGauss, method="c")))
    print(summary(dppm(jpines ~ 1, dppGauss, method="p")))
    print(summary(dppm(jpines ~ 1, dppGauss, method="a")))
  }
  #' dppeigen code blocks
  if(ALWAYS) {
    mod <- dppMatern(lambda=2, alpha=0.01, nu=1, d=2)
    uT <- dppeigen(mod, trunc=1.1,  Wscale=c(1,1), stationary=TRUE)
  }
  if(FULLTEST) {
    uF <- dppeigen(mod, trunc=1.1,  Wscale=c(1,1), stationary=FALSE)
    vT <- dppeigen(mod, trunc=0.98, Wscale=c(1,1), stationary=TRUE)
    vF <- dppeigen(mod, trunc=0.98, Wscale=c(1,1), stationary=FALSE)
  }
})
