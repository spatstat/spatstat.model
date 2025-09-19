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
# tests/kppm.R
#
# $Revision: 1.42 $ $Date: 2025/09/19 02:53:00 $
#
# Test functionality of kppm that once depended on RandomFields
# Test update.kppm for old style kppm objects

if(!FULLTEST)
  spatstat.options(npixel=32, ndummy.min=16)

local({

 fit <- kppm(redwood ~1, "Thomas") # sic
 fitx <- kppm(redwood ~x, "Thomas", verbose=TRUE)
 if(FULLTEST) {
   fitx <- update(fit, ~ . + x)
   fitM <- update(fit, clusters="MatClust")
   fitC <- update(fit, cells)
   fitCx <- update(fit, cells ~ x)
   #'
   Wsub <- owin(c(0, 0.5), c(-0.5, 0))
   Zsub <- (bdist.pixels(Window(redwood)) > 0.1)
   fitWsub <- kppm(redwood ~1, "Thomas", subset=Wsub)
   fitZsub <- kppm(redwood ~1, "Thomas", subset=Zsub)
   fitWsub
 
   #' various methods
   ff <- as.fv(fitx)
   uu <- unitname(fitx)
   unitname(fitCx) <- "furlong"
   mo <- model.images(fitCx)
   p <- psib(fit)
   px <- psib(fitx)
 }
 if(ALWAYS) {
   Y <- simulate(fitx, seed=42, saveLambda=TRUE)[[1]]
   stopifnot(is.im(attr(Y, "Lambda")))
 }

 if(FULLTEST) {
   #' vcov.kppm different algorithms
   vc  <- vcov(fitx)
   vc2 <- vcov(fitx, fast=TRUE)
   vc3 <- vcov(fitx, fast=TRUE, splitup=TRUE)
   vc4 <- vcov(fitx,            splitup=TRUE)

   ## other code blocks
   a <- varcount(fitx, function(x,y){x+1}) # always positive
   a <- varcount(fitx, function(x,y){y-1}) # always negative
   a <- varcount(fitx, function(x,y){x+y}) # positive or negative

   #' improve.kppm
   fitI <- update(fit, improve.type="quasi")
   fitxI <- update(fitx, improve.type="quasi")
   fitxIs <- update(fitx, improve.type="quasi", fast=FALSE) 
   #' vcov.kppm
   vcI <- vcov(fitxI)
 }

  ## plot.kppm including predict.kppm
 if(ALWAYS) {
   fitMC <- kppm(redwood ~ x, "Thomas")
   plot(fitMC)
 }
 if(FULLTEST) {
   fitCL <- kppm(redwood ~ x, "Thomas", method="c")
   fitPA <- kppm(redwood ~ x, "Thomas", method="p")
   plot(fitCL)
   plot(fitPA)

   ## fit with composite likelihood method [thanks to Abdollah Jalilian]
   fut <- kppm(redwood ~ x, "VarGamma", method="clik2", nu.ker=-3/8)
   kfut <- as.fv(fut)
 }
 
 if(ALWAYS) {
   fit0 <- kppm(redwood ~1, "LGCP")
   is.poisson(fit0)
   Y0 <- simulate(fit0, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y0))
   stopifnot(is.im(attr(Y0, "Lambda")))
   p0 <- psib(fit0) # issues a warning

   if(FULLTEST) {
     ## fit LGCP using K function: slow
     fit1 <- kppm(redwood ~x, "LGCP",
                  covmodel=list(model="matern", nu=0.3),
                  control=list(maxit=3))
     Y1 <- simulate(fit1, saveLambda=TRUE)[[1]]
     stopifnot(is.ppp(Y1))
     stopifnot(is.im(attr(Y1, "Lambda")))
   }

   ## fit LGCP using pcf
   fit1p <- kppm(redwood ~x, "LGCP",
                 covmodel=list(model="matern", nu=0.3),
                 statistic="pcf")
   Y1p <- simulate(fit1p, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1p))
   stopifnot(is.im(attr(Y1p, "Lambda")))

   ## .. and using different fitting methods
   if(FULLTEST) {
     fit1pClik <- update(fit1p, method="clik")
     fit1pPalm <- update(fit1p, method="palm")
   }

   ## shortcut evaluation of pcf
   ##  (the code being tested is in spatstat.random::clusterinfo.R)
   if(FULLTEST) {
     putSpatstatVariable("RFshortcut", TRUE)
     fitGshort <- kppm(redwood ~ 1, "LGCP", covmodel=list(model="gauss"))
     fitSshort <- kppm(redwood ~ 1, "LGCP", covmodel=list(model="stable", alpha=1))
     putSpatstatVariable("RFshortcut", FALSE)
     fitGlong <-  kppm(redwood ~ 1, "LGCP", covmodel=list(model="gauss"))
     fitSlong <-  kppm(redwood ~ 1, "LGCP", covmodel=list(model="stable", alpha=1))
     discrepG <- unlist(parameters(fitGshort)) - unlist(parameters(fitGlong))
     discrepS <- unlist(parameters(fitSshort)) - unlist(parameters(fitSlong))
     print(discrepG)
     print(discrepS)
     if(max(abs(discrepG) > 0.01))
       stop("Discrepancy in short-cut fitting of Gaussian LGCP")
     if(max(abs(discrepS) > 0.01))
       stop("Discrepancy in short-cut fitting of stable LGCP")
   }
   
   ## image covariate (a different code block) 
   xx <- as.im(function(x,y) x, Window(redwood))
   fit1xx <- update(fit1p, . ~ xx, data=solist(xx=xx))
   Y1xx <- simulate(fit1xx, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1xx))
   stopifnot(is.im(attr(Y1xx, "Lambda")))
   if(FULLTEST) {
     fit1xxVG <- update(fit1xx, clusters="VarGamma", nu=-1/4)
     Y1xxVG <- simulate(fit1xxVG, saveLambda=TRUE)[[1]]
     stopifnot(is.ppp(Y1xxVG))
     stopifnot(is.im(attr(Y1xxVG, "Lambda")))
   }
   fit1xxLG <- update(fit1xx, clusters="LGCP",
                      covmodel=list(model="matern", nu=0.3),
                      statistic="pcf")
   Y1xxLG <- simulate(fit1xxLG, saveLambda=TRUE, drop=TRUE)
   stopifnot(is.ppp(Y1xxLG))
   stopifnot(is.im(attr(Y1xxLG, "Lambda")))
   
   # ... and Abdollah's code
   if(FULLTEST) {
     fit2 <- kppm(redwood ~x, cluster="Cauchy", statistic="K")
     Y2 <- simulate(fit2, saveLambda=TRUE)[[1]]
     stopifnot(is.ppp(Y2))
     stopifnot(is.im(attr(Y2, "Lambda")))
   }
 }

})

if(FULLTEST) {
local({
  #' various code blocks
  fut <- kppm(redwood, ~x)
  fet <- update(fut, redwood3)
  fot <- update(fut, trend=~y)
  fit <- kppm(redwoodfull ~ x)
  YR <- simulate(fit, window=redwoodfull.extra$regionII, saveLambda=TRUE)[[1]]
  stopifnot(is.im(attr(YR, "Lambda")))
  gut <- improve.kppm(fit, type="wclik1")
  gut <- improve.kppm(fit, vcov=TRUE, fast.vcov=TRUE, save.internals=TRUE)
  hut <- kppm(redwood ~ x, method="clik", weightfun=NULL)
  hut <- kppm(redwood ~ x, method="palm", weightfun=NULL)
  mut <- kppm(redwood)
  nut <- update(mut, YR)
})
}

if(FULLTEST) {
local({
  #' minimum contrast code
  K <- Kest(redwood)
  a <- matclust.estK(K)
  a <- thomas.estK(K)
  a <- cauchy.estK(K)
  a <- vargamma.estK(K)
  a <- lgcp.estK(K)

  print(a)
  u <- unitname(a)
  
  g <- pcf(redwood)
  a <- matclust.estpcf(g)
  a <- thomas.estpcf(g)
  a <- cauchy.estpcf(g)
  a <- vargamma.estpcf(g)
  a <- lgcp.estpcf(g)

  #' auxiliary functions
  b <- resolve.vargamma.shape(nu.pcf=1.5)
  Z <- clusterfield("Thomas", kappa=1, scale=0.2)
  
  aa <- NULL
  aa <- accumulateStatus(simpleMessage("Woof"), aa)
  aa <- accumulateStatus(simpleMessage("Sit"),  aa)
  aa <- accumulateStatus(simpleMessage("Woof"), aa)
  printStatusList(aa)

  RMIN <- 0.01
  fit <- kppm(redwood ~ 1, ctrl=list(rmin=RMIN,q=1/2))
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("kppm did not handle parameter 'rmin' in argument 'ctrl' ")
  fit <- kppm(redwood ~ 1, ctrl=list(rmin=0,q=1/2), rmin=RMIN)
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("kppm did not handle parameter 'rmin' in argument 'ctrl'")

  RMIN <- 2
  fit <- dppm(swedishpines~1, dppGauss(), ctrl=list(rmin=RMIN,q=1))
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("dppm did not handle parameter 'rmin' in argument 'ctrl'")
  fit <- dppm(swedishpines~1, dppGauss(), ctrl=list(rmin=0,q=1), rmin=RMIN)
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("dppm did not handle argument 'rmin'")
})
}



if(FULLTEST) {
local({
  #' cover a few code blocks
  fut <- kppm(redwood ~ x, method="clik")
  print(summary(fut))
  a <- residuals(fut)
  fut2 <- kppm(redwood ~ x, "LGCP", method="palm")
  print(summary(fut2))
  b <- residuals(fut2)
  #'
  po <- ppm(redwood ~ 1)
  A <- kppmComLik(redwood, Xname="redwood", po=po, clusters="Thomas",
                  statistic="pcf", statargs=list(), control=list(),
                  weightfun=NULL, rmax=0.1)
  A <- kppmPalmLik(redwood, Xname="redwood", po=po, clusters="Thomas",
                   statistic="pcf", statargs=list(), control=list(),
                   weightfun=NULL, rmax=0.1)
})
}

local({
  if(FULLTEST) {
    ## Simulate a model on a completely different window
    wleft <- owin(c(0, 300), c(0, 500))
    wright <- owin(c(700, 1000), c(0, 500))
    bleft <- bei[wleft]
    dleft <- solapply(bei.extra, "[", i=wleft)
    dright <- solapply(bei.extra, "[", i=wright)
    mleft <- kppm(bleft ~ elev, data=dleft)
    ## simulate left model on right window
    X <- simulate(mleft, w=wright, covariates=dright)
  }
})

reset.spatstat.options()
  
#'
#'   tests/Kfuns.R
#'
#'   Various K and L functions and pcf
#'
#'   $Revision: 1.45 $  $Date: 2025/03/15 11:29:33 $
#'
#'   Assumes 'EveryStart.R' was run

if(FULLTEST) {
  Cells <- cells
  Amacrine <- amacrine
  Redwood <- redwood
} else {
  ## reduce numbers of data + dummy points
  spatstat.options(npixel=32, ndummy.min=16)
  Cells <- cells[c(FALSE, TRUE)]
  Amacrine <- amacrine[c(FALSE, TRUE)]
  Redwood <- redwood[c(FALSE, TRUE)]
}


if(FULLTEST) {
  local({
    #' code blocks using fitted model to determine intensity
    #' Kinhom
    X <- rpoispp(function(x,y) { 100 * x }, 100, square(1))
    fut <- ppm(X ~ x)
    Kio <- Kinhom(X, fut, update=FALSE)
    Kiu <- Kinhom(X, fut, update=TRUE, diagonal=FALSE)
    fit <- ppm(Amacrine ~ marks)
    #' lohboot Linhom
    Zred <- predict(ppm(Redwood ~ x+y))
    Lred <- lohboot(Redwood, Linhom, lambda=Zred)
    #' Kmulti.inhom
    K1 <- Kcross.inhom(Amacrine, lambdaX=fit)
    On <- split(Amacrine)$on
    Off <- split(Amacrine)$off
    K4 <- Kcross.inhom(Amacrine, lambdaI=ppm(On), lambdaJ=ppm(Off))
    #' local K functions
    fut <- ppm(swedishpines ~ polynom(x,y,2))
    Z <- predict(fut)
    Lam <- fitted(fut, dataonly=TRUE)
    a <- localLinhom(swedishpines, lambda=fut)
    a <- localLinhom(swedishpines, lambda=Z)
    a <- localLinhom(swedishpines, lambda=Lam)
    a <- localLinhom(swedishpines, lambda=Z, correction="none")
    a <- localLinhom(swedishpines, lambda=Z, correction="translate")
    #' local cross K functions
    fat <- ppm(Amacrine ~ x * marks)
    Zed <- predict(fat)
    Lum <- fitted(fat, dataonly=TRUE)
    moff <- (marks(Amacrine) == "off")
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Zed)
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Lum)
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=fat)
    a <- localLcross.inhom(Amacrine, from="off", to="on",
                           lambdaFrom=Lum[moff], lambdaTo=Lum[!moff])
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Zed,
                           correction="none")
    a <- localLcross.inhom(Amacrine, from="off", to="on", lambdaX=Zed,
                           correction="translate")
    #' cases of resolve.lambdacross
    h <- resolve.lambdacross(Amacrine, moff, !moff, lambdaX=fat)
    h <- resolve.lambdacross(Amacrine, moff, !moff, lambdaX=fat, update=FALSE)
    h <- resolve.lambdacross(Amacrine, moff, !moff,
                              lambdaI=fat, lambdaJ=fat)
    h <- resolve.lambdacross(Amacrine, moff, !moff,
                              lambdaI=fat, lambdaJ=fat,
                              update=FALSE)
    #' lohboot
    b <- lohboot(Amacrine, Lcross.inhom, from="off", to="on", lambdaX=Zed)
    b <- lohboot(Amacrine, Lcross.inhom, from="off", to="on", lambdaX=Lum)
    b <- lohboot(Amacrine, Lcross.inhom, from="off", to="on", lambdaX=fat)
    b <- lohboot(Amacrine, Lcross.inhom, from="off", to="on", 
                 lambdaFrom=Lum[moff], lambdaTo=Lum[!moff])
    #'
    #'  residual K functions etc
    #'
    rco <- compareFit(Cells, Kcom,
                      interaction=anylist(P=Poisson(), S=Strauss(0.08)),
                      same="trans", different="tcom")
    fit <- ppm(Cells ~ x, Strauss(0.07))
    K <- Kcom(Cells, model=fit, restrict=TRUE)
  })
}

reset.spatstat.options()
