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
#'   spatstat.model/tests/sidefx.R
#'
#' Test whether plot(do.plot=FALSE) has no side effects on graphics system
#'
#'  $Revision: 1.4 $  $Date: 2025/12/22 08:29:32 $

local({
  if(FULLTEST) {
    ## test whether a graphics device has been started
    deviceExists <- function() { length(dev.list()) != 0 }
    ## check whether executing 'expr' causes creation of a graphics device
    chk <- function(expr) {
      ename <- sQuote(deparse(substitute(expr)))
      if(deviceExists()) {
        ## try switching off the graphics
        graphics.off()
        if(deviceExists()) {
          ## Dang
          warning(paste("Cannot check", ename, 
                        "as a graphics device already exists"),
                  call.=FALSE)
          return(FALSE)
        }
      }
      eval(expr)
      if(deviceExists()) {
        stop(paste("Evaluating", ename, 
                   "caused a graphics device to be started"),
             call.=FALSE)
      }
      return(TRUE)
    }
    



    fit <- ppm(redwood ~ x)
    ## class 'msr'
    chk(plot(residuals(fit, type="raw"), do.plot=FALSE))
    ## class 'leverage.ppm'
    chk(plot(leverage(fit), do.plot=FALSE))
    ## class 'influence.ppm'
    chk(plot(influence(fit), do.plot=FALSE))


  }
})
#
# tests/slrm.R
#
# $Revision: 1.3 $ $Date: 2020/05/01 09:59:59 $
#
# Test slrm fitting and prediction when there are NA's
#

if(ALWAYS) {
local({
  X <- copper$SouthPoints
  W <- owin(poly=list(x=c(0,35,35,1),y=c(1,1,150,150)))
  Y <- X[W]
  fit <- slrm(Y ~ x+y)
  pred <- predict(fit)
  extractAIC(fit)
  fitx <- update(fit, . ~ x)
  simulate(fitx, seed=42)
  if(FULLTEST) {
    unitname(fitx)
    unitname(fitx) <- "km"

    mur <- solapply(murchison,rescale, 1000, "km")
    mur$dfault <- distfun(mur$faults)
    fut <- slrm(gold ~ dfault, data=mur, splitby="greenstone")
    A <- model.images(fut)
  }
})
}


#
#   tests/step.R
#
#   $Revision: 1.5 $  $Date: 2020/05/01 09:59:59 $
#
# test for step() operation
#
if(FULLTEST) {
local({
  Z <- as.im(function(x,y){ x^3 - y^2 }, nztrees$window)
  fitP <- ppm(nztrees ~x+y+Z, covariates=list(Z=Z))
  step(fitP)
  fitS <- update(fitP, Strauss(7))
  step(fitS)
  fitM <- ppm(amacrine ~ marks*(x+y),
              MultiStrauss(types=levels(marks(amacrine)), radii=matrix(0.04, 2, 2)))
  step(fitM)
})
}



