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



