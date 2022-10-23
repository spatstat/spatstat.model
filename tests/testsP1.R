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
## 
##    tests/percy.R
##
## Tests of Percus-Yevick approximations
##
##    $Revision: 1.3 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({
  fit <- ppm(swedishpines ~1, DiggleGatesStibbard(6))
  K <- Kmodel(fit)
})
}

##
## tests/pixelgripes.R
##     Problems related to pixellation of windows
##
## $Revision: 1.8 $ $Date: 2022/10/23 06:21:10 $

if(FULLTEST) {
  local({
    
  ## From Philipp Hunziker: bug in rNeymanScott (etc)
  ## Create an irregular window
  PM <- matrix(c(1,0,0.5,1,0,0), 3, 2, byrow=TRUE)
  P <- owin(poly=PM)
  ## Generate Matern points
  X <- rMatClust(50, 0.05, 5, win=P)
  ## Some distance function as a covariate
  distorigin <- function(x, y) { sqrt(x^2 + y^2) }
  ## No covariates: works fine
  fit0 <- kppm(X ~ 1, clusters="MatClust")
  Y0 <- simulate(fit0, retry=0)
  ## Covariates: Simulation fails
  fit1 <- kppm(X ~ distorigin, clusters="MatClust")
  Y1 <- simulate(fit1, retry=0)

})
}
