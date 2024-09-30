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
# tests/NAinCov.R
#
# Testing the response to the presence of NA's in covariates
#
# $Revision: 1.9 $ $Date: 2023/11/05 01:45:36 $

if(FULLTEST) {
local({
  X <- runifpoint(42)
  Y <- as.im(function(x,y) { x+y }, owin())
  Y[owin(c(0.2,0.4),c(0.2,0.4))] <- NA
  ## ..... ppm ......................................
  ## fit model: should produce a warning but no failure
  misfit <- ppm(X ~ Y)
  ## prediction 
  Z <- predict(misfit, type="trend", se=TRUE)
  ## covariance matrix: all should be silent
  v <- vcov(misfit)
  ss <- vcov(misfit, what="internals")
  ## ..... kppm ......................................
  ## should produce warnings but no failures
  misfit <- kppm(X ~Y)
  V <- predict(misfit, type="trend", se=TRUE)
  refit <- improve.kppm(misfit, dimyx=20)
})
}
