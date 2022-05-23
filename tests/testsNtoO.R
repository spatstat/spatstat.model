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
# $Revision: 1.8 $ $Date: 2022/05/22 08:37:38 $

if(FULLTEST) {
local({
  X <- runifpoint(42)
  Y <- as.im(function(x,y) { x+y }, owin())
  Y[owin(c(0.2,0.4),c(0.2,0.4))] <- NA
  # fit model: should produce a warning but no failure
  misfit <- ppm(X ~Y, covariates=list(Y=Y))
  # prediction 
  Z <- predict(misfit, type="trend", se=TRUE)
  # covariance matrix: all should be silent
  v <- vcov(misfit)
  ss <- vcov(misfit, what="internals")
})
}
