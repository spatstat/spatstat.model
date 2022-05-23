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
## $Revision: 1.5 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({

})
}
