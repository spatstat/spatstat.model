# eem.R
#
# Computes the Stoyan-Grabarnik "exponential energy weights" 
#
# $Revision: 1.9 $ $Date: 2026/01/21 06:26:39 $
#

eem <- function(fit, ...) {
  UseMethod("eem")
}

eem.ppm <- function(fit, check=TRUE, ...) {
  verifyclass(fit, "ppm")
  lambda <- fitted.ppm(fit, dataonly=TRUE, check=check)
  eemarks <- 1/lambda
  attr(eemarks, "type") <- "eem"
  attr(eemarks, "typename") <- "exponential energy marks"
  return(eemarks)
}

## Note: class 'exactppm' is defined in spatstat.explore
## but eem() is a generic defined in spatstat.model

eem.exactppm <- function(fit, ...) {
  verifyclass(fit, "exactppm")
  ## lambda <- fitted(fit)
  lambda <- predict(fit, locations=fit$X)
  eemarks <- 1/lambda
  attr(eemarks, "type") <- "eem"
  attr(eemarks, "typename") <- "exponential energy marks"
  return(eemarks)
}

eem.slrm <- function(fit, check=TRUE, ...) {
  verifyclass(fit, "slrm")
  Y <- response(fit)
  lambdaY <- predict(fit, type="intensity")[Y, drop=FALSE]
  eemarks <- 1/lambdaY
  attr(eemarks, "type") <- "eem"
  attr(eemarks, "typename") <- "exponential energy marks"
  return(eemarks)
}

