# eem.R
#
# Computes the Stoyan-Grabarnik "exponential energy weights" 
#
# $Revision: 1.7 $ $Date: 2025/11/24 04:01:55 $
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

eem.exactppm <- function(fit, ...) {
  verifyclass(fit, "exactppm")
  lambda <- fitted(fit)
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

