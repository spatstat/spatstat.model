#'
#'   response.R
#'
#'   Extract the values of the response, from a fitted model
#' 
#'

response <- function(object) {
  UseMethod("response")
}

response.glm <- response.lm <- function(object) {
  mo <- object$model
  if(is.null(mo)) return(NULL)
  te <- terms(object)
  rn <- attr(te, "response")
  if(is.null(rn)) return(NULL)
  y <- mo[,rn]
  return(y)
}

response.ppm <- function(object) { data.ppm(object) }

response.dppm <- response.kppm <- function(object) { data.ppm(as.ppm(object)) }

response.slrm <- function(object) { object$Data$response }

response.mppm <- function(object) { data.mppm(object) }
