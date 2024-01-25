#'
#'  getglmfit.R
#'
#'  Generics to extract fitted glm/gam model objects
#'  from point process model objects
#'
#'  $Revision: 1.3 $ $Date: 2024/01/25 09:17:46 $
#'

getglmdata <- function(object, ...) {
  UseMethod("getglmdata")
}

getglmfit <- function(object, ...) {
  UseMethod("getglmfit")
}

getglmsubset <- function(object, ...) {
  UseMethod("getglmsubset")
}

hasglmfit <- function(object) {
  UseMethod("hasglmfit")
}


