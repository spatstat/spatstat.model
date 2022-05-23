#'
#'   ppmcovariates.R
#'
#'   Utilities for wrangling covariates in ppm
#'
#'   $Revision: 1.1 $ $Date: 2022/05/20 03:59:14 $


ppmCovariates <- function(model) {
  # generate list of all covariates in ppm (excluding marks)
  stopifnot(is.ppm(model))
  co <- as.list(model$covariates)
  xy <- list(x=function(x,y){x}, y=function(x,y){y})
  coplus <- append(co, xy)
  return(as.anylist(coplus))
}

findCovariate <- function(covname, scope, scopename=NULL) {
  # find the named covariate in the given ppm object or list
  if(is.ppm(scope)) {
    covlist <- ppmCovariates(scope)
    if(missing(scopename)) scopename <- "covariates in model"
  } else if(is.list(scope)) {
    covlist <- scope
  } else stop("scope should be a named list of covariates, or a ppm object")
  if(!(covname %in% names(covlist))) 
    stop(paste("covariate", dQuote(covname), "not found",
               if(!is.null(scopename)) paste("amongst", scopename) else NULL))
  covlist[[covname]]
}

