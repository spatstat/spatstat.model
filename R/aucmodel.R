##
##
##  aucmodel.R
##
##  Calculate Area Under ROC curve
##       (in spatstat.model)
##
##
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner
##





auc.kppm <- function(X, ..., subset=NULL) {
  auc(as.ppm(X), ..., subset=subset)
}

  auc.ppm <-
  function(X, ..., subset=NULL) {
    model <- X
    use.roc <- is.multitype(model) || needROC(...)
    if(use.roc) {
      ro <- roc(model, ..., subset=subset)
      aobs <- with(ro, mean(.y))
      atheo <- with(ro, mean(theo))
    } else if(is.ppm(model) && is.stationary(model)) {
      aobs <- atheo <- 1/2
    } else {
      lambda <- predict(model, ..., type="trend")
      X <- if(is.ppm(model)) data.ppm(model) else model$X
      if(!is.null(subset)) {
        #' restrict to subset
        lambda <- lambda[subset, drop=FALSE]
        X <- X[subset]
      }
      lamX <- lambda[X]
      lamW <- lambda[]
      Fl <- ecdf(lamW)
      aobs <- mean(Fl(lamX))
      atheo <- mean(lamW * Fl(lamW))/mean(lamW)
    }
    result <- c(aobs, atheo)
    names(result) <- c("obs", "theo")
    return(result)
  }

auc.slrm <- function(X, ..., subset=NULL) {
  model <- X
  if(is.stationary(model)) {
    aobs <- atheo <- 1/2
  } else {
    ro <- roc(model, ..., subset=subset)
    aobs <- with(ro, mean(.y))
    atheo <- with(ro, mean(theo))
  }
  result <- c(aobs, atheo)
  names(result) <- c("obs", "theo")
  return(result)
}



