##
## auc.R
##
##  Calculate ROC curve or area under it
##
## $Revision: 1.15 $ $Date: 2022/05/22 04:08:02 $

## Code for roc() and roc.ppp() and internals is moved to spatstat.explore

roc.ppm <- function(X, ...) {
  stopifnot(is.ppm(X))
  model <- X
  lambda <- predict(model, ...)
  Y <- data.ppm(model)
  nullmodel <- ppm(Y)
  result <- rocModel(lambda, nullmodel, ...)
  return(result)
}

roc.kppm <- function(X, ...) {
  stopifnot(is.kppm(X))
  model <- as.ppm(X)
  lambda <- predict(model, ...)
  Y <- data.ppm(model)
  nullmodel <- ppm(Y)
  result <- rocModel(lambda, nullmodel, ...)
  return(result)
}

roc.slrm <- function(X, ...) {
  stopifnot(is.slrm(X))
  model <- X
  lambda <- predict(model, ..., type="probabilities")
  Y <- response(model)
  nullmodel <- slrm(Y ~ 1)
  dont.complain.about(Y)
  result <- rocModel(lambda, nullmodel, ..., lambdatype="probabilities")
  return(result)
}

#    ......................................................

## Code for auc(), auc.ppp() is moved to spatstat.explore

auc.kppm <- function(X, ...) { auc(as.ppm(X), ...) }

auc.ppm <- function(X, ...) {
  model <- X
  if(is.multitype(model)) {
    # cheat
    ro <- roc(model, ...)
    aobs <- with(ro, mean(fobs))
    atheo <- with(ro, mean(ftheo))
  } else if(is.stationary(model)) {
    aobs <- atheo <- 1/2
  } else {
    lambda <- intensity(model)
    Fl <- ecdf(lambda[])
    lambda <- as.im(lambda, Window(model))
    X <- data.ppm(model)
    lamX <- lambda[X]
    aobs <- mean(Fl(lamX))
    atheo <- mean(lambda[] * Fl(lambda[]))/mean(lambda)
  }
  result <- c(aobs, atheo)
  names(result) <- c("obs", "theo")
  return(result)
}

auc.slrm <- function(X, ...) {
  ro <- roc(X, ...)
  result <- with(ro, list(obs=mean(fobs), theo=mean(ftheo)))
  return(unlist(result))
}


