#'
#'   ic.kppm.R
#'
#'   Original by Rasmus Waagepetersen, 10 april 2021
#'
#'   Tweaks by Adrian Baddeley

ic <- function(object) { UseMethod("ic") }

ic.ppm <- function(object) {
  loglike <- logLik(object)
  ll <- as.numeric(loglike)
  df <- attr(loglike, "df")
  ## betahat <- coef(object)
  n <- npoints(data.ppm(object))
  cbic <- -2*ll+df*log(n)
  cic  <- -2*ll+df*2
  ## cbic is BIC and cic is AIC in case of Poisson process
  return(list(loglike=ll, cbic=cbic, cic=cic, df=df))
}

ic.kppm <- function(object){
  po <- as.ppm(object)
  loglike <- logLik(po)
  ll <- as.numeric(loglike)
  betahat <- coef(object)
  p <- length(betahat)
  n <- npoints(data.ppm(po))
  co <- vcov(object, what="internals")
  df <- p + sum(diag(as.matrix(co$J.inv %*% co$E))) #compute p_approx
  cbic = -2*loglike+df*log(n)
  cic = -2*loglike+df*2
  cbic <- -2*ll+df*log(n)
  cic  <- -2*ll+df*2
  ## cbic is BIC and cic is AIC in case of Poisson process
  return(list(loglike=ll, cbic=cbic, cic=cic, df=df))
}
