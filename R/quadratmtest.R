#
#   method for 'quadrat.test' for class mppm
#
#   $Revision: 1.9 $   $Date: 2025/11/17 08:59:18 $
#
quadrat.test.mppm <- function(X, ...) {
  Xname <- short.deparse(substitute(X))
  if(!is.poisson.mppm(X))
    stop("Model is not a Poisson point process")
  
  subs <- subfits(X)
  tests <- anylapply(subs, quadrat.test, ..., fitname=Xname)

  df.est <- length(coef(X))
  return(pool.quadrattest(tests, Xname=Xname, df.est=df.est))
}

