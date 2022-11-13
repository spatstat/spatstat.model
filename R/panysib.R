#'          panysib.R
#'  Probability that a point has any siblings
#'
#'  $Revision: 1.3 $ $Date: 2022/11/13 03:13:18 $
#'
#'  Copyright (c) Adrian Baddeley 2022
#'  GNU Public Licence >= 2.0

panysib <- function(object) {
  stopifnot(is.poissonclusterprocess(object))
  if(is.stationary(object))
    return(1-exp(-object$mu))
  R <- reach(object)
  W <- Window(object)
  if(R > 5 * diameter(Frame(W))) {
    ## treat as stationary, but return image
    value <- 1 - exp(-mean(object$mu))
    result <- as.im(value, W=W)
  } else {
    par <- parameters(object)
    lam <- predict(object, window=grow.rectangle(Frame(W), R))
    mu <- lam/par[["kappa"]]
    h <- clusterkernel(object)
    EM <- blur(mu, kernel=h)
    EM <- eval.im(pmax(0, EM))
    P <- blur(exp(-EM), kernel=h)
    P <- eval.im(pmax(0, P))
    result <- 1 - P[W, drop=FALSE, tight=TRUE]
  }
  return(result)
}

