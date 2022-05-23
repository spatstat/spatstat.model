#'
#'    varcount.R
#'
#'   Variance of N(B)
#'
#'  $Revision: 1.15 $  $Date: 2021/08/12 09:28:12 $
#'

varcount <- function(model, B=Window(model), ..., dimyx=NULL) {
  stopifnot(is.owin(B) || is.im(B) || is.function(B))
  if(is.owin(B)) {
    f <- NULL
    lambdaB <- predict(model, locations=B, ngrid=dimyx, type="intensity")
  } else {
    f <- if(is.im(B)) B else as.im(B, W=as.owin(model), ..., dimyx=dimyx)
    B <- as.owin(f)
    lambdaB <- predict(model, locations=B, type="intensity")
  }
  ## important range of distances
  ## need to integrate over [0, R]
  R <- reach(model, epsilon=0.001) 
  if(!isTRUE(is.finite(R))) R <- NULL
  if(!is.null(R)) {
    ## detect very small cluster radius (or very large window)
    diam <- diameter(Frame(B))
    if((R < 0.001 * diam) && (area(erosion(B, R))/area(B) > 0.999)) {
      #' integrate by parts and ignore edge effect
      K <- Kmodel(model)
      if(is.function(K)) {
        excess <- K(diam) - pi * diam^2 
        if(is.null(f)) {
          E <- integral(lambdaB)
          V <- integral(lambdaB^2)
        } else {
          E <- integral(lambdaB * f)
          V <- integral((lambdaB * f)^2)
        }
        v <- E + V * excess
        if(is.finite(v)) 
          return(v)
      }
    }
  } 
  g <- pcfmodel(model)
  if(!is.function(g))
    stop("Pair correlation function is not available")
  v <- varcountEngine(g, B, lambdaB, f, R=R)
  return(v)
}

varcountEngine <- local({

  varcountEngine <- function(g, B, lambdaB, f=1, R=NULL,
                             what=c("variance","excess","pairs","squared")) {
    ## variance = var(N)
    ## excess   = var(N) - E(N)
    ## pairs    = E[ N(N-1) ]
    ## squared  = E[N^2]
    what <- match.arg(what)
    g1 <- function(r) { g(r) - 1 }
    if(missing(f) || is.null(f) || identical(f, 1)) {
      v <- switch(what,
                  variance = integral(lambdaB) + dublin(g1, B, lambdaB, R=R),
                  excess   =                     dublin(g1, B, lambdaB, R=R),
                  pairs    =                     dublin(g,  B, lambdaB, R=R),
                  squared  = integral(lambdaB) + dublin(g,  B, lambdaB, R=R))
    } else if(min(f) >= 0) {
      ## nonnegative integrand
      v <- switch(what,
        variance = integral(lambdaB * f^2) + dublin(g1, B, lambdaB * f, R=R),
        excess   =                           dublin(g1, B, lambdaB * f, R=R),
        pairs    =                           dublin(g,  B, lambdaB * f, R=R),
        squared  = integral(lambdaB * f^2) + dublin(g,  B, lambdaB * f, R=R))
    } else if(max(f) <= 0) {
      ## nonpositive integrand
      v <- switch(what,
        variance=integral(lambdaB * f^2) + dublin(g1, B, lambdaB * (-f), R=R),
        excess  =                          dublin(g1, B, lambdaB * (-f), R=R),
        pairs   =                          dublin(g,  B, lambdaB * (-f), R=R),
        squared =integral(lambdaB * f^2) + dublin(g,  B, lambdaB * (-f), R=R))
    } else {
      ## integrand has both positive and negative parts
      lamfplus <- eval.im(lambdaB * pmax(0, f))
      lamfminus <- eval.im(lambdaB * pmax(0, -f))
      h <- switch(what,
                  variance = g1,
                  excess   = g1,
                  pairs    = g,
                  squared  = g)
      co <- (dublin(h, B, lamfplus, R=R) 
            + dublin(h, B, lamfminus, R=R)
            - dublin(h, B, lamfplus, lamfminus, R=R)
            - dublin(h, B, lamfminus, lamfplus, R=R))
      v <- switch(what,
                  variance = integral(lambdaB * f^2) + co,
                  excess   =                           co,
                  pairs    =                           co,
                  squared  = integral(lambdaB * f^2) + co)
    }
    return(v)
  }

  dublin <- function(h, B, f, f2, R=NULL) {
    ## Double integral
    ##          \int_B \int_B h(|u-v|) f(u) f(v) du dv    
    ## or       \int_B \int_B h(|u-v|) f(u) f2(v) du dv
    ## Assume h, f, f2 are nonnegative
    ## R = reach of model
    dr <- R/100
    if(missing(f2)) {
      ## \int_B \int_B h(|u-v|) f(u) f(v) du dv
      M <- distcdf(B, dW=f, nr=NULL, delta=dr)
      ## integrate 
      a <- integral(f)^2 * as.numeric(stieltjes(h, M))
    } else {
      ## \int_B \int_B h(|u-v|) f(u) f2(v) du dv
      M <- distcdf(B, dW=f, dV=f2, nr=NULL, delta=dr)
      ## integrate 
      a <- integral(f) * integral(f2) * as.numeric(stieltjes(h, M))
    }
    return(a)
  }


  varcountEngine
})


