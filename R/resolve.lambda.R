#' resolve.lambda.R
#'
#' Common code to evaluate the intensity 'lambda' in Kinhom, pcfinhom
#'    (and multitype counterparts)
#'
#'    resolve.lambda
#'    resolve.lambdacross
#'    resolve.reciplambda
#'    validate.weights
#'    updateData (generic)
#'
#' $Revision: 1.18 $ $Date: 2022/08/09 06:46:50 $

resolve.lambda <- function(X, lambda=NULL, ...) {
  UseMethod("resolve.lambda")
}

resolve.lambda.ppp <- function(X, lambda=NULL, ...,
                               sigma=NULL, varcov=NULL,
                               leaveoneout=TRUE, update=TRUE,
                               check=TRUE) {
  dangerous <- "lambda"
  danger <- TRUE
  if(npoints(X) == 0) {
    danger <- FALSE
    lambda <- numeric(0)
  } else if(is.null(lambda)) {
    ## No intensity data provided
    ## Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=leaveoneout)
    lambda <- as.numeric(lambda)
    if(check) validate.weights(lambda, how="density estimation")
    danger <- FALSE
  } else if(is.im(lambda)) {
    lambda <- safelookup(lambda, X)
  } else if(is.function(lambda)) {
    lambda <- lambda(X$x, X$y)
  } else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) {
    check.nvector(lambda, npoints(X), vname="lambda")
  } else if(inherits(lambda, c("ppm", "kppm", "dppm", "slrm"))) {
    ## model provides intensity
    model <- lambda
    if(update) {
      model <- updateData(model, X)
      danger <- FALSE
    }
    if(inherits(model, "slrm")) {
      #' predict.slrm has different syntax, 
      #' and does not support leave-one-out prediction
      lambda <- predict(model)[X]
    } else {
      lambda <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)

    }
  } else stop(paste(sQuote("lambda"),
                    "should be a vector, a pixel image,",
                    "a fitted model, or a function"))
  if(check) validate.weights(lambda)
  return(list(lambda=lambda,
              danger=danger,
              dangerous=if(danger) dangerous else NULL))
}

resolve.reciplambda <- function(X, lambda=NULL, reciplambda=NULL, ...) {
  UseMethod("resolve.reciplambda")
}

resolve.reciplambda.ppp <- function(X, lambda=NULL, reciplambda=NULL,
                                    ...,
                                    sigma=NULL, varcov=NULL,
                                    leaveoneout=TRUE, update=TRUE,
                                    check=TRUE) {
  dangerous <- c("lambda", "reciplambda")
  danger <- TRUE
  if(npoints(X) == 0) {
    danger <- FALSE
    lambda <- reciplambda <- numeric(0)
  } else if(is.null(lambda) && is.null(reciplambda)) {
    ## No intensity data provided
    danger <- FALSE
    ## Estimate density by leave-one-out kernel smoothing
    lambda <- density(X, ..., sigma=sigma, varcov=varcov,
                      at="points", leaveoneout=leaveoneout)
    lambda <- as.numeric(lambda)
    if(check) validate.weights(lambda, how="density estimation")
    reciplambda <- 1/lambda
  } else if(!is.null(reciplambda)) {
    ## 1/lambda values provided
    lambda <- NULL
    if(is.im(reciplambda)) {
      reciplambda <- safelookup(reciplambda, X)
      if(check) validate.weights(reciplambda, recip=TRUE, how="image lookup")
    } else if(is.function(reciplambda)) {
      reciplambda <- reciplambda(X$x, X$y)
      if(check) validate.weights(reciplambda, recip=TRUE, how="function evaluation")
    } else if(is.numeric(reciplambda) && is.vector(as.numeric(reciplambda))) {
      check.nvector(reciplambda, npoints(X), vname="reciplambda")
      if(check) validate.weights(reciplambda, recip=TRUE)
    } else stop(paste(sQuote("reciplambda"),
                    "should be a vector, a pixel image, or a function"))
  } else {
    #' lambda values provided
    if(is.im(lambda)) {
      lambda <- safelookup(lambda, X)
      if(check) validate.weights(lambda, how="image lookup")
    } else if(is.function(lambda)) {
      lambda <- lambda(X$x, X$y)
      if(check) validate.weights(lambda, how="function evaluation")
    } else if(is.numeric(lambda) && is.vector(as.numeric(lambda))) {
      check.nvector(lambda, npoints(X), vname="lambda")
      if(check) validate.weights(lambda)
    } else if(inherits(lambda, c("ppm", "kppm", "dppm", "slrm"))) {
      ## model provides intensity
      model <- lambda
      if(update) {
        model <- updateData(model, X)
        danger <- FALSE
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        lambda <- predict(model)[X]
      } else {
        lambda <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
      }
      if(check) validate.weights(lambda, how="model prediction")
    } else stop(paste(sQuote("lambda"),
                      "should be a vector, a pixel image,",
                      "a fitted model, or a function"))
    ## evaluate reciprocal
    reciplambda <- 1/lambda
  }
  return(list(lambda=lambda, reciplambda=reciplambda,
              danger=danger,
              dangerous=if(danger) dangerous else NULL))
}

resolve.lambdacross <- function(X, I, J,
                                lambdaI=NULL, lambdaJ=NULL, ...) {
  UseMethod("resolve.lambdacross")
}

resolve.lambdacross.ppp <- function(X, I, J,
                                    lambdaI=NULL, lambdaJ=NULL,
                                    ...,
                                    lambdaX=NULL,
                                    sigma=NULL, varcov=NULL,
                                    leaveoneout=TRUE, update=TRUE,
                                    lambdaIJ=NULL,
                                    Iexplain="points satisfying condition I",
                                    Jexplain="points satisfying condition J") {
  dangerous <- c("lambdaI", "lambdaJ")
  dangerI <- dangerJ <- TRUE
  XI <- X[I]
  XJ <- X[J]
  nI <- npoints(XI)
  nJ <- npoints(XJ)

  lamIname <- short.deparse(substitute(lambdaI))
  lamJname <- short.deparse(substitute(lambdaJ))
  bothnames <- c(lamIname, lamJname)
  givenI <- !is.null(lambdaI)
  givenJ <- !is.null(lambdaJ)
  givenX <- !is.null(lambdaX)

  if(givenI != givenJ) {
    givenone <- bothnames[c(givenI, givenJ)]
    missedone <- setdiff(bothnames, givenone)
    stop(paste("If", givenone, "is given, then",
               missedone, "should also be given"),
         call.=FALSE)
  }
  if(givenX && givenI && givenJ)
    warning(paste(paste(sQuote(bothnames), collapse=" and "),
                  "were ignored because", sQuote("lambdaX"),
                  "was given"),
            call.=FALSE)

  if(givenX) {
    ## Intensity values for all points of X
    if(is.im(lambdaX)) {
      ## Look up intensity values
      lambdaI <- safelookup(lambdaX, XI)
      lambdaJ <- safelookup(lambdaX, XJ)
    } else if(is.imlist(lambdaX) &&
              is.multitype(X) &&
              length(lambdaX) == length(levels(marks(X)))) {
      ## Look up intensity values
      Y <- split(X)
      lamY <- mapply("[", x=lambdaX, i=Y, SIMPLIFY=FALSE)
      lamX <- unsplit(lamY, marks(X))
      lambdaI <- lamX[I]
      lambdaJ <- lamX[J]
    } else if(is.function(lambdaX)) {
      ## evaluate function at locations
      if(!is.marked(X) || length(formals(lambdaX)) == 2) {
        lambdaI <- lambdaX(XI$x, XI$y)
        lambdaJ <- lambdaX(XJ$x, XJ$y)
      } else {
        lambdaI <- lambdaX(XI$x, XI$y, marks(XI))
        lambdaJ <- lambdaX(XJ$x, XJ$y, marks(XJ))
      }
    } else if(is.numeric(lambdaX) && is.vector(as.numeric(lambdaX))) {
      ## vector of intensity values
      if(length(lambdaX) != npoints(X))
        stop(paste("The length of", sQuote("lambdaX"),
                   "should equal the number of points of X"))
      lambdaI <- lambdaX[I]
      lambdaJ <- lambdaX[J]
    } else if(inherits(lambdaX, c("ppm", "kppm", "dppm", "slrm"))) {
      ## point process model provides intensity
      model <- lambdaX
      if(update) {
        model <- update(model, X)
        dangerI <- dangerJ <- FALSE
        dangerous <- "lambdaIJ"
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        Lambda <- predict(model)
        lambdaI <- Lambda[XI]
        lambdaJ <- Lambda[XJ]
      } else {
        ## re-fit model to data X
        lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        lambdaI <- lambdaX[I]
        lambdaJ <- lambdaX[J]
      }
    } else stop(paste("Argument lambdaX is not understood:",
                      "it should be a numeric vector,",
                      "an image, a function(x,y)",
                      "or a fitted point process model (ppm, kppm or dppm)"))
  } else {
    ## lambdaI, lambdaJ expected
    if(!givenI) {
      ## estimate intensity
      dangerI <- FALSE
      dangerous <- setdiff(dangerous, "lambdaI")
      lambdaI <- density(XI, ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=leaveoneout)
    } else if(is.im(lambdaI)) {
      ## look up intensity values
      lambdaI <- safelookup(lambdaI, XI)
    } else if(is.function(lambdaI)) {
      ## evaluate function at locations
      lambdaI <- lambdaI(XI$x, XI$y)
    } else if(is.numeric(lambdaI) && is.vector(as.numeric(lambdaI))) {
      ## validate intensity vector
      check.nvector(lambdaI, nI, things=Iexplain, vname="lambdaI")
    } else if(inherits(lambdaI, c("ppm", "kppm", "dppm", "slrm"))) {
      ## point process model provides intensity
      model <- lambdaI
      if(update) {
        model <- updateData(model, X)
        dangerI <- FALSE
        dangerous <- setdiff(dangerous, "lambdaI")
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        lambdaI <- predict(model)[XI]
      } else {
        lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        lambdaI <- lambdaX[I]
      }
    } else stop(paste(sQuote("lambdaI"), "should be a vector or an image"))
    
    if(!givenJ) {
      ## estimate intensity
      dangerJ <- FALSE
      dangerous <- setdiff(dangerous, "lambdaJ")
      lambdaJ <- density(XJ, ..., sigma=sigma, varcov=varcov,
                         at="points", leaveoneout=leaveoneout)
    } else if(is.im(lambdaJ)) {
      ## look up intensity values
      lambdaJ <- safelookup(lambdaJ, XJ)
    } else if(is.function(lambdaJ)) {
      ## evaluate function at locations
      lambdaJ <- lambdaJ(XJ$x, XJ$y)
    } else if(is.numeric(lambdaJ) && is.vector(as.numeric(lambdaJ))) {
      ## validate intensity vector
      check.nvector(lambdaJ, nJ, things=Jexplain, vname="lambdaJ")
    } else if(inherits(lambdaJ, c("ppm", "kppm", "dppm", "slrm"))) {
      ## point process model provides intensity
      model <- lambdaJ
      if(update) {
        model <- updateData(model, X)
        dangerJ <- FALSE
        dangerous <- setdiff(dangerous, "lambdaJ")
      }
      if(inherits(model, "slrm")) {
        #' predict.slrm has different syntax, 
        #' and does not support leave-one-out prediction
        lambdaJ <- predict(model)[XJ]
      } else {
        lambdaX <- fitted(model, dataonly=TRUE, leaveoneout=leaveoneout)
        lambdaJ <- lambdaX[J]
      }
    } else 
      stop(paste(sQuote("lambdaJ"), "should be a vector or an image"))
  }
  
  ## Weight for each pair
  if(!is.null(lambdaIJ)) {
    dangerIJ <- TRUE
    dangerous <- union(dangerous, "lambdaIJ")
    if(!is.matrix(lambdaIJ))
      stop("lambdaIJ should be a matrix")
    if(nrow(lambdaIJ) != nI)
      stop(paste("nrow(lambdaIJ) should equal the number of", Iexplain))
    if(ncol(lambdaIJ) != nJ)
      stop(paste("ncol(lambdaIJ) should equal the number of", Jexplain))
  } else {
    dangerIJ <- FALSE
  }
    
  danger <- dangerI || dangerJ || dangerIJ
    
  return(list(lambdaI = lambdaI, lambdaJ = lambdaJ, lambdaIJ=lambdaIJ,
                danger = danger, dangerous = dangerous))
}


validate.weights <- function(x, recip=FALSE, how = NULL,
                             allowzero = recip,
                             allowinf  = !recip) {
  ra <- range(x)
  offence <-
    if(!allowinf && !all(is.finite(ra)))  "infinite" else
    if(ra[1] < 0)                         "negative" else
    if(!allowzero && ra[1] == 0)          "zero" else NULL
  if(!is.null(offence)) {
    xname <- deparse(substitute(x))
    offenders <- paste(offence, "values of", sQuote(xname))
    if(is.null(how))
      stop(paste(offenders, "are not allowed"), call.=FALSE)
    stop(paste(how, "yielded", offenders), call.=FALSE)
  }
  return(TRUE)
}

updateData <- function(model, X, ...) {
  UseMethod("updateData")
}

updateData.default <- function(model, X, ..., warn=TRUE) {
  ## for some bizarre reason, method dispatch often fails for this function
  ## so we do it by hand as a backup
  if(warn) warning("Reached 'updateData.default'", call.=FALSE)
  if(inherits(model, c("ppm", "kppm", "dppm", "slrm"))) {
    if(inherits(model, "ppm")) {
      model <- updateData.ppm(model, X)
    } else if(inherits(model, "kppm")) {
      model <- updateData.kppm(model, X)
    } else if(inherits(model, "dppm")) {
      model <- updateData.dppm(model, X)
    } else if(inherits(model, "slrm")) {
      model <- updateData.slrm(model, X)
    }
  } else if(inherits(model, "lppm")) {
    if(requireNamespace("spatstat.linnet")) {
      model <- spatstat.linnet::updateData.lppm(model, X)
    } else
      if(warn)
        warning("Model was not updated; this requires a recent version of spatstat.linnet", call.=FALSE)
  } else
    if(warn)
      warning("Unrecognised kind of 'model'; no update performed",
              call.=FALSE)
  return(model)
}
  
