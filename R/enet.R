#'
#'  enet.R
#'
#'  Fit ppm using elastic net
#'
#'  Original by Achmad Choiruddin
#'  plus original by Suman Rakshit
#'  modified by Adrian Baddeley
#'
#'  Copyright (c) Achmad Choiruddin, Suman Rakshit and Adrian Baddeley 2022
#'  GNU Public Licence >= 2.0
#' 
#'  $Revision: 1.14 $ $Date: 2022/06/20 06:23:10 $


enet.engine <-function(model, ...,
                       standardize=TRUE, lambda=NULL, alpha=1, adaptive=FALSE) {
  #' alpha=1, enet
  #' alpha=0, ridge
  #' 0<alpha<1,enet
  verifyclass(model, "ppm")
  kraever("glmnet")
  X <- response(model)
  nX <- npoints(X)
  #' which composite likelihood
  CL <- switch(model$method,
               mpl = "mpl",
               ho = "mpl",
               logi = "logi",
               VBlogi = "logi",
               stop(paste("Unrecognised fitting method", sQuote(model$method))))
  #' ensure model was actually fitted by glm or gam
  switch(model$fitter,
         exact = {
           model <- update(model, forcefit=TRUE)
         },
         gam = {
           warning("use.gam=TRUE was ignored because we are using glmnet")
           model <- update(model, use.gam=FALSE, forcefit=TRUE)
         },
         {})
  ## extract Berman-Turner response and weights
  gd <- getglmdata(model)
  switch(CL,
         mpl = {
           yy  <- gd$.mpl.Y
           wts <- gd$.mpl.W
           sub <- gd$.mpl.ok
         },
         logi = {
           yy  <- gd$.logi.Y
           wts <- gd$.logi.w
           sub  <- gd$.logi.ok
         })
  ## extract canonical variables
  mm <- model.matrix(model)
  if(nrow(mm) != nrow(gd))
    stop("Internal error: mismatch in glmdata and model.matrix")
  ## get offset 
  gm <- getglmfit(model)
  off <- model.offset(gm) %orifnull% numeric(nrow(gd))
  ## restrict to subset
  if(!all(sub)) {
    yy <- yy[sub]
    wts <- wts[sub]
    mm <- mm[sub, , drop=FALSE]
    gd <- gd[sub, , drop=FALSE]
    off <- off[sub]
  }
  ## initial estimates of coefficients
  beta.ini <- coef(model)
  intercept.position <- match("(Intercept)", names(beta.ini)) # can be NA
  has.intercept <- !is.na(intercept.position)
  beta.ini.strip <- if(has.intercept) beta.ini[-intercept.position] else beta.ini
  mm.strip <- if(has.intercept) mm[,-intercept.position, drop=FALSE] else mm
  ## penalty factor
  penalty.factor <- if(adaptive) 1/abs(beta.ini.strip) else rep.int(1, ncol(mm.strip))
  penalty.factor <- pmin(penalty.factor, 1/.Machine$double.eps)
  ## threshold
  thresh <- if(is.null(lambda)) 1e-7 else 1e-30
  ## go
  gfit <- glmnet::glmnet(mm.strip, yy,
                         weights=wts, offset=off,
                         family= switch(CL,
                                        mpl=quasipoisson(link = "log"),
                                        logi=binomial()),
                         standardize=standardize,
                         alpha=alpha, lambda=lambda,
                         penalty.factor=penalty.factor,
                         thresh=thresh, intercept=has.intercept)
  lambda <- as.matrix(gfit$lambda)
  beta.estimates <- as.matrix(gfit$beta)
  if(has.intercept) {
    ## reinstate intercept in original position
    ii <- 1:nrow(beta.estimates)
    beta.estimates <- rbind(beta.estimates[ii < intercept.position, , drop=FALSE],
                            "(Intercept)"=gfit$a0,
                            beta.estimates[ii >= intercept.position, , drop=FALSE])
  }
  ## Tuning parameter selection
  stuff <- apply(beta.estimates, 2,
                 function(beta) {
                   eta <- mm %*% beta + off
                   v <- switch(CL, mpl = exp(eta), logi = log(1+exp(eta)))
                   loglike <- sum(wts * (yy * eta - v)) 
                   bic <- -2*loglike + (sum(beta != 0) - has.intercept) * log(nX)
                   return(c(loglike=loglike, bic=bic))
                 })
  jopt <- which.min(stuff["bic", ])
  optim.bic <- stuff["bic", jopt]
  optim.loglike <- stuff["loglike", jopt]
  optim.beta <- beta.estimates[,jopt]
  optim.lambda <- lambda[jopt]
  ## update fitted model object
  newmodel <- model
  newmodel$coef <- optim.beta
  newmodel$maxlogpl <- optim.loglike
  newmodel$improve.type <- "enet"
  ## recompute/update internals
  newmodel$fisher <- NULL
  newmodel$varcov <- NULL
  newmodel$fitin <- NULL
  newmodel$fitin <- fitin(newmodel)
  ## save glmnet information
  newmodel$internal$glmnet <- list(fit=gfit,
                                   lambda=lambda, alpha=alpha,
                                   standardize=standardize, adaptive=adaptive,
                                   intercept=has.intercept,
                                   estimates=beta.estimates,
                                   criteria=stuff)
  ## remember original estimates
  newmodel$coef.orig <- beta.ini
  newmodel$maxlogpl.orig <- model$maxlogpl
  return(newmodel)
}

