##
##
##  rocmodel.R
##
##  Calculate ROC curve
##       (in spatstat.model)
##
## Copyright (c) 2017-2025 Adrian Baddeley/Ege Rubak/Rolf Turner
##
##





roc.ppm <- roc.kppm <- roc.slrm <-
  function(X, covariate=NULL, ...,
           baseline=NULL, high=TRUE,
           method = "raw", CI = "none", alpha=0.05,
           leaveoneout=FALSE, subset=NULL) {
  model <- X
  X <- response(model)
  traditional <- is.null(covariate)
  lambda <- do.call.matched(predict,
                            list(object=model, ...),
                            c("object", "ngrid", "dimyx", "eps",
                              "correction", "new.coef"))
  if(traditional) {
    #' discriminant is the predicted intensity
    covariate <- lambda
    covtype <- if(is.slrm(model)) "probability" else "intensity"
  } else {
    #' discriminant is a user-specified covariate
    covariate <- digestCovariates(covariate, W=Window(model))
    if(length(covariate) == 1) covariate <- covariate[[1L]]
    covtype <- "covariate"
  }
  if(is.ppp(baseline) || is.lpp(baseline)) {
    if(is.slrm(model) && is.null(model$Data$dataAtPoints)) {
      W <- model$Data$W
      X <- discretise(X, xy=W, move.points=TRUE)
      baseline <- discretise(baseline, xy=W, move.points=TRUE)
    }
    result <- rocDummy(X, baseline, covariate, ..., high=high, method=method,
                       CI=CI, alpha=alpha, subset=subset)
  } else {
    nullmodel <- resolveNullModel(baseline, model)
    leaveoneout <- if(!traditional) FALSE else as.logical(leaveoneout)
    ## leaveoneout can be TRUE, FALSE or c(TRUE, FALSE)
    if(any(leaveoneout)) {
      fittedX <- fitted(model, dataonly=TRUE, leaveoneout=TRUE)
      Rminus <- rocEngine(covariate, nullmodel, method = method,
                          covtype = covtype,
                          fittedmodel = if(traditional) NULL else model,
                          discrimAtPoints = fittedX, # NB
                          leftoneout=TRUE, # for labelling
                          high=high,
                          CI = CI, alpha=alpha,
                          subset=subset,
                          ...)
    } else {
      Rminus <- NULL
    }
    if(any(!leaveoneout)) {
      Rplus <- rocEngine(covariate, nullmodel, method = method,
                         covtype = covtype,
                         fittedmodel = if(traditional) NULL else model,
                         discrimAtPoints = NULL,  # NB
                         high=high,
                         CI = CI, alpha=alpha,
                         subset=subset,
                         ...)
    } else {
      Rplus <- NULL
    }
    if(all(leaveoneout)) {
      ## return only the leave-one-out calculation
      result <- Rminus
    } else if(all(!leaveoneout)) {
      ## return only the vanilla calculation
      result <- Rplus
    } else {
      ## return both results, combined
      common <- c("null", "theo", "fun", "thresh")
      if(leaveoneout[1]) {
        Rfirst <- Rminus
        Rsecond <- Rplus
      } else {
        Rfirst <- Rplus
        Rsecond <- Rminus
      }
      result <- bind.fv(Rfirst,
                        Rsecond[, setdiff(colnames(Rsecond), common)])
      ## preferred column
      fvnames(result, ".y") <- fvnames(Rfirst, ".y")
      ## shade columns if any
      if(!is.null(sn <- fvnames(Rfirst, ".s")))
        fvnames(result, ".s") <- sn
      ## dotnames
      dn <- fvnames(Rplus, ".")
      dnew <- character(0)
      for(d in dn) {
        if(!leaveoneout[1]) dnew <- c(dnew, d)
        if(!(d %in% common))
          dnew <- c(dnew, makeRocTag(d, TRUE))
        if(leaveoneout[1]) dnew <- c(dnew, d)
      }
      fvnames(result, ".") <- dnew
      ## return as roc
      class(result) <- union("roc", class(result))
    }
  }
  return(result)
}









