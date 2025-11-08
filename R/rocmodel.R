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
## Source file roc.R      $Revision$ $Date$
## ..................................................................






## ------------ roc method for ppm ------------------
#' [rocmodelcode.R]
#' Code template included multiple times in roc.R
#' $Revision: 1.2 $ $Date: 2025/11/08 03:08:06 $

roc.ppm <- function(X, covariate=NULL, ...,
           baseline=NULL, high=TRUE,
           method = "raw", CI = "none", alpha=0.05,
           leaveoneout=FALSE, subset=NULL) {
  model <- X
  X <- response(model)
  leaveoneout <- unique(as.logical(leaveoneout))
  meantype <- "trend"
  if(traditional <- is.null(covariate)) {
    #' discriminant is the fitted model intensity or presence probability
    covtype <- "intensity"
    covariate <- do.call.matched(predict,
                                 list(object=model, type=meantype, ...),
                                 c("object", "type",
                                   "ngrid", "dimyx", "eps",
                                   "correction", "new.coef"))
  } else {
    #' discriminant is a user-specified covariate
    covtype <- "covariate"
    covariate <- digestCovariates(covariate, W=Window(model))
    if(length(covariate) == 1) covariate <- covariate[[1L]]
    if(any(leaveoneout)) {
      warning("Argument leaveoneout=TRUE ignored because covariate is provided",
              call.=FALSE)
      leaveoneout <- FALSE
    }
  }
  if(is.ppp(baseline) || is.lpp(baseline)) {
    result <- rocDummy(X, baseline, covariate, ..., high=high, method=method,
                       CI=CI, alpha=alpha, subset=subset)
  } else {
    nullmodel <- resolveNullModel(baseline, model)
    ## leaveoneout can be TRUE, FALSE or c(TRUE, FALSE) or c(FALSE,TRUE)
    if(any(!leaveoneout)) {
      ## calculate vanilla estimate
      Rplus <- rocEngine(covariate, nullmodel, method = method,
                         covtype = covtype,
                         fittedmodel = if(traditional) NULL else model,
                         high=high,
                         CI = CI, alpha=alpha,
                         subset=subset,
                         lambdatype=meantype, # for prediction from null model
                         ...)
    } else {
      Rplus <- NULL
    }
    if(any(leaveoneout)) {
      ## calculate leave-one-out estimate
      ## use leave-one-POINT-out intensity at DATA POINTS only
      fittedX <- fitted(model, dataonly=TRUE, leaveoneout=TRUE, type=meantype)
      Rminus <- rocEngine(covariate, nullmodel, method = method,
                          covtype = covtype,
                          fittedmodel = if(traditional) NULL else model,
                          discrimAtPoints = fittedX, # <<<< NB >>>>>
                          leftoneout=TRUE, # for labelling
                          high=high,
                          CI = CI, alpha=alpha,
                          subset=subset,
                          lambdatype=meantype, # for prediction from null
                          ...)
    } else {
      Rminus <- NULL
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
  
## ------------ roc method for kppm ------------------
#' [rocmodelcode.R]
#' Code template included multiple times in roc.R
#' $Revision: 1.2 $ $Date: 2025/11/08 03:08:06 $

roc.kppm <- function(X, covariate=NULL, ...,
           baseline=NULL, high=TRUE,
           method = "raw", CI = "none", alpha=0.05,
           leaveoneout=FALSE, subset=NULL) {
  model <- X
  X <- response(model)
  leaveoneout <- unique(as.logical(leaveoneout))
  meantype <- "trend"
  if(traditional <- is.null(covariate)) {
    #' discriminant is the fitted model intensity or presence probability
    covtype <- "intensity"
    covariate <- do.call.matched(predict,
                                 list(object=model, type=meantype, ...),
                                 c("object", "type",
                                   "ngrid", "dimyx", "eps",
                                   "correction", "new.coef"))
  } else {
    #' discriminant is a user-specified covariate
    covtype <- "covariate"
    covariate <- digestCovariates(covariate, W=Window(model))
    if(length(covariate) == 1) covariate <- covariate[[1L]]
    if(any(leaveoneout)) {
      warning("Argument leaveoneout=TRUE ignored because covariate is provided",
              call.=FALSE)
      leaveoneout <- FALSE
    }
  }
  if(is.ppp(baseline) || is.lpp(baseline)) {
    result <- rocDummy(X, baseline, covariate, ..., high=high, method=method,
                       CI=CI, alpha=alpha, subset=subset)
  } else {
    nullmodel <- resolveNullModel(baseline, model)
    ## leaveoneout can be TRUE, FALSE or c(TRUE, FALSE) or c(FALSE,TRUE)
    if(any(!leaveoneout)) {
      ## calculate vanilla estimate
      Rplus <- rocEngine(covariate, nullmodel, method = method,
                         covtype = covtype,
                         fittedmodel = if(traditional) NULL else model,
                         high=high,
                         CI = CI, alpha=alpha,
                         subset=subset,
                         lambdatype=meantype, # for prediction from null model
                         ...)
    } else {
      Rplus <- NULL
    }
    if(any(leaveoneout)) {
      ## calculate leave-one-out estimate
      ## use leave-one-POINT-out intensity at DATA POINTS only
      fittedX <- fitted(model, dataonly=TRUE, leaveoneout=TRUE, type=meantype)
      Rminus <- rocEngine(covariate, nullmodel, method = method,
                          covtype = covtype,
                          fittedmodel = if(traditional) NULL else model,
                          discrimAtPoints = fittedX, # <<<< NB >>>>>
                          leftoneout=TRUE, # for labelling
                          high=high,
                          CI = CI, alpha=alpha,
                          subset=subset,
                          lambdatype=meantype, # for prediction from null
                          ...)
    } else {
      Rminus <- NULL
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
  
## ------------ roc method for slrm ------------------
#' [rocmodelcode.R]
#' Code template included multiple times in roc.R
#' $Revision: 1.2 $ $Date: 2025/11/08 03:08:06 $

roc.slrm <- function(X, covariate=NULL, ...,
           baseline=NULL, high=TRUE,
           method = "raw", CI = "none", alpha=0.05,
           leaveoneout=FALSE, subset=NULL) {
  model <- X
  X <- response(model)
  leaveoneout <- unique(as.logical(leaveoneout))
  meantype <- "probabilities" 
  if(traditional <- is.null(covariate)) {
    #' discriminant is the fitted model intensity or presence probability
    covtype <- "probability" 
    covariate <- do.call.matched(predict,
                                 list(object=model, type=meantype, ...),
                                 c("object", "type",
                                   "ngrid", "dimyx", "eps",
                                   "correction", "new.coef"))
  } else {
    #' discriminant is a user-specified covariate
    covtype <- "covariate"
    covariate <- digestCovariates(covariate, W=Window(model))
    if(length(covariate) == 1) covariate <- covariate[[1L]]
    if(any(leaveoneout)) {
      warning("Argument leaveoneout=TRUE ignored because covariate is provided",
              call.=FALSE)
      leaveoneout <- FALSE
    }
  }
  if(is.ppp(baseline) || is.lpp(baseline)) {
    if(is.null(model$Data$dataAtPoints)) {
      W <- model$Data$W
      X <- discretise(X, xy=W, move.points=TRUE)
      baseline <- discretise(baseline, xy=W, move.points=TRUE)
    }
    result <- rocDummy(X, baseline, covariate, ..., high=high, method=method,
                       CI=CI, alpha=alpha, subset=subset)
  } else {
    nullmodel <- resolveNullModel(baseline, model)
    ## leaveoneout can be TRUE, FALSE or c(TRUE, FALSE) or c(FALSE,TRUE)
    if(any(!leaveoneout)) {
      ## calculate vanilla estimate
      Rplus <- rocEngine(covariate, nullmodel, method = method,
                         covtype = covtype,
                         fittedmodel = if(traditional) NULL else model,
                         high=high,
                         CI = CI, alpha=alpha,
                         subset=subset,
                         lambdatype=meantype, # for prediction from null model
                         ...)
    } else {
      Rplus <- NULL
    }
    if(any(leaveoneout)) {
      ## calculate leave-one-out estimate
      ## covariate = leave-one-PIXEL-out probability at every pixel
      loocov <- do.call.matched(predict,
                                list(object=model,
                                     type=meantype,
                                     leaveoneout=TRUE,
                                     ...),
                                c("object", "type",
                                  "ngrid", "dimyx", "eps",
                                  "correction", "new.coef",
                                  "leaveoneout"))
      Rminus <- rocEngine(loocov, nullmodel, method = method,
                          covtype = covtype,
                          fittedmodel = if(traditional) NULL else model,
                          leftoneout=TRUE, # for labelling
                          high=high,
                          CI = CI, alpha=alpha,
                          subset=subset,
                          lambdatype=meantype, # for prediction from null
                          ...)
    } else {
      Rminus <- NULL
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
  










