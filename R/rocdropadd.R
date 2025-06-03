#'
#'    rocdropadd.R
#'
#'  ROC for adding single terms to a point process model
#'          /dropping single terms from a point process model
#'
#' Copyright (c) 2017-2025 Adrian Baddeley, Ege Rubak and Rolf Turner


addROC <- function(object, scope, high=TRUE, ...) {
  addapply(object, scope=scope, action="roc", high=high, ...)
}

addapply <- function(object, 
                     action=c("berman.test", "cdf.test", "rhohat", "roc"),
                     scope, ..., high=TRUE) {
  action <- match.arg(action)
  modelclasses <- c("ppm", "kppm", "dppm", "slrm", "lppm")
  om <- inherits(object, modelclasses)
  sm <- inherits(scope, modelclasses)
  if(!(om || sm))
    stop(paste("Either 'object' or 'scope' should be a fitted model of class",
               commasep(modelclasses, "or")),
         call.=FALSE)
  ## determine reference model and variables (refvars, maxvars)
  if(om) {
    ## 'object' is a fitted model defining the reference model
    reference <- object
    X <- response(reference)
    refvars <- termsinformula(formula(reference))
    #' determine maximal variables from 'scope' 
    if(length(scope) == 0) 
      stop("Argument `scope` contains no data") 
    if(sm) {
      ## 'scope' is a model too
      maxvars <- attr(terms(scope), "term.labels")
    } else if(inherits(scope, "formula")) {
      maxvars <- termsinformula(update.formula(reference, scope))
    } else if(is.character(scope)) {
      maxvars <- scope
    } else stop("Format of argument `scope` is not understood")
  } else {
    ## 'scope' is a fitted model defining the maximal model
    maxmodel <- scope
    X <- response(maxmodel)
    maxvars <- termsinformula(formula(maxmodel))
    #' determine the reference model from 'object'
    if(length(object) == 0) 
      stop("Argument `object` contains no data")
    if(inherits(object, "formula")) {
      reference <- update(maxmodel, object)
    } else stop("Format of argument `object` is not understood")
    refvars <- termsinformula(formula(reference))
  }
  ## finally determine scope
  scope <- setdiff(maxvars, refvars)
  ## handle 'high'
  if(length(high) == 1) high <- rep(high, length(scope)) else
  stopifnot(length(high) == length(scope))
  #'
  #' pixel resolution arguments only
  rezargs <- list(...)[c("eps", "dimyx", "xy", "rule.eps")]
  rezargs <- rezargs[!sapply(rezargs, is.null)]
  ## go
  ns <- length(scope)
  result <- vector(mode="list", length=ns)
  for(i in seq_len(ns)) {
    tt <- scope[i]
    add.i <- as.formula(paste(". ~ . + ", tt))
    fit.i <- update(reference, add.i, env=environment(formula(reference)))
    mim.i <- do.call(model.images,
                     append(list(object=fit.i),
                            rezargs))
    cov.i <- mim.i[[tt]]
    result[[i]] <-
      switch(action,
             roc = {
               roc(X, cov.i, ..., baseline=reference, high=high[i])
             },
             rhohat = {
               rhohat(reference, cov.i, ..., covname=tt)
             },
             berman.test = {
               berman.test(reference, cov.i, ..., covname=tt)
             },
             cdf.test = {
               cdf.test(reference, cov.i, ..., covname=tt)
             })
  }
  names(result) <- scope
  return(as.anylist(result))
}
  
dropROC <- function(object, scope=NULL, high=TRUE, ...) {
  dropply(object, action="roc", scope=scope, high=high, ...)
}

dropply <- function(object,
                    action=c("berman.test", "cdf.test", "rhohat", "roc"),
                    scope=NULL, ..., high=TRUE) {
  modelclasses <- c("ppm", "kppm", "dppm", "slrm", "lppm")
  if(!inherits(object, modelclasses))
    stop(paste("'object' should be a fitted model of class",
               commasep(modelclasses, "or")),
         call.=FALSE)
  action <- match.arg(action)
  X <- response(object)
  tl <- termsinformula(formula(object))
  if(is.null(scope)) {
    scope <- drop.scope(object)
  } else {
    if (!is.character(scope)) 
      scope <- termsinformula(update.formula(object, scope))
    if (!all(match(scope, tl, 0L) > 0L)) 
      stop("scope is not a subset of term labels")
  }
  if(length(high) == 1) high <- rep(high, length(scope)) else
  stopifnot(length(high) == length(scope))
  #' pixel resolution arguments only
  rezargs <- list(...)[c("eps", "dimyx", "xy", "rule.eps")]
  rezargs <- rezargs[!sapply(rezargs, is.null)]
  #' recognised covariates
  covars <- do.call(model.images,
                    append(list(object=object),
                           rezargs))
  W <- if(length(covars)) Window(covars[[1]]) else Window(X)
  wanted <- setdiff(scope, names(covars))
  if("x" %in% wanted) {
    covars$x <- as.im(function(x,y){x}, W)
    wanted <- setdiff(wanted, "x")
  }
  if("y" %in% wanted) {
    covars$y <- as.im(function(x,y){y}, W)
    wanted <- setdiff(wanted, "y")
  }
  if((nw <- length(wanted)) > 0) {
    #' find correspondence between original and canonical covariates
    dmat <- model.depends(object)
    orig <- colnames(dmat)
    canon <- rownames(dmat)
    mapisunique <- (colSums(dmat) == 1)
    kk <- match(wanted, orig)
    matched <- !is.na(kk)
    for(k in kk[matched]) {
      if(mapisunique[k]) {
        origname <- orig[k]
        canoname <- canon[which(dmat[,k])]
        covars[[origname]] <- covars[[canoname]]
      } 
    }
    #' deal with interactions etc
    if(any(hard <- !mapisunique[kk[matched]])) {
      nhard <- sum(hard)
      hardnames <- (wanted[matched])[hard]
      warning(paste(ngettext(nhard, "Variable", "Variables"),
                    commasep(sQuote(hardnames)),
                    ngettext(nhard, "is", "are"),
                    "difficult to remove from the model"),
            call.=FALSE)
    }
    wanted <- setdiff(wanted, names(covars))
    nw <- length(wanted)
    if(nw > 0) {
      warning(paste(ngettext(nw, "Variable", "Variables"),
                    commasep(sQuote(wanted)),
                    ngettext(nw, "is", "are"),
                    "not available in the model"),
              call.=FALSE)
      #' give up and remove these variables from the scope.
      scope <- setdiff(scope, wanted)
    }
  }
  #' 
  ns <- length(scope)
  result <- vector(mode="list", length=ns)
  for(i in seq_len(ns)) {
    tt <- scope[i]
    drop.i <- as.formula(paste(". ~ . - ", tt))
    fit.i <- update(object, drop.i, env=environment(formula(object)))
    cov.i <- covars[[tt]]
    result[[i]] <-
      switch(action,
             roc = {
               roc(X, cov.i, ..., baseline=fit.i, high=high[i])
             },
             rhohat = {
               rhohat(fit.i, cov.i, ..., covname=tt)
             },
             berman.test = {
               berman.test(fit.i, cov.i, ..., covname=tt)
             },
             cdf.test = {
               cdf.test(fit.i, cov.i, ..., covname=tt)
             })
  }
  names(result) <- scope
  return(as.anylist(result))
}
  
