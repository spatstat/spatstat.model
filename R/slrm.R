#
#  slrm.R
#
#  Spatial Logistic Regression
#
#  $Revision: 1.63 $   $Date: 2022/08/11 03:53:48 $
#

slrm <- function(formula, ..., data=NULL, offset=TRUE, link="logit",
                 dataAtPoints=NULL, splitby=NULL) {
  
  # remember call
  CallInfo <- list(callstring = short.deparse(sys.call()),
                   cl = match.call(),
                   formula = formula,
                   offset=offset,
                   link=link,
                   splitby=splitby,
                   dotargs=list(...))
  if(!(link %in% c("logit", "cloglog")))
    stop(paste("Unrecognised link", dQuote(link)))

  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(formula, "formula"))
    stop(paste("Argument", dQuote("formula"), "should be a formula"))

  # check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Argument", sQuote("formula"),
               "must have a left hand side"))
  Yname <- formula[[2]]
  trend <- rhs <- formula[c(1,3)]
  if(!is.name(Yname))
    stop("Left hand side of formula should be a single name")
  Yname <- paste(Yname)
  if(!inherits(trend, "formula"))
    stop("Internal error: failed to extract RHS of formula")

  varnames <- unique(variablesinformula(trend))
  if(isTRUE(CallInfo$dotargs$save.all.vars))
    varnames <- union(varnames, names(data))
  specials <- c("x", "y", "logpixelarea")
  covnames <- varnames[!(varnames %in% specials)]

  # add 'splitby' to covariate names
  if(!is.null(splitby)) {
    if(!is.character(splitby) || length(splitby) != 1)
      stop("splitby should be a single character string")
    covnames <- unique(c(covnames, splitby))
  }

  CallInfo$responsename <- Yname
  CallInfo$varnames     <- varnames
  CallInfo$covnames     <- covnames
  
  # Parent environment
  parenv <- environment(formula)

  ########  FIND DATA AND RESHAPE #######################

  Data <- slr.prepare(CallInfo, parenv, data, dataAtPoints, splitby)

#  W  <- Data$W
  df <- Data$df
  nY <- npoints(Data$response)
  
  ########  FIT  ###############################

  dformula <- formula
  if(offset) {
    # insert offset term in formula
    rhs <- paste(as.character(rhs), collapse=" ")
    rhs <- paste(c(rhs, "offset(logpixelarea)"), collapse="+")
    dformula <- as.formula(paste(Yname, rhs))
  }

  linkname <- link
  FIT  <- glm(dformula, family=binomial(link=linkname),
              data=df, na.action=na.exclude)

  result <- list(call     = CallInfo$cl,
                 CallInfo = CallInfo,
                 Data     = Data,
                 Fit      = list(FIT=FIT, dformula=dformula),
                 terms    = terms(formula),
                 nobs     = nY)

  class(result) <- c("slrm", class(result))
  return(result)
}

################ UTILITY TO FIND AND RESHAPE DATA #################

slr.prepare <- local({

  slr.prepare <- function(CallInfo, envir, data,
                        dataAtPoints=NULL, splitby=NULL,
                        clip=TRUE) {
    ## CallInfo is produced by slrm()
    ## envir is parent environment of model formula
    ## data  is 'data' argument that takes precedence over 'envir'
    ## 'clip' is TRUE if the data should be clipped to the domain of Y
    Yname    <- CallInfo$responsename
    ##  varnames <- CallInfo$varnames
    covnames <- CallInfo$covnames
    dotargs  <- CallInfo$dotargs
    ##
    ## Get the response point pattern Y 
    Y <- getobj(Yname, envir, data)
    if(!is.ppp(Y))
      stop(paste("The response", sQuote(Yname), "must be a point pattern"))
    ##
    if(!is.null(dataAtPoints)) {
      dataAtPoints <- as.data.frame(dataAtPoints)
      if(nrow(dataAtPoints) != npoints(Y))
        stop(paste("dataAtPoints should have one row for each point in",
                   dQuote(Yname)))
    }
    ## Find the covariates
    ncov <- length(covnames)
    covlist <- lapply(as.list(covnames), getobj, env = envir, dat=data)
    names(covlist) <- covnames
    ## Each covariate should be an image, a window, a function, or single number
    if(ncov == 0) {
      isim <- isowin <- ismask <- isfun <- isnum <-
        isspatial <- israster <- logical(0)
    } else {
      isim  <- sapply(covlist, is.im)
      isowin  <- sapply(covlist, is.owin)
      ismask  <- sapply(covlist, is.mask)
      isfun  <- sapply(covlist, is.function)
      isspatial <- isim | isowin | isfun
      israster <- isim | ismask
      isnum <- sapply(covlist, is.numeric) & (lengths(covlist) == 1)
    }
    if(!all(ok <- (isspatial | isnum))) {
      n <- sum(!ok)
      stop(paste(ngettext(n, "The argument", "Each of the arguments"),
                 commasep(sQuote(covnames[!ok])),
                 "should be either an image, a window, or a single number"))
    }
    ## 'splitby' 
    if(!is.null(splitby)) {
      splitwin <- covlist[[splitby]]
      if(!is.owin(splitwin))
        stop("The splitting covariate must be a window")
      ## ensure it is a polygonal window
      covlist[[splitby]] <- splitwin <- as.polygonal(splitwin)
      ## delete splitting covariate from lists to be processed
      issplit <- (covnames == splitby)
      isspatial[issplit] <- FALSE
      israster[issplit] <- FALSE
    }
    ## 
    ##  nnum <- sum(isnum)
    ##  nspatial <- sum(isspatial)
    nraster <- sum(israster)
    ##
    numlist <- covlist[isnum]
    spatiallist <- covlist[isspatial]
    rasterlist <- covlist[israster]
    ##
    numnames <- names(numlist)
    spatialnames <- names(spatiallist)
    ##  rasternames <- names(rasterlist)
    ##
  
    ########  CONVERT TO RASTER DATA  ###############################

    ## determine spatial domain & common resolution: convert all data to it
    if(length(dotargs) > 0 || nraster == 0) {
      ## Pixel resolution is determined by explicit arguments
      if(clip) {
        ## Window extent is determined by response point pattern
        D <- as.owin(Y)
      } else {
        ## Window extent is union of domains of data
        domains <- lapply(append(spatiallist, list(Y)), as.owin)
        D <- do.call(union.owin, domains)
      }
      ## Create template mask
      W <- do.call.matched(as.mask, append(list(w=D), dotargs))
      ## Convert all spatial objects to this resolution
      spatiallist <- lapply(spatiallist, convert, W=W)
    } else {
      ## Pixel resolution is determined implicitly by covariate data
      W <- do.call(commonGrid, rasterlist)
      if(clip) {
        ## Restrict data to spatial extent of response point pattern
        W <- intersect.owin(W, as.owin(Y))
      }
      ## Adjust spatial objects to this resolution
      spatiallist <- lapply(spatiallist, convert, W=W)
    }
    ## images containing coordinate values
    xcoordim <- as.im(function(x,y){x}, W=W)
    ycoordim <- as.im(function(x,y){y}, W=W)
    ##
    ## create a list of covariate images, with names as in formula
    covimages <- append(list(x=xcoordim, y=ycoordim), spatiallist)

    basepixelarea <- W$xstep * W$ystep

    ########  ASSEMBLE DATA FRAME  ###############################

    if(is.null(splitby)) {
      df <- slrAssemblePixelData(Y, Yname, W,
                                 covimages, dataAtPoints, basepixelarea)
      sumYloga <- Y$n * log(basepixelarea)
      serial  <- attr(df, "serial")
      Yserial <- attr(df, "Yserial")
    } else {
      ## fractional pixel areas
      pixsplit <- pixellate(splitwin, W)
      splitpixelarea <- as.vector(as.matrix(pixsplit))
      ## determine which points of Y are inside/outside window
      ins <- inside.owin(Y$x, Y$y, splitwin)
      ## split processing
      dfIN <- slrAssemblePixelData(Y[ins], Yname, W, covimages,
                                   dataAtPoints[ins, ], splitpixelarea)
      serialIN <- attr(dfIN, "serial")
      YserialIN   <- attr(dfIN, "Yserial")
      dfIN[[splitby]] <- TRUE
      dfOUT <- slrAssemblePixelData(Y[!ins], Yname, W, covimages,
                                    dataAtPoints[!ins, ],
                                    basepixelarea - splitpixelarea)
      serialOUT    <- attr(dfOUT, "serial")
      YserialOUT   <- attr(dfOUT, "Yserial")
      dfOUT[[splitby]] <- FALSE
      df <- rbind(dfIN, dfOUT)
      serial <- c(serialIN, serialOUT)
      Yserial <- c(YserialIN, YserialOUT)
      ## sum of log pixel areas associated with points
      Ysplit <- pixsplit[Y]
      sumYloga <- sum(log(ifelseXY(ins, Ysplit, basepixelarea - Ysplit)))
    }
  
    ## tack on any numeric values
    df <- do.call(cbind, append(list(df), numlist))
  
    ### RETURN ALL 
    Data <- list(response=Y,
                 covariates=covlist,
                 spatialnames=spatialnames,
                 numnames=numnames,
                 W=W,
                 df=df,
                 serial=serial,
                 Yserial=Yserial,
                 sumYloga=sumYloga,
                 dataAtPoints=dataAtPoints)
    return(Data)
  }

  getobj <- function(nama, env, dat) {
    if(!is.null(dat) && !is.null(x <- dat[[nama]]))
      return(x)
    else return(get(nama, envir=env))
  }

  convert <- function(x,W) {
    if(is.im(x) || is.function(x)) return(as.im(x,W))
    if(is.owin(x)) return(as.im(x, W, value=TRUE, na.replace=FALSE))
    return(NULL)
  }

  slr.prepare
})

## .............................................................

slrAssemblePixelData <- local({
  
  slrAssemblePixelData <- function(Y, Yname, W,
                                 covimages, dataAtPoints, pixelarea) {
    #' pixellate point pattern
    PY <- pixellate(Y, W=W, savemap=TRUE)
    IY <- eval.im(as.integer(PY>0))
    #'
    if(!is.null(dataAtPoints)) {
      #' overwrite pixel entries for data points using exact values
      #' spatial coordinates
      covimages[["x"]][Y] <- Y$x
      covimages[["y"]][Y] <- Y$y
      #' other values provided
      enames <- colnames(dataAtPoints)
      relevant <- enames %in% names(covimages)
      for(v in enames[relevant]) {
        cova <- covimages[[v]]
        cova[Y] <- dataAtPoints[, v, drop=TRUE]
        covimages[[v]] <- cova
      }
    }
    #' assemble list of all images
    Ylist <- list(IY)
    names(Ylist) <- Yname
    allimages <- append(Ylist, covimages)
    #' extract pixel values of each image, convert to data frame
    pixdata <- lapply(allimages, pixelvalues)
    df <- as.data.frame(pixdata)
    serial <- seq_len(nrow(df))
    ## add log(pixel area) column
    if(length(pixelarea) == 1) {
      df <- cbind(df, logpixelarea=log(pixelarea))
    } else {
      ok <- (pixelarea > 0)
      df <- cbind(df[ok, ], logpixelarea=log(pixelarea[ok]))
      serial <- serial[ok]
    }
    #' map rows of data frame to pixels
    attr(df, "serial") <- serial
    #' map original data points to pixels
    Yrowcol <- attr(PY, "map")
    attr(df, "Yserial") <- Yrowcol[,"row"] + (nrow(PY) - 1L) * Yrowcol[,"col"]
    return(df)
  }

  pixelvalues <- function(z) {
    v <- as.vector(as.matrix(z))
    if(z$type != "factor") return(v)
    lev <- levels(z)
    return(factor(v, levels=seq_along(lev), labels=lev))
  }

  slrAssemblePixelData

})

## ................. Methods ...................................


is.slrm <- function(x) {
  inherits(x, "slrm")
}

coef.slrm <- function(object, ...) {
  coef(object$Fit$FIT)
}

print.slrm <- function(x, ...) {
  lk <- x$CallInfo$link
  switch(lk,
         logit= {
           splat("Fitted spatial logistic regression model")
         },
         cloglog= {
           splat("Fitted spatial regression model (complementary log-log)")
         },
         {
           splat("Fitted spatial regression model")
           splat("Link =", dQuote(lk))
         })
  cat("Formula:\t")
  print(x$CallInfo$formula)
  splat("Fitted coefficients:")
  print(coef(x))
  return(invisible(NULL))
}

summary.slrm <- function(object, ...) {
  y <- object$CallInfo[c("link", "formula", "callstring")]
  co <- coef(object)
  se <- sqrt(diag(vcov(object)))
  two <- qnorm(0.975)
  lo <- co - two * se
  hi <- co + two * se
  zval <- co/se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  psig <- cut(pval, c(0, 0.001, 0.01, 0.05, 1),
              labels = c("***", "**", "*", "  "),
              include.lowest = TRUE)
  y$coefs.SE.CI <- data.frame(Estimate = co,
                              S.E.     = se, 
                              CI95.lo  = lo,
                              CI95.hi  = hi,
                              Ztest    = psig,
                              Zval     = zval)
  class(y) <- c(class(y), "summary.slrm")
  return(y)
}

print.summary.slrm <- function(x, ...) {
  switch(x$link,
         logit= {
           splat("Fitted spatial logistic regression model")
         },
         cloglog= {
           splat("Fitted spatial regression model (complementary log-log)")
         },
         {
           splat("Fitted spatial regression model")
           splat("Link =", dQuote(x$link))
         })
  cat("Call:\t")
  print(x$callstring)
  cat("Formula:\t")
  print(x$formula)
  splat("Fitted coefficients:\t")
  print(x$coefs.SE.CI)
  return(invisible(NULL))
}

coef.summary.slrm <- function(object, ...) { object$coefs.SE.CI }

logLik.slrm <- function(object, ..., adjust=TRUE) {
  FIT  <- object$Fit$FIT
  ll <- -deviance(FIT)/2
  if(adjust) {
    sumYloga <- object$Data$sumYloga
    ll <- ll - sumYloga
  }
  attr(ll, "df") <- length(coef(object))
  class(ll) <- "logLik"
  return(ll)
}

fitted.slrm <- function(object, ...) {
  if(length(list(...)) > 0)
    warning("second argument (and any subsequent arguments) ignored")
  predict(object, type="probabilities")
}

intensity.slrm <- function(X, ...) {
  Z <- predict(X, type="intensity", ..., newdata=NULL, window=NULL)
  if(is.stationary(X)) Z <- mean(Z)
  return(Z)
}

predict.slrm <- function(object, ..., type="intensity",
                         newdata=NULL, window=NULL) {
  type <- pickoption("type", type,
                     c(probabilities="probabilities",
                       link="link",
                       intensity="intensity",
                       lambda="intensity"))
  
  FIT     <- object$Fit$FIT
  link    <- object$CallInfo$link
  splitby <- object$CallInfo$splitby
  Yname   <- object$CallInfo$responsename
  W       <- object$Data$W
  df      <- object$Data$df
  loga    <- df$logpixelarea

  if(!is.null(window)) window <- as.owin(window)
  
  if(is.null(newdata) && is.null(window) && is.null(splitby)) {
    # fitted pixel values from existing fit
    switch(type,
           probabilities={
             values <- fitted(FIT)
           },
           link={
             values <- predict(FIT, type="link")
           },
           intensity={
             # this calculation applies whether an offset was included or not
             if(link == "cloglog") {
               linkvalues <- predict(FIT, type="link")
               values <- exp(linkvalues - loga)
             } else {
               probs <- fitted(FIT)
               values <- -log(1-probs)/exp(loga)
             }
           }
           )
    out <- im(values, xcol=W$xcol, yrow=W$yrow, unitname=unitname(W))
    return(out)
  } else {
    ## prediction from new data and/or at new locations
    if(is.null(newdata)) {
      ## prediction using existing covariates, at new locations
      newdata <- object$Data$covariates
    } else {
      ## prediction with completely new data
      stopifnot(is.list(newdata))
    }
    ## ensure newdata includes response pattern to placate internal code
    if(!(Yname %in% names(newdata)))
      newdata[[Yname]] <- ppp(window=window %orifnull% W)

    ## Update arguments that may affect pixel resolution
    CallInfo <- object$CallInfo
    CallInfo$dotargs <- resolve.defaults(list(...), CallInfo$dotargs)
    ## prevent pixel splitting
    CallInfo$splitby <- NULL

    ## process new data
    newData <- slr.prepare(CallInfo, environment(CallInfo$formula), newdata,
                           clip=!is.null(window))
    newdf   <- newData$df
    newW    <- newData$W
    newloga <- newdf$logpixelarea
    ## avoid NA etc
    npixel <- nrow(newdf)
    ok <- complete.cases(newdf)
    if(!all(ok)) {
      newdf   <- newdf[ok, , drop=FALSE]
      newloga <- newloga[ok]
    }
    ## compute link values
    linkvalues <- predict(FIT, newdata=newdf, type="link")
    ## transform to desired scale
    linkinv <- family(FIT)$linkinv
    switch(type,
           probabilities={
             values <- linkinv(linkvalues)
           },
           link={
             values <- linkvalues
           },
           intensity={
             # this calculation applies whether an offset was included or not
             if(link == "cloglog") {
               values <- exp(linkvalues - newloga)
             } else {
               probs <- linkinv(linkvalues)
               values <- -log(1-probs)/exp(newloga)
             }
           }
           )
    ## form image
    v <- rep.int(NA_real_, npixel)
    v[ok] <- values
    out <- im(v, xcol=newW$xcol, yrow=newW$yrow, unitname=unitname(W))
    return(out)
  }
}

plot.slrm <- function(x, ..., type="intensity") {
  xname <- short.deparse(substitute(x))
  y <- predict(x, type=type)
  dont.complain.about(y)
  do.call(plot.im, resolve.defaults(list(x=quote(y)), 
				    list(...), 
				    list(main=xname)))
}

formula.slrm <- function(x, ...) {
  f <- x$CallInfo$formula
  return(f)
}

terms.slrm <- function(x, ...) {
  terms(formula(x), ...)
}

labels.slrm <- function(object, ...) {
  # extract fitted trend coefficients
  co <- coef(object)
  # model terms
  tt <- terms(object)
  lab <- attr(tt, "term.labels")
  if(length(lab) == 0)
    return(character(0))
  # model matrix
  mm <- model.matrix(object)
  ass <- attr(mm, "assign")
  # 'ass' associates coefficients with model terms
  # except ass == 0 for the Intercept
  coef.ok <- is.finite(co)
  relevant <- (ass > 0) 
  okterms <- unique(ass[coef.ok & relevant])
  return(lab[okterms])
}

deviance.slrm <- function(object, ...) {
  deviance(object$Fit$FIT, ...)
}


extractAIC.slrm <- function (fit, scale = 0, k = 2, ...)
{
    edf <- length(coef(fit))
    aic <- AIC(fit)
    c(edf, aic + (k - 2) * edf)
}

model.frame.slrm <- function(formula, ...) {
  FIT <- formula$Fit$FIT
  mf <- model.frame(FIT, ...)
  return(mf)
}

model.matrix.slrm <- function(object,..., keepNA=TRUE) {
  FIT <- object$Fit$FIT
  mm <- model.matrix(FIT, ...)
  if(!keepNA)
    return(mm)
  df <- object$Data$df
  comp <- complete.cases(df)
  if(all(comp))
    return(mm)
  if(sum(comp) != nrow(mm))
      stop("Internal error in patching NA's")
  mmplus <- matrix(NA, nrow(df), ncol(mm))
  mmplus[comp, ] <- mm
  colnames(mmplus) <- colnames(mm)
  return(mmplus)
}

model.images.slrm <- function(object, ...) {
  mm <- model.matrix(object, ...)
  mm <- as.data.frame(mm)
  Data <- object$Data
  W      <- Data$W
  serial <- Data$serial
  splitby <- object$CallInfo$splitby
  blank   <- as.im(NA_real_, W)
  assignbyserial <- function(values, serial, template) {
    Z <- template
    Z$v[serial] <- values
    return(Z)
  }
  if(is.null(splitby)) {
    result <- lapply(as.list(mm), assignbyserial, serial=serial, template=blank)
  } else {
    df <- Data$df
    IN <- as.logical(df[[splitby]])
    OUT <- !IN
    mmIN <- mm[IN, , drop=FALSE]
    mmOUT <- mm[OUT, , drop=FALSE]
    resultIN <- lapply(as.list(mmIN), assignbyserial,
                       serial=serial[IN], template=blank)
    resultOUT <- lapply(as.list(mmOUT), assignbyserial,
                       serial=serial[OUT], template=blank)
    names(resultIN) <- paste(names(resultIN), splitby, "TRUE", sep="")
    names(resultOUT) <- paste(names(resultOUT), splitby, "FALSE", sep="")
    result <- c(resultIN, resultOUT)
  }
  return(as.solist(result))
}

update.slrm <- function(object, fmla, ..., evaluate=TRUE, env=parent.frame()) {
  z <- getCall(object)
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(fmla)) 
    z$formula <- update(formula(object), fmla)
  if(!missing(env))
    environment(z$formula) <- env
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(z)))
    for(a in names(extras)[existing]) z[[a]] <- extras[[a]]
    if (any(!existing)) {
      z <- c(as.list(z), extras[!existing])
      z <- as.call(z)
    }
  }
  if(evaluate) eval(z, envir=env) else z
}

updateData.slrm <- function(model, X, ...) {
  covlist <- model$Data$covariates
  Yname <- model$CallInfo$responsename
  newdata <- append(covlist, setNames(list(X), Yname))
  update(model, data=newdata)
}
  
anova.slrm <- local({

  anova.slrm <- function(object, ..., test=NULL) {
    objex <- append(list(object), list(...))
    if(!all(unlist(lapply(objex, is.slrm))))
      stop("Some arguments are not of class slrm")
    fitz <- lapply(objex, getFIT)
    do.call(anova, append(fitz, list(test=test)))
  }

  getFIT <- function(z) {z$Fit$FIT}

  anova.slrm
})

vcov.slrm <- function(object, ..., what=c("vcov", "corr", "fisher", "Fisher")) {
  stopifnot(is.slrm(object))
  what <- match.arg(what)
  vc <- vcov(object$Fit$FIT)
  result <- switch(what,
                   vcov = vc,
                   corr = {
                     sd <- sqrt(diag(vc))
                     vc / outer(sd, sd, "*")
                   },
                   fisher=,
                   Fisher={
                     solve(vc)
                   })
  return(result)
}

unitname.slrm <- function(x) {
  return(unitname(x$Data$response))
}

"unitname<-.slrm" <- function(x, value) {
  unitname(x$Data$response) <- value
  return(x)
}

domain.slrm <- Window.slrm <- function(X, ..., from=c("points", "covariates")) {
  from <- match.arg(from)
  as.owin(X, ..., from=from)
}

as.owin.slrm <- function(W, ..., from=c("points", "covariates")) {
  from <- match.arg(from)
  U <- switch(from,
              points     = W$Data$response,
              covariates = W$Data$W)
  V <- as.owin(U, ...)
  return(V)
}
  
is.stationary.slrm <- function(x) {
  trend <- rhs.of.formula(formula(x))
  return(identical.formulae(trend, ~1))
}

is.poisson.slrm <- function(x) { TRUE }

is.marked.slrm <- is.multitype.slrm <- function(X, ...) { FALSE }

reach.slrm <- function(x, ...) { 0 }

## pseudoR2.slrm is defined in ppmclass.R

Kmodel.slrm <- function(model, ...) { function(r) { pi * r^2 } }

pcfmodel.slrm <- function(model, ...) { function(r) { rep.int(1, length(r)) } }

parameters.slrm <- function(model, ...) { list(trend=coef(model)) }  

## ............ SIMULATION ..............................
  
simulate.slrm <- function(object, nsim=1, seed=NULL, ...,
                          window=NULL, covariates=NULL, 
                          verbose=TRUE, drop=FALSE) {
  check.1.integer(nsim)
  stopifnot(nsim >= 0)
  if(nsim == 0) return(simulationresult(list()))
  
  # .... copied from simulate.lm ....
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  starttime <- proc.time()

  # determine simulation window and compute intensity
  if(!is.null(window))
    stopifnot(is.owin(window))
  lambda <- predict(object, type="intensity", newdata=covariates, window=window)

  # max lambda (for efficiency)
  summ <- summary(lambda)
  lmax <- summ$max + 0.05 * diff(summ$range)

  # run
  out <- list()
  verbose <- verbose && (nsim > 1)
  if(verbose) {
    cat(paste("Generating", nsim, "simulations... "))
    pstate <- list()
  }
  for(i in 1:nsim) {
    out[[i]] <- rpoispp(lambda, lmax=lmax)
    if(verbose) pstate <- progressreport(i, nsim, state=pstate)
  }
  #' pack up
  out <- simulationresult(out, nsim, drop)
  out <- timed(out, starttime=starttime)
  attr(out, "seed") <- RNGstate
  return(out)
}

## ------------------ residuals --------------------------------

residuals.slrm <- function(object,
                           type=c("raw", "deviance", "pearson", "working", 
                                  "response", "partial", "score"),
                           ...) {
  type <- match.arg(type)
  otype <- if(type %in% c("raw", "score")) "response" else type
  FIT <- object$Fit$FIT  
  W <- object$Data$W
  res <- residuals(FIT, type=otype, ...)
  if(type == "score") {
    M <- model.matrix(object)
    res <- res * M
    colnames(res) <- colnames(M)
  }
  R <- wrangle2image(res, W)
  return(R)
}


## ------------------ leverage and influence -------------------


leverage.slrm <- function(model, ...) {
  slrmInfluence(model, "leverage", ...)[["leverage"]]
}

influence.slrm <- function(model, ...) {
  slrmInfluence(model, "influence", ...)[["influence"]]
}

dfbetas.slrm <- function(model, ...) {
  slrmInfluence(model, "dfbetas", ...)[["dfbetas"]]
}

dffit.slrm <- function(object, ...) {
  slrmInfluence(object, "dffit", ...)[["dffit"]]
}


slrmInfluence <- function(model,
                          what=c("all", "leverage", "influence",
                                 "dfbetas", "dffit"),
                          ...) {
  stopifnot(is.slrm(model))
  what <- match.arg(what, several.ok=TRUE)
  if("all" %in% what)
    what <- c("leverage", "influence", "dfbetas", "dffit")
  FIT <- model$Fit$FIT
  W <- model$Data$W
  ## nr <- nrow(W)
  ## nc <- ncol(W)
  result <- list()
  if("leverage" %in% what) {
    h <- hatvalues(FIT, ...)
    result$leverage <- wrangle2image(h, W)
  }
  if("influence" %in% what) {
    h <- hatvalues(FIT, ...)
    rP <- rstandard(FIT, type="pearson", ...)
    p <- length(coef(model))
    s <- (1/p) * rP^2 * h/(1-h)
    result$influence <- wrangle2image(s, W)
  }
  if("dfbetas" %in% what) {
    dfb <- dfbetas(FIT, ...)
    result$dfbetas <- wrangle2image(dfb, W)
  }
  if("dffit" %in% what) {
    dfb <- dfbeta(FIT, ...)  #sic
    X <- model.matrix(model) # sic
    if(is.null(dim(X)) || is.null(dim(dfb)) || !all(dim(X) == dim(dfb)))
      stop("Internal error: model.matrix dimensions incompatible with dfbeta")
    dff <- rowSums(X * dfb)
    result$dffit <- wrangle2image(dff, W)
  }

  return(result)
}

valid.slrm <- function(object, warn=TRUE, ...) {
  verifyclass(object, "slrm")
  coeffs <- coef(object)
  ok <- all(is.finite(coeffs))
  return(ok)
}
  
emend.slrm <- local({
  tracemessage <- function(depth, ...) {
    if(depth == 0) return(NULL)
    spacer <- paste(rep.int("  ", depth), collapse="")
    marker <- ngettext(depth, "trace", paste("trace", depth))
    marker <- paren(marker, "[")
    splat(paste0(spacer, marker, " ", paste(...)))
  }
  leaving <- function(depth) {
    tracemessage(depth, ngettext(depth, "Returning.", "Exiting level."))
  }
  emend.slrm <- function(object, ..., fatal=FALSE, trace=FALSE) {
    verifyclass(object, "slrm")
    fast <- spatstat.options("project.fast")
    # user specifies 'trace' as logical
    # but 'trace' can also be integer representing trace depth
    td <- as.integer(trace)
    trace <- (td > 0)
    tdnext <- if(trace) td+1 else 0
    if(valid.slrm(object)) {
      tracemessage(td, "Model is valid.")
      leaving(td)
      return(object)
    }
    # Fitted coefficients
    coef.orig <- coeffs <- coef(object)
    ## coefnames  <- names(coeffs)
    # Trend terms in trend formula
    trendterms <- attr(terms(object), "term.labels")
    # Mapping from coefficients to terms of GLM
    coef2term  <- attr(model.matrix(object), "assign")
    ## istrend <- (coef2term > 0)
    # Identify non-finite trend coefficients
    bad <-  !is.finite(coeffs)
    if(!any(bad)) {
      tracemessage(td, "Trend terms are valid.")
    } else {
      nbad <- sum(bad)
      tracemessage(td,
                   "Non-finite ",
                   ngettext(nbad,
                            "coefficient for term ",
                            "coefficients for terms "),
                   commasep(sQuote(trendterms[coef2term[bad]])))
      if(fast) {
        # remove first illegal term
        firstbad <- min(which(bad))
        badterm <- trendterms[coef2term[firstbad]]
        # remove this term from model
        tracemessage(td, "Removing term ", sQuote(badterm))
        removebad <- as.formula(paste("~ . - ", badterm), env=object$callframe)
        newobject <- update(object, removebad)
        if(trace) {
          tracemessage(td, "Updated model:")
          print(newobject)
        }
        # recurse
        newobject <- emend.slrm(newobject, fatal=fatal, trace=tdnext)
        # return
        leaving(td)
        return(newobject)
      } else {
        # consider all illegal terms
        bestobject <- NULL
        for(i in which(bad)) {
          badterm <- trendterms[coef2term[i]]
          # remove this term from model
          tracemessage(td, "Considering removing term ", sQuote(badterm))
          removebad <- as.formula(paste("~ . - ", badterm),
                                  env=object$callframe)
          object.i <- update(object, removebad)
          if(trace) {
            tracemessage(td, "Considering updated model:")
            print(object.i)
          }
          # recurse
          object.i <- emend.slrm(object.i, fatal=fatal, trace=tdnext)
          # evaluate log likelihood
          logL.i   <- logLik(object.i, warn=FALSE)
          tracemessage(td, "max log likelihood = ", logL.i)
          # optimise
          if(is.null(bestobject) || (logLik(bestobject, warn=FALSE) < logL.i))
            bestobject <- object.i
        }
        if(trace) {
          tracemessage(td, "Best submodel:")
          print(bestobject)
        }
        # return
        leaving(td)
        return(bestobject)
      }
    }
    object$projected <- TRUE
    object$coef.orig  <- coef.orig
    leaving(td)
    return(object)
  }
  emend.slrm
})

spatialCovariateUnderModel.slrm <- function(model, covariate, ...) {
  P <- predict(model, type="probabilities")
  loc <- rasterxy.im(P, drop=TRUE)
  Zvalues <- evaluateCovariateAtPoints(covariate, loc, ...)
  wtvalues <- P[drop=TRUE]
  df <- data.frame(Z=Zvalues, wt=wtvalues)
  return(df)
}
