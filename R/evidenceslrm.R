#'
#'   evidenceslrm.R
#'
#'   method for 'spatialCovariateEvidence' for class 'slrm'
#'
#'   $Revision: 1.6 $ $Date: 2022/05/20 10:13:08 $

spatialCovariateEvidence.slrm <- function(model, covariate, ...,
                           lambdatype=c("probabilities", "intensity"),
                           jitter=TRUE, jitterfactor=1,
                           modelname=NULL, covname=NULL,
                           dataname=NULL, subset=NULL) {
  lambdatype <- match.arg(lambdatype)
  #' trap misuse
  badargs <- intersect(c("eps", "dimyx"), names(list(...)))
  nbad <- length(badargs)
  if(nbad > 0)
    warning(paste(ngettext(nbad, "Argument", "Arguments"),
                  commasep(sQuote(badargs)),
                  ngettext(nbad, "is", "are"),
                  "ignored by rhohat.slrm"),
            call.=FALSE)
  #' evaluate covariate values at presence pixels and all pixels
  #' determine names
  if(is.null(modelname))
    modelname <- short.deparse(substitute(model))
  if(covNotNamed <- is.null(covname)) {
    covname <- singlestring(short.deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
  }
  if(is.null(dataname))
    dataname <- model$CallInfo$responsename

  csr <- is.stationary(model)
  
  info <-  list(modelname=modelname, covname=covname,
                dataname=dataname, csr=csr, ispois=TRUE,
                spacename="two dimensions")

  FIT  <- model$Fit$FIT
  link <- model$CallInfo$link
  ## original point pattern
  X <- model$Data$response
  W <- Window(X)

  ## extract data from each pixel (or split pixel)
  df <- model$Data$df
    
  ## restrict to subset if required
  if(!is.null(subset)) {
    ok <- inside.owin(df$x, df$y, subset)
    df <- df[ok, drop=FALSE]
    X <- X[subset]
    W <- W[subset, drop=FALSE]
  }

  ## presence/absence values
  responsename <- model$CallInfo$responsename
  presence <- as.logical(df[[responsename]])
  ## areas of pixels or split pixels
  pixelareas <- exp(df$logpixelarea)

  ## pixel centres as a point pattern
  P <- ppp(df$x, df$y, window=W)

  #' parse covariate argument
  if(is.character(covariate)) {
    #' One of the characters 'x' or 'y'
    #' Turn it into a function.
    ns <- length(covariate)
    if(ns == 0) stop("covariate is empty")
    if(ns > 1) stop("more than one covariate specified")
    covname <- covariate
    covNotNamed <- FALSE
    covariate <- switch(covname,
                        x=function(x,y) { x },
                        y=function(x,y) { y },
                        stop(paste("Unrecognised covariate",
                                   dQuote(covariate))))
  } 

  if(is.im(covariate)) {
    type <- "im"
    ZP <- safelookup(covariate, P)
    Z <- covariate[W, drop=FALSE]
    W <- as.owin(Z)
  } else if(is.function(covariate)) {
    type <- "function"
    ZP <- covariate(P$x, P$y)
    if(!all(is.finite(ZP)))
      warning("covariate function returned NA or Inf values")
    #' window
    W <- as.mask(W)
    #' covariate in window
    Z <- as.im(covariate, W=W)
    #' collapse function body to single string
    if(covNotNamed) covname <- singlestring(covname)
  } else if(is.null(covariate)) {
    stop("The covariate is NULL", call.=FALSE)
  } else stop(paste("The covariate should be",
                    "an image, a function(x,y)",
                    "or one of the characters",
                    sQuote("x"), "or", sQuote("y")),
              call.=FALSE)

  #' values of covariate at pixels or split pixels
  Zvalues <- ZP
  #'values of covariate at 'presence' pixels
  ZX <- Zvalues[presence]

  #' fitted probability/intensity values at all pixels or split pixels
  switch(lambdatype,
         probabilities = {
           lambda <- predict(FIT, newdata=df, type="response")
         },
         intensity = {
           if(link == "cloglog") {
             linkvalues <- predict(FIT, newdata=df, type="link")
             lambda <- exp(linkvalues)/pixelareas
           } else {
             probs <- predict(FIT, newdata=df, type="response")
             lambda <- -log(1-probs)/pixelareas
           }
         })

  #' apply jittering to avoid ties
  if(jitter) {
    ZX <- jitter(ZX, factor=jitterfactor)
    Zvalues <- jitter(Zvalues, factor=jitterfactor)
  }

  lambdaname <- paste("the fitted", lambdatype)
  check.finite(lambda, xname=lambdaname, usergiven=FALSE)
  check.finite(Zvalues, xname="the covariate", usergiven=TRUE)
  
  #' lambda values at data points
  lambdaX <- lambda[presence]

  #' lambda image(s)
  lambdaimage <- predict(model, window=W, type=lambdatype)
    
  #' wrap up 
  values <- list(Zimage      = Z,
                 lambdaimage = lambdaimage,
                 Zvalues     = Zvalues,
                 lambda      = lambda,
                 lambdaX     = lambdaX,
                 weights     = pixelareas,
                 ZX          = ZX,
                 type        = type)
  return(list(values=values, info=info))
}

