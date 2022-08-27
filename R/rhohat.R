#'
#'  rhohat.R
#'
#'  $Revision: 1.113 $  $Date: 2022/08/27 05:59:54 $
#'
#'  Non-parametric estimation of a function rho(z) determining
#'  the intensity function lambda(u) of a point process in terms of a
#'  spatial covariate Z(u) through lambda(u) = rho(Z(u)).
#'  More generally allows offsets etc.

#' Copyright (c) Adrian Baddeley 2015-2022
#' GNU Public Licence GPL >= 2.0

#' Code for generic rhohat() and rhohat.ppp() 
#' is now moved to spatstat.explore


rhohat.ppm <- function(object, covariate, ...,
                       weights=NULL,
                       method=c("ratio", "reweight", "transform"),
                       horvitz=FALSE,
                       smoother=c("kernel", "local",
                                  "decreasing", "increasing",
                                  "mountain", "valley",
                                  "piecewise"),
                       subset=NULL,
                       do.CI=TRUE,
                       jitter=TRUE, jitterfactor=1, interpolate=TRUE,
                       dimyx=NULL, eps=NULL,
                       n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
                       bwref=bw, covname, confidence=0.95,
                       positiveCI, breaks=NULL) {
  callstring <- short.deparse(sys.call())
  smoother <- match.arg(smoother)
  method <- match.arg(method)
  if(missing(positiveCI))
    positiveCI <- (smoother == "local")
  if(missing(covname)) 
    covname <- sensiblevarname(short.deparse(substitute(covariate)), "X")
  if(is.null(adjust))
    adjust <- 1

  if("baseline" %in% names(list(...)))
    warning("Argument 'baseline' ignored: not available for rhohat.ppm")

  ## validate model
  model <- object
  reference <- "model"
  modelcall <- model$call

  if(is.character(covariate) && length(covariate) == 1) {
    covname <- covariate
    switch(covname,
           x={
             covariate <- function(x,y) { x }
           }, 
           y={
             covariate <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
    covunits <- unitname(data.ppm(model))
  } else if(inherits(covariate, "distfun")) {
    covunits <- unitname(covariate)
  } else {
    covunits <- NULL
  }

  W <- Window(data.ppm(model))
  if(!is.null(subset)) W <- W[subset, drop=FALSE]
  areaW <- area(W)
  
  rhohatEngine(model, covariate, reference, areaW, ...,
               subset=subset,
               do.CI=do.CI,
               weights=weights,
               method=method,
               horvitz=horvitz,
               smoother=smoother,
               resolution=list(dimyx=dimyx, eps=eps),
               spatCovarArgs=list(clip.predict=FALSE,
                                  jitter=jitter,
                                  jitterfactor=jitterfactor,
                                  interpolate=interpolate),
               n=n, bw=bw, adjust=adjust, from=from, to=to,
               bwref=bwref, covname=covname, covunits=covunits,
               confidence=confidence, positiveCI=positiveCI,
               breaks=breaks,
               modelcall=modelcall, callstring=callstring)
}




#' Code for rhohat infrastructure is now moved to spatstat.explore
