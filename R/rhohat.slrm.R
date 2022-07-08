#'
#'      rhohat.slrm.R
#'
#'   $Revision: 1.5 $ $Date: 2022/07/08 13:39:13 $
#' 

rhohat.slrm <- function(object, covariate, ...,
                        weights=NULL,
                        method=c("ratio", "reweight", "transform"),
                        horvitz=FALSE,
                        smoother=c("kernel", "local",
                                   "decreasing", "increasing",
                                   "mountain", "valley",
                                   "piecewise"),
                        subset=NULL,
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
    warning("Argument 'baseline' ignored: not available for rhohat.slrm")

  ## validate model
  model <- object
  reference <- "model"
  modelcall <- model$call

  if(!is.null(object$CallInfo$splitby))
    stop("Sorry, rhohat.slrm is not yet implemented for split pixels",
         call.=FALSE)

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
    covunits <- unitname(response(model))
  } else {
    covunits <- NULL
  }

  W <- Window(response(model))
  if(!is.null(subset)) W <- W[subset, drop=FALSE]
  areaW <- area(W)
  
  rhohatEngine(model, covariate, reference, areaW, ...,
               spatCovarArgs=list(lambdatype="intensity"),
               weights=weights,
               method=method,
               horvitz=horvitz,
               smoother=smoother,
               n=n, bw=bw, adjust=adjust, from=from, to=to,
               bwref=bwref, covname=covname, covunits=covunits,
               confidence=confidence, positiveCI=positiveCI,
               breaks=breaks,
               modelcall=modelcall, callstring=callstring)
}
