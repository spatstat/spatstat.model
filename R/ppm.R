#
#	$Revision: 1.64 $	$Date: 2025/09/11 03:56:53 $
#
#    ppm()
#          Fit a point process model to a two-dimensional point pattern
#
#

ppm <- function(Q, ...) {
  UseMethod("ppm")
}


ppm.formula <- function(Q, interaction=NULL, ..., data=NULL, subset) {
  ## remember call
  callstring <- short.deparse(sys.call())
  cl <- match.call()

  ## trap a common error to give a more informative message
  if(is.sob(data) || is.function(data)) 
    stop(paste("The argument", sQuote("data"),
               "should not be a spatial object;",
               "it should be a list of spatial objects"),
         call.=FALSE)
  
  ########### INTERPRET FORMULA ##############################
  
  if(!inherits(Q, "formula"))
    stop(paste("Argument 'Q' should be a formula"))
  formula <- Q
  
  ## check formula has LHS and RHS. Extract them
  if(length(formula) < 3)
    stop(paste("Formula must have a left hand side"))
  Yexpr <- formula[[2]]
  trend <- formula[c(1,3)]
  
  ## FIT #######################################
  thecall <- if(missing(subset)) {
    call("ppm", Q=Yexpr, trend=trend, data=data, interaction=interaction)
  } else {
    call("ppm", Q=Yexpr, trend=trend, data=data, interaction=interaction,
         subset=substitute(subset))
  }
  ncall <- length(thecall)
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    thecall[ncall + 1:nargh] <- argh
    names(thecall)[ncall + 1:nargh] <- names(argh)
  }
  callenv <- list2env(as.list(data), parent=parent.frame())
  result <- eval(thecall, envir=callenv)

  result$call <- cl
  result$callstring <- callstring
  result$callframe <- parent.frame()
  
  return(result)
}


ppm.quad <- ppm.ppp <- ppm.default <- 
function(Q,
         trend = ~1,
	 interaction = Poisson(),
         ..., 
         covariates = data,
         data = NULL,
         covfunargs = list(),
         subset,
         clipwin,
	 correction="border",
	 rbord = reach(interaction),
         use.gam=FALSE,
         method = c("mpl", "logi", "VBlogi"),
         forcefit=FALSE,
         improve.type = c("none", "ho", "enet"),
         improve.args=list(),
         emend=project,
         project=FALSE,
         prior.mean = NULL,
         prior.var = NULL,
         nd = NULL,
         eps = NULL,
         quad.args = list(),
         gcontrol=list(),
         nsim=100,
         nrmh=1e5,
         start=NULL,
         control=list(nrep=nrmh),
         verb=TRUE,
         callstring=NULL
) {
  Qname <- short.deparse(substitute(Q))

  if(is.NAobject(Q)) return(NAobject("ppm"))
  
  subsetexpr <- if(!missing(subset)) substitute(subset) else NULL
  clipwin    <- if(!missing(clipwin)) clipwin else NULL

  datalistname <- if(missing(covariates)) "data" else "covariates"

  if(!missing(emend) && !missing(project) && emend != project)
    stop("Conflicting options given: emend != project")

  ## Parse fitting method
  method.given <- !missing(method)
  improvetype.given <- !missing(improve.type)

  method <- match.arg(method)
  improve.type <- match.arg(improve.type)

  if(!is.null(prior.mean) | !is.null(prior.var)){
    if(!method.given) {
      method <- "VBlogi"
    } else if(method != "VBlogi") {
      stop(paste("Prior specification only works with method",
                 sQuote("VBlogi")))
    }
  }
  
  switch(method,
         mpl = { },
         logi = {
           VB <- FALSE
         },
         VBlogi = {
           method <- "logi"
           VB <- TRUE
         },
         ho = {
           ## old syntax
           method <- "mpl"
           if(!improvetype.given) {
             warning(paste("Syntax 'method=\"ho\"' is deprecated;",
                           "use 'improve.type=\"ho\"'"))
             improve.type <- "ho"
           } else if(improve.type != "ho") {
             stop(paste("Old syntax 'method=\"ho\"' is inconsistent with",
                        sQuote(paste("improve.type=", dQuote(improve.type)))))
           }
         },
         stop(paste("Unrecognised option", sQuote(paste("method=", dQuote(method))))))
             
           
  if(is.sob(covariates) || is.function(covariates))
    stop(paste("The argument", sQuote(datalistname),
               "should not be a spatial object;",
               "it should be a list of spatial objects"),
         call.=FALSE)
    
  if(inherits(Q, "logiquad")){
    if(missing(method))
      method <- "logi"
    if(method != "logi")
      stop(paste("Only method =", sQuote("logi"),
                 "makes sense when Q is of type", sQuote("logiquad")))
  }
  cl <- match.call()
  if(is.null(callstring)) 
    callstring <- paste(short.deparse(sys.call()), collapse="")

  if(is.ppp(Q) && is.marked(Q) && !is.multitype(Q)) 
    stop(paste("ppm is not yet implemented for marked point patterns,",
               "other than multitype patterns."))
  if(!(is.ppp(Q) || is.quad(Q) || checkfields(Q, c("data", "dummy")))) {
    stop("Argument Q must be a point pattern or a quadrature scheme")
  }
  X <- if(is.ppp(Q)) Q else Q$data

  ## Validate interaction
  if(is.null(interaction)) {
    interaction <- Poisson()
  } else if(inherits(interaction, "intermaker")) {
    ## e.g. 'interaction=Hardcore': invoke it without arguments
    interaction <- (f <- interaction)()
    dont.complain.about(f)
  } else if(!is.interact(interaction))
    stop("Argument 'interaction' must be an object of class 'interact'")
  
  ## Ensure interaction is fully defined  
  if(!is.null(ss <- interaction$selfstart)) {
    # invoke selfstart mechanism to fix all parameters
    interaction <- ss(X, interaction)
  }

  if(inherits(trend, "formula")) {
    ## handle "." in formula, representing all variables in 'data'
    if("." %in% variablesinformula(trend)) {
      if(is.null(covariates))
        stop("Cannot expand '.' since 'data' is not present", call.=FALSE)
      rhs <- paste(names(covariates), collapse=" + ")
      allmaineffects <- as.formula(paste("~", rhs))
      environment(allmaineffects) <- environment(trend)
      trend <- update(allmaineffects, trend)
    }
    ## expand polynom() in formula
    if(spatstat.options("expand.polynom"))
      trend <- expand.polynom(trend)
  }
  
  # validate choice of edge correction
  correction <- pickoption("correction", correction,
                           c(border="border",
                             periodic="periodic",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             trans="translate",
                             translate="translate",
                             translation="translate",
                             none="none"))
  
  # validate rbord 
  if(correction == "border") {
    # rbord for border correction
    rbord.given <- !missing(rbord) && !is.null(rbord)
    if(is.null(rbord))
      rbord <- reach(interaction)
    infin <- is.infinite(rbord)
    too.large <- infin || (eroded.areas(as.owin(X), rbord) == 0)
    if(too.large) {
      whinge <-
        paste(if(rbord.given) "rbord" else "the reach of this interaction",
              if(infin) "is infinite or unknown;"
              else "is too large for this window;",
              "please specify",
              if(rbord.given) "a smaller value of",
              "rbord, or use a different edge correction")
      stop(whinge)
    }
  } else {
    # rbord must be numeric to satisfy mpl.engine
    if(is.null(rbord))
      rbord <- 0
  }

  ##..............  Fit model ...............................
  switch(method,
         logi = {
           ## Fit by logistic composite likelihood
           fit <- logi.engine(Q=Q, trend=trend,
                              interaction=interaction,
                              covariates=covariates,
                              covfunargs=covfunargs,
                              subsetexpr=subsetexpr,
                              clipwin=clipwin,
                              correction=correction,
                              rbord=rbord,
                              use.gam=use.gam,
                              forcefit=forcefit,
                              nd = nd,
                              gcontrol=gcontrol,
                              callstring=callstring,
                              prior.mean=prior.mean,
                              prior.var=prior.var,
                              VB=VB,
                              quad.args = quad.args,
                              ...)
         },
         mpl = {
           ## fit by maximum pseudolikelihood
           fit <- mpl.engine(Q=Q, trend=trend,
                             interaction=interaction,
                             covariates=covariates,
                             covfunargs=covfunargs,
                             subsetexpr=subsetexpr,
                             clipwin=clipwin,
                             correction=correction,
                             rbord=rbord,
                             use.gam=use.gam,
                             forcefit=forcefit,
                             nd = nd,
                             eps = eps,
                             quad.args = quad.args,
                             gcontrol=gcontrol,
                             callstring=callstring,
                             ...)
         },
         stop(paste("Internal error - method =", sQuote(method)))
         )

  ## Fill in model details
  fit$Qname <- Qname
  if(!is.ppm(fit)) {
    ## internal use only - returns some other data
    return(fit)
  }
  fit$call <- cl
  fit$callstring <- callstring
  fit$callframe <- parent.frame()

  ## Detect invalid coefficients
  if(emend && !valid.ppm(fit))
    fit <- emend.ppm(fit)

  ##..............  Improve fit ...............................
  switch(improve.type,
         none = {
           fit$improve.type <- "none"
         },
         ho = {
           fit <- do.call(ho.engine,
                          resolve.defaults(list(quote(fit),
                                                nsim=nsim, nrmh=nrmh, start=start,
                                                control=control, verb=verb),
                                           improve.args))
         },
         enet = {
           fit <- do.call(enet.engine,
                          append(list(quote(fit)), improve.args))
         })

  if(emend && !valid.ppm(fit))
    fit <- emend.ppm(fit)
  
  return(fit)
}

