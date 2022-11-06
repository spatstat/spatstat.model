## palmdiagnose.R
##
## Palm intensity diagnostic
## proposed by Tanaka, Ogata and Stoyan (2008)
##
## Copyright (c) Adrian Baddeley 2022
## GNU Public Licence >= 2.0
##
##  $Revision: 1.2 $ $Date: 2022/11/06 09:14:21 $

palmdiagnose <- function(object, ..., breaks=30, trim=30, rmax=Inf) {
  if(missing(object)) {
    models <- list(...)
    if(length(models) == 0) stop("No fitted models were supplied")
    if(!all(sapply(models, is.kppm)))
      stop("Each argument should be a kppm object, or a list of kppm objects")
  } else if(is.kppm(object)) {
    models <- list(object, ...)
    if(length(models) > 1 && !all(sapply(models, is.kppm)))
      stop("Each argument should be a kppm object, or a list of kppm objects")
  } else if(is.list(object) && all(sapply(object, is.kppm))) {
    models <- object
  } else stop("Argument 'object' should be a kppm object, or a list of kppm objects")
  ## must be stationary
  if(!all(sapply(models, is.stationary)))
    stop("Sorry, not yet implemented for inhomogeneous models")
  ## names
  nmodels <- length(models)
  if(sum(nzchar(names(models))) < nmodels) 
    names(models) <- if(nmodels == 1) "fit" else paste0("Model", 1:nmodels)
  modelnames <- names(models) <- make.names(names(models), unique=TRUE, allow_ = FALSE)
  ## extract point pattern
  Xlist <- unique(lapply(models, response))
  if(length(Xlist) > 1) stop("Models were fitted to different point patterns")
  X <- response(models[[1]])
  nX <- npoints(X)
  W <- Window(X)
  aW <- area(W)
  lamX <- nX/aW
  ## translation edge correction weights
  tran <- aW/setcov(W)
  tran <- eval.im(pmin(100, tran))
  tran <- tran/nX
  ## nonparametric estimate of Palm intensity
  r <- function(x,y) { sqrt(x^2 + y^2) }
  Z <- frypoints(X, dmax=rmax)
  R <- rhohat(Z, r, weights=tran, smoother="piecewise",
              breaks=breaks, from = 0, to = if(is.finite(rmax)) rmax else NULL)
  breaks <- attr(R, "stuff")$breaks
  ## replace 'ave' by \bar\lambda (using knowledge of Fry points)
  R$ave <- lamX
  ## rename function
  fnam <- c("lambda", "P")
  lmap <- c(rho=makefvlabel(NULL, "hat", fnam, "nonpar"),
            ave="bar(lambda)",
            var=makefvlabel("Var", "hat", fnam, "nonpar"),
            hi=makefvlabel(NULL, NULL, fnam, "hi"),
            lo=makefvlabel(NULL, NULL, fnam, "hi"))
  R <- rebadge.fv(R,
                  new.ylab=quote(lambda[P](r)),
                  new.fname=fnam,
                  tags=names(lmap),
                  new.labl=unname(lmap))
  R <- rebadge.fv(R,
                  tags=c("r", "rho"),
                  new.tags=c("r", "est"),
                  new.desc=c("Interpoint distance",
                             "Nonparametric estimate"))
  shadenames <- fvnames(R, ".s")
  dotnames <- fvnames(R, ".")
  ## parametric estimate(s)
  palmlam <- lapply(models,
                    function(fit, r, lam) { lam * (pcfmodel(fit))(r) },
                    r=R$r,
                    lam=lamX)
  if(nmodels == 1) {
    new.df <- data.frame(fit=palmlam[[1]])
    new.labl <-  makefvlabel(NULL, NULL, fnam, "fit")
    new.desc <-  "Parametric estimate"
  } else {
    new.df <- as.data.frame(palmlam)
    new.labl <- sapply(modelnames,
                       function(x, fn) makefvlabel(NULL, NULL, fn, x), fn=fnam)
    new.desc <- paste("Parametric estimate from", modelnames)
  }
  R <- bind.fv(R, new.df, new.labl, new.desc)
  ## tidy up
  dotnames <- fvnames(R, ".")  <- c(dotnames[1], modelnames, dotnames[-1])
  if(length(shadenames)) fvnames(R, ".s") <- shadenames
  attr(R, "alim") <- intersect.ranges(attr(R, "alim"),
                                      c(0, rmax.rule("K", W, lamX)))
  unitname(R) <- unitname(X)
  ## add a new class so that we can imitate the plot style of Tanaka et al
  class(R) <- union("palmdiag", class(R))
  attr(R, "breaks") <- breaks
  attr(R, "fitnames") <- setdiff(dotnames, c("est", "var", "hi", "lo"))
  return(R)
}

plot.palmdiag <- function(x, ..., style=c("intervals", "dots", "bands"),
                          args.dots=list(pch=16), args.intervals=list(),
                          main) {
  if(missing(main))
    main <- short.deparse(substitute(x))
  style <- match.arg(style)
  switch(style,
         bands = {
           z <- plot.fv(x, ..., main=main)
         },
         dots = ,
         intervals = {
           fitnames <- attr(x, "fitnames")
           fmla <- as.formula(paste("cbind", paren(paste(fitnames, collapse=", ")), "~ r"))
           z <- do.call(plot.fv,
                        resolve.defaults(list(quote(x)),
                                         list(fmla,
                                              shade=NULL,
                                              main=main),
                                         list(...),
                                         list(ylim=range(x, na.rm=TRUE))
                                         ))
           b <- attr(x, "breaks")
           rmid <- (b[-1] + b[-length(b)])/2
           f <- as.function(x, value=c("est", "lo", "hi"))
           ymid <- f(rmid)
           do.call(points, append(list(rmid, ymid), args.dots))
           if(style == "intervals") {
             yhi  <- f(rmid, "hi")
             ylo  <- f(rmid, "lo")
             do.call(segments,
                     append(list(rmid, ylo, rmid, yhi),
                            args.intervals))
           }
         })
  return(invisible(z))
}
