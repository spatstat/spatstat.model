#
# mppm.R
#
#  $Revision: 1.108 $   $Date: 2025/06/26 01:52:24 $
#

mppm <- local({

  mppm <- function(formula, data, interaction=Poisson(), ...,
                   iformula=NULL,
                   random=NULL,
                   weights=NULL,
                   use.gam=FALSE,
                   reltol.pql=1e-3,
                   gcontrol=list()
                   ) {
    ## remember call
    cl <- match.call()
    callstring <- paste(short.deparse(sys.call()), collapse="")

    ## Validate arguments
    if(!inherits(formula, "formula"))
      stop(paste("Argument", dQuote("formula"), "should be a formula"))
    stopifnot(is.hyperframe(data))
    data.sumry <- summary(data, brief=TRUE)
    npat <- data.sumry$ncases
    if(npat == 0)
      stop(paste("Hyperframe", sQuote("data"), "has zero rows"))
    if(!is.null(iformula) && !inherits(iformula, "formula"))
      stop(paste("Argument", sQuote("iformula"), "should be a formula or NULL"))
    if(has.random <- !is.null(random)) {
      if(!inherits(random, "formula"))
        stop(paste(sQuote("random"), "should be a formula or NULL"))
      if(use.gam)
        stop("Sorry, random effects are not available in GAMs")
    }
    if(! (is.interact(interaction) || is.hyperframe(interaction)))
      stop(paste("The argument", sQuote("interaction"),
                 "should be a point process interaction object (class",
                 dQuote("interact"), 
                 "), or a hyperframe containing such objects", sep=""))

    if(is.null(weights)) {
      weights <- rep(1, npat)
    } else {
      check.nvector(weights, npat, things="rows of data", oneok=TRUE, vname="weights")
      if(length(weights) == 1L) weights <- rep(weights, npat)
    } 
    
    backdoor <- list(...)$backdoor
    if(is.null(backdoor) || !is.logical(backdoor))
      backdoor <- FALSE

    ############## HANDLE FORMULAS ############################
  
    ##------  Trend Formula ------------------
    
    ## check all variables in trend formula are recognised
    checkvars(formula, data.sumry$col.names,
              extra=c("x","y","id","marks"), bname="data")
    ## check formula has LHS and RHS. Extract them
    if(length(formula) < 3)
      stop(paste("Argument", sQuote("formula"),
                 "must have a left hand side"))
    Yname <- formula[[2]]
    trend <- formula[c(1,3)]
    if(!is.name(Yname))
      stop("Left hand side of formula should be a single name")
    Yname <- paste(Yname)
    if(!inherits(trend, "formula"))
      stop("Internal error: failed to extract RHS of formula")
    allvars <- variablesinformula(trend)
  
    ## --- Interaction formula -----
    
    ## names of interactions as they may appear in formulae
    itags <- 
      if(is.hyperframe(interaction)) names(interaction) else "Interaction"
    ninteract <- length(itags)
    ## ensure `iformula' is a formula without a LHS
    ## and determine which columns of `interaction' are actually used
    if(is.null(iformula)) {
      if(ninteract > 1)
        stop(paste("interaction hyperframe has more than 1 column;",
                   "you must specify the choice of interaction",
                   "using argument",  sQuote("iformula")))
      iused <- TRUE
      iformula <-  as.formula(paste("~", itags))
    } else {
      if(length(iformula) > 2)
        stop(paste("The interaction formula",
                   sQuote("iformula"),
                   "should not have a left hand side"))
      ## valid variables in `iformula' are interactions and data frame columns
      permitted <- paste(sQuote("interaction"),
                         "or permitted name in", sQuote("data"))
      checkvars(iformula, itags, extra=c(data.sumry$dfnames, "id"),
                bname=permitted)
      ivars <- variablesinformula(iformula)
      ## check which columns of `interaction' are actually used
      iused <- itags %in% ivars
      if(sum(iused) == 0)
        stop("No interaction specified in iformula")
      ## OK
      allvars <- c(allvars, ivars)
    } 

  
    ## --- Random effects formula ----
    if(!is.null(random))  {
      if(length(random) > 2)
        stop(paste("The random effects formula",
                   sQuote("random"),
                   "should not have a left hand side"))
      checkvars(random, itags, extra=c(data.sumry$col.names, "x", "y", "id"),
                bname="either data or interaction")
      allvars <- c(allvars, variablesinformula(random))
    }

  
    ## ---- variables required (on RHS of one of the above formulae) -----
    allvars <- unique(allvars)

  
    ########  EXTRACT DATA  #####################################
  
    ## Insert extra variable 'id'
    data <- cbind.hyperframe(data, id=factor(1:npat))
    data.sumry <- summary(data, brief=TRUE)
    allvars <- unique(c(allvars, "id"))

    ## Extract the list of responses (point pattern/quadscheme)
    Y <- data[, Yname, drop=TRUE]
    if(npat == 1) Y <- solist(Y)
    Yclass <- data.sumry$classes[Yname]
    if(Yclass == "ppp") {
      ## convert to quadrature schemes, for efficiency's sake
      Y <- solapply(Y, quadschemeplus, ...)
      ## Ydescrip <- "point patterns" ## not used
    } else if(Yclass == "quad") {
      Y <- as.solist(Y)
      ## Ydescrip <- "quadrature schemes" ## not used
    } else {
      stop(paste("Column", dQuote(Yname), "of data",
                 "does not consist of point patterns (class ppp)",
                 "nor of quadrature schemes (class quad)"),
           call.=FALSE)
    }
    ## Extract sub-hyperframe of data named in formulae
    datanames <- names(data)
    used.cov.names <- allvars[allvars %in% datanames]
    has.covar <- (length(used.cov.names) > 0) 
    if(has.covar) {
      dfvar <- used.cov.names %in% data.sumry$dfnames
      imvar <- data.sumry$types[used.cov.names] == "im"
      if(any(nbg <- !(dfvar | imvar)))
        stop(paste("Inappropriate format for",
                   ngettext(sum(nbg), "covariate", "covariates"),
                   paste(sQuote(used.cov.names[nbg]), collapse=", "),
                   ": should contain image objects or vector/factor"))
      covariates.hf <- data[, used.cov.names, drop=FALSE]
      has.design <- any(dfvar)
      dfvarnames <- used.cov.names[dfvar]
      datadf <-
        if(has.design)
          as.data.frame(covariates.hf, discard=TRUE, warn=FALSE)
        else NULL
      if(has.design) {
        ## check for NA's in design covariates
#        if(any(nbg <- apply(is.na(datadf), 2, any)))
        if(any(nbg <- matcolany(is.na(datadf))))
          stop(paste("There are NA's in the",
                     ngettext(sum(nbg), "covariate", "covariates"),
                     commasep(dQuote(names(datadf)[nbg]))))
      }
    } else {
      has.design <- FALSE
      datadf     <- NULL
    }
  
    ############### INTERACTION ###################################
    ## ensure `interaction' is a hyperframe of `interact' objects
    ## with the right number of rows.
    ## All entries in a column must represent the same process
    ## (possibly with different values of the irregular parameters).
    ## Extract the names of the point processes.
    if(is.interact(interaction)) {
      ninteract <- 1
      processes <- list(Interaction=interaction$name)
      interaction <- hyperframe(Interaction=interaction, id=1:npat)[,1]
      constant <- c(Interaction=TRUE)
    } else if(is.hyperframe(interaction)) {
      inter.sumry <- summary(interaction)
      ninteract <- inter.sumry$nvars
      ## ensure it has the same number of rows as 'data'
      nr <- inter.sumry$ncases
      if(nr == 1 && npat > 1) {
        interaction <- cbind.hyperframe(id=1:npat, interaction)[,-1]
        inter.sumry <- summary(interaction)
      } else if(nr != npat)
        stop(paste("Number of rows in", sQuote("interaction"),
                   "=", nr, "!=", npat, "=",
                   "number of rows in", sQuote("data")))
      ## check all columns contain interaction objects
      ok <- (inter.sumry$classes == "interact")
      if(!all(ok)) {
        nbg <- names(interaction)[!ok]
        nn <- sum(!ok)
        stop(paste(ngettext(nn, "Column", "Columns"),
                   paste(sQuote(nbg), collapse=", "),
                   ngettext(nn, "does", "do"),
                   "not consist of interaction objects"))
      }
      ## all entries in a column must represent the same process type
      ## (with possibly different values of the irregular parameters)
      ok <- unlist(lapply(as.list(interaction), consistentname))
      if(!all(ok)) {
        nbg <- names(interaction)[!ok]
        stop(paste("Different interactions may not appear in a single column.",
                   "Violated by",
                   paste(sQuote(nbg), collapse=", ")))
      }
      processes <- lapply(as.list(interaction), firstname)
    
      ## determine whether all entries in a column are EXACTLY the same
      ## (=> have the same parameters)
      constant <- (inter.sumry$storage == "hyperatom")
      if(any(!constant)) {
        others <- interaction[,!constant]
        constant[!constant] <- sapply(lapply(as.list(others), unique),
                                      length) == 1
      }
    }
    ## check for trivial (Poisson) interactions
    trivial <- unlist(lapply(as.list(interaction), allpoisson))
    
    ## check that iformula does not combine two interactions on one row
    nondfnames <- datanames[!(datanames %in% data.sumry$dfnames)]
    ip <- impliedpresence(itags, iformula, datadf, nondfnames)
    if(any(rowSums(ip) > 1))
      stop("iformula invokes more than one interaction on a single row")
  
    ##
    #################### BERMAN-TURNER DEVICE #########################
    ##
    ## set up list to contain the glm variable names for each interaction.
    Vnamelist <- rep(list(character(0)), ninteract)
    names(Vnamelist) <- itags
    ## set up list to contain 'IsOffset'
    Isoffsetlist <- rep(list(logical(0)), ninteract)
    names(Isoffsetlist) <- itags
    ####
    ## ---------------- L O O P ---------------------------------
    for(i in 1:npat) {
      ## extract responses and covariates for presentation to ppm()
      Yi <- Y[[i]]
      covariates <-
        if(has.covar) covariates.hf[i, , drop=TRUE, strip=FALSE] else NULL
      if(has.design) {
        ## convert each data frame value to an image
        covariates[dfvarnames] <-
          lapply(as.list(as.data.frame(covariates[dfvarnames])),
                 as.im, W=Yi$data$window)
      }
    
      ## Generate data frame and glm info for this point pattern
      ## First the trend covariates
      prep0 <- bt.frame(Yi, trend, Poisson(), ..., covariates=covariates,
                        allcovar=TRUE, use.gam=use.gam)
      glmdat <- prep0$glmdata
    
      ## now the nontrivial interaction terms
      for(j in (1:ninteract)[iused & !trivial]) {
        inter <- interaction[i,j,drop=TRUE]
        if(!is.null(ss <- inter$selfstart)) 
          interaction[i,j] <- inter <- ss(Yi$data, inter)
        prepj <- bt.frame(Yi, ~1, inter, ..., covariates=covariates,
                          allcovar=TRUE, use.gam=use.gam,
                          vnamebase=itags[j], vnameprefix=itags[j])
        ## store GLM variable names & check consistency
        vnameij <- prepj$Vnames
        if(i == 1)
          Vnamelist[[j]] <- vnameij
        else if(!identical(vnameij, Vnamelist[[j]]))
          stop("Internal error: Unexpected conflict in glm variable names")
        ## store offset indicator vectors
        isoffset.ij <- prepj$IsOffset
        if(i == 1)
          Isoffsetlist[[j]] <- isoffset.ij
        else if(!identical(isoffset.ij, Isoffsetlist[[j]]))
          stop("Internal error: Unexpected conflict in offset indicators")
        ## GLM data frame for this interaction
        glmdatj <- prepj$glmdata
        if(nrow(glmdatj) != nrow(glmdat))
          stop("Internal error: differing numbers of rows in glm data frame")
        iterms.ij <- glmdatj[vnameij]
        subset.ij <- glmdatj$.mpl.SUBSET
        ## tack on columns of interaction terms
        glmdat <- cbind(glmdat, iterms.ij)
        ## update subset (quadrature points where cif is positive)
        glmdat$.mpl.SUBSET <- glmdat$.mpl.SUBSET & subset.ij
      }

      ## assemble the Mother Of All Data Frames
      if(i == 1) {
        moadf <- glmdat
      } else {
        ## There may be new or missing columns
        recognised <- names(glmdat) %in% names(moadf)
        if(any(!recognised)) {
          newnames <- names(glmdat)[!recognised]
          zeroes <- as.data.frame(matrix(0, nrow(moadf), length(newnames)))
          names(zeroes) <- newnames
          moadf <- cbind(moadf, zeroes)
        }
        provided   <- names(moadf)  %in% names(glmdat)
        if(any(!provided)) {
          absentnames <- names(moadf)[!provided]
          zeroes <- as.data.frame(matrix(0, nrow(glmdat), length(absentnames)))
          names(zeroes) <- absentnames
          glmdat <- cbind(glmdat, zeroes)
        }
        ## Ensure factor columns are consistent
        m.isfac <- sapply(as.list(glmdat), is.factor)
        g.isfac <- sapply(as.list(glmdat), is.factor)
        if(any(uhoh <- (m.isfac != g.isfac)))
          errorInconsistentRows("values (factor and non-factor)",
                                colnames(moadf)[uhoh])
        if(any(m.isfac)) {
          m.levels <- lapply(as.list(moadf)[m.isfac], levels)
          g.levels <- lapply(as.list(glmdat)[g.isfac], levels)
          clash <- !mapply(identical, x=m.levels, y=g.levels)
          if(any(clash))
            errorInconsistentRows("factor levels",
                                  (colnames(moadf)[m.isfac])[clash])
        }
        ## Finally they are compatible
        moadf <- rbind(moadf, glmdat)
      }
    }
    ## ---------------- E N D   o f    L O O P  --------------------------
    ##
    ## add case weights
    moadf$caseweight <- weights[moadf$id]
    ##
    ## backdoor exit - Berman-Turner frame only - used by predict.mppm 
    if(backdoor)
      return(moadf)
    ##
    ##
    ## --------------------------------------------------------------------
    ## 
    ## Construct the glm formula for the Berman-Turner device
    ## 
    ## Get trend part from the last-computed prep0
    fmla  <- prep0$trendfmla
    ## Tack on the RHS of the interaction formula
    if(!all(trivial))
      fmla <- paste(fmla, "+", as.character(iformula)[[2]])
    ## Make it a formula
    fmla <- as.formula(fmla)

    ## Ensure that each interaction name is recognised.
    ##
    ## To the user, an interaction is identified by its `tag' name
    ## (default tag: "Interaction")
    ##
    ## Internally, an interaction is fitted using its sufficient statistic
    ## which may be 0, 1 or k-dimensional. 
    ## The column names of the sufficient statistic are the Vnames
    ## returned from ppm.
    ## The Poisson process is a special case: it is 0-dimensional (no Vnames).
    ##
    ## For k-dimensional sufficient statistics, we modify the formulae,
    ## replacing the interaction name by (vname1 + vname2 + .... + vnamek)
    ## 
    for(j in (1:ninteract)[iused]) {
      vnames <- Vnamelist[[j]]
      tag    <- itags[j]
      isoffset <- Isoffsetlist[[j]]
      if(any(isoffset)) {
        ## enclose names of offset variables in 'offset()'
        vnames[isoffset] <- paste("offset(", vnames[isoffset], ")", sep="")
      }
      if(trivial[j]) 
        ## Poisson case: add a column of zeroes
        moadf[[tag]] <- 0
      else if(!identical(vnames, tag)) {
        if(length(vnames) == 1) 
          ## tag to be replaced by vname
          vn <- paste("~", vnames[1])
        else 
          ## tag to be replaced by (vname1 + vname2 + .... + vnamek)
          vn <- paste("~(", paste(vnames, collapse=" + "), ")")
        ## pull out formula representation of RHS
        vnr <- as.formula(vn)[[2]]
        ## make substitution rule: list(<tag>=<vnr>)
        vnsub <- list(vnr)
        names(vnsub) <- tag
        ## perform substitution in trend formula
        fmla <- eval(substitute(substitute(fom, vnsub), list(fom=fmla)))
        ## perform substitution in random effects formula
        if(has.random && tag %in% variablesinformula(random))
          random <- eval(substitute(substitute(fom, vnsub), list(fom=random)))
      }
    }

    fmla <- as.formula(fmla)
    ## Fix scoping problem
    assign("glmmsubset", moadf$.mpl.SUBSET, envir=environment(fmla))
    for(nama in colnames(moadf))
      assign(nama,       moadf[[nama]],     envir=environment(fmla))
    ## Satisfy package checker
    glmmsubset <- .mpl.SUBSET <- moadf$.mpl.SUBSET
    .mpl.W      <- moadf$.mpl.W
    caseweight  <- moadf$caseweight
    
    ## ---------------- FIT THE  ------------------------------------
    want.trend <- prep0$info$want.trend
    if(want.trend && use.gam) {
      fitter <- "gam"
      ctrl <- do.call(gam.control, resolve.defaults(gcontrol, list(maxit=50)))
      FIT  <- gam(fmla, family=quasi(link=log, variance=mu),
                  weights=.mpl.W * caseweight,
                  data=moadf, subset=(.mpl.SUBSET=="TRUE"),
                  control=ctrl)
      deviants <- deviance(FIT)
    } else if(!is.null(random)) {
      fitter <- "glmmPQL"
      ctrl <- do.call(lmeControl, resolve.defaults(gcontrol, list(maxIter=50)))
      attr(fmla, "ctrl") <- ctrl # very strange way to pass argument
      fixed <- 42 # to satisfy package checker
      FIT  <- hackglmmPQL(fmla, random=random,
                          family=quasi(link=log, variance=mu),
                          weights=.mpl.W * caseweight,
                          data=moadf, subset=glmmsubset,
                          control=attr(fixed, "ctrl"),
                          reltol=reltol.pql)
      deviants <-  -2 * logLik(FIT)
    } else {
      fitter <- "glm"
      ctrl <- do.call(glm.control, resolve.defaults(gcontrol, list(maxit=50)))
      FIT  <- glm(fmla, family=quasi(link="log", variance="mu"),
                  weights=.mpl.W * caseweight,
                  data=moadf, subset=(.mpl.SUBSET=="TRUE"),
                  control=ctrl)
      deviants <- deviance(FIT)
    }
    env <- list2env(moadf, parent=sys.frame(sys.nframe()))
    environment(FIT$terms) <- env
    ## maximised log-pseudolikelihood
    W <- with(moadf, .mpl.W * caseweight)
    SUBSET <- moadf$.mpl.SUBSET
    Z <- (moadf$.mpl.Y != 0)
    maxlogpl <- -(deviants/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
    ##
    ## ---------------- PACK UP THE RESULT --------------------------------
    ##
    result <- list(Call = list(callstring=callstring, cl=cl),
                   Info =
                   list(
                     has.random=has.random,
                     has.covar=has.covar,
                     has.design=has.design,
                     Yname=Yname,
                     used.cov.names=used.cov.names,
                     allvars=allvars,
                     names.data=names(data),
                     is.df.column=(data.sumry$storage == "dfcolumn"),
                     rownames=row.names(data),
                     correction=prep0$info$correction,
                     rbord=prep0$info$rbord
                     ),
                   Fit=
                   list(
                     fitter=fitter,
                     use.gam=use.gam,
                     fmla=fmla,
                     FIT=FIT,
                     moadf=moadf,
                     Vnamelist=Vnamelist,
                     Isoffsetlist=Isoffsetlist
                     ),
                   Inter =
                   list(
                     ninteract=ninteract,
                     interaction=interaction,
                     iformula=iformula,
                     iused=iused,
                     itags=itags,
                     processes=processes,
                     trivial=trivial,
                     constant=constant
                     ),
                   formula=formula,
                   trend=trend,
                   iformula=iformula,
                   random=random,
                   npat=npat,
                   data=data,
                   Y=Y,
                   maxlogpl=maxlogpl,
                   datadf=datadf)

    class(result) <- c("mppm", class(result))
    return(result)
  }

  # helper functions
  checkvars <- function(f, b, extra=NULL, bname=short.deparse(substitute(b))){
    fname <- short.deparse(substitute(f))
    fvars <- variablesinformula(f)
    bvars <- if(is.character(b)) b else names(b)
    bvars <- c(bvars, extra)
    nbg <- !(fvars %in% bvars)
    if(any(nbg)) {
      nn <- sum(nbg)
      stop(paste(ngettext(nn, "Variable", "Variables"),
                 commasep(dQuote(fvars[nbg])),
                 "in", fname,
                 ngettext(nn, "is not one of the", "are not"),
                 "names in", bname))
    }
    return(NULL)
  }
  
  consistentname <- function(x) {
    xnames <- unlist(lapply(x, getElement, name="name"))
    return(length(unique(xnames)) == 1)
  }

  firstname <- function(z) { z[[1]]$name }

  allpoisson <- function(x) all(sapply(x, is.poisson.interact))

  marklevels <- function(x) { levels(marks(x)) }
  
  errorInconsistentRows <- function(what, offending) {
    stop(paste("There are inconsistent",
               what,
               "for the",
               ngettext(length(offending), "variable", "variables"),
               commasep(sQuote(offending)),
               "between different rows of the hyperframe 'data'"),
         call.=FALSE)
  }

  quadschemeplus <- function(..., quad.args) {
    #' catch argument 'quad.args' recognised by ppm
    if(missing(quad.args)) quadscheme(...) else do.call(quadscheme, append(list(...), quad.args))
  }

  
  mppm
})


is.mppm <- function(x) {
  inherits(x, "mppm")
}

coef.mppm <- function(object, ...) {
  coef(object$Fit$FIT)
}

getglmdata.mppm <- function(object, ...) {
  verifyclass(object, "mppm")
  return(object$Fit$moadf)
}

getglmfit.mppm <- function(object, ...) {
  verifyclass(object, "mppm")
  return(object$Fit$FIT)
}

getglmsubset.mppm <- function(object, ...) {
  verifyclass(object, "mppm")
  return(object$Fit$moadf$.mpl.SUBSET)
}

hasglmfit.mppm <- function(object) { TRUE }



fixef.mppm <- function(object, ...) {
  if(object$Fit$fitter == "glmmPQL")
    fixef(object$Fit$FIT)
  else
    coef(object$Fit$FIT)
}

ranef.mppm <- function(object, ...) {
  if(object$Fit$fitter == "glmmPQL")
    ranef(object$Fit$FIT)
  else
    as.data.frame(matrix(, nrow=object$npat, ncol=0))
}


print.mppm <- function(x, ...) {
  print(summary(x, ..., brief=TRUE))
}

is.poisson.mppm <- function(x) {
  trivial <- x$Inter$trivial
  iused <- x$Inter$iused
  all(trivial[iused])
}

quad.mppm <- function(x) {
  as.solist(x$Y)
}

data.mppm <- function(x) {
  solapply(x$Y, getElement, name="data")
}

is.marked.mppm <- function(X, ...) {
  any(sapply(data.mppm(X), is.marked))
}
  
is.multitype.mppm <- function(X, ...) {
  any(sapply(data.mppm(X), is.multitype))
}
  
windows.mppm <- function(x) {
  solapply(data.mppm(x), Window)
}

logLik.mppm <- function(object, ..., warn=TRUE) {
  if(warn && !is.poisson.mppm(object))
    warning(paste("log likelihood is not available for non-Poisson model;",
                  "log-pseudolikelihood returned"))
  ll <- object$maxlogpl
  attr(ll, "df") <- length(fixef(object))
  class(ll) <- "logLik"
  return(ll)
}

AIC.mppm <- function(object, ..., k=2, takeuchi=TRUE) {
  ll <- logLik(object, warn=FALSE)
  pen <- attr(ll, "df")
  if(takeuchi && !is.poisson(object)) {
    vv <- vcov(object, what="all")
    J  <- vv$fisher
    H  <- vv$internals$A1
    ## Takeuchi penalty = trace of J H^{-1} = trace of H^{-1} J
    JiH <- try(solve(H, J), silent=TRUE)
    if(!inherits(JiH, "try-error")) 
      pen <- sum(diag(JiH))
  } 
  return(- 2 * as.numeric(ll) + k * pen)
}

extractAIC.mppm <- function(fit, scale = 0, k = 2, ..., takeuchi = TRUE) 
{
  edf <- length(coef(fit))
  aic <- AIC(fit, k = k, takeuchi = takeuchi)
  c(edf, aic)
}

getCall.mppm <- function(x, ...) { x$Call$cl }

terms.mppm <- function(x, ...) { terms(formula(x)) }

nobs.mppm <- function(object, ...) { sum(sapply(data.mppm(object), npoints)) }

simulate.mppm <- function(object, nsim=1, ..., verbose=TRUE) {
  subs <- subfits(object)
  nr <- length(subs)
  sims <- list()
  if(verbose) {
    splat("Generating simulated realisations of", nr, "models..")
    state <- list()
  }
  for(irow in seq_len(nr)) {
    model.i <- subs[[irow]]
    dont.complain.about(model.i)
    sims[[irow]] <- do.call(simulate,
                            resolve.defaults(list(object=quote(model.i),
                                                  nsim=nsim, drop=FALSE),
                                             list(...),
                                             list(progress=FALSE)))
    if(verbose) state <- progressreport(irow, nr, state=state)
  }
  sim1list <- lapply(sims, "[[", i=1)
  h <- hyperframe("Sim1"=sim1list)
  if(nsim > 1) {
    for(j in 2:nsim) {
      simjlist <- lapply(sims, "[[", i=j)
      hj <- hyperframe(Sim=simjlist)
      names(hj) <- paste0("Sim", j)
      h <- cbind(h, hj)
    }
  }
  return(h)
}

model.matrix.mppm <- function(object, ..., keepNA=TRUE, separate=FALSE) {
  df <- object$Fit$moadf
  FIT <- object$Fit$FIT
  environment(FIT) <- list2env(df)
  #' track subset
  ok <- complete.cases(df) & df$.mpl.SUBSET
  nok <- sum(ok)
  nfull <- nrow(df)
  #'
  if(keepNA) {
    mm <- model.matrix(FIT, ..., subset=NULL, na.action=NULL)
    nr <- nrow(mm)
    if(nr != nfull) {
      if(nr == nok) {
        ## some methods for model.matrix ignore 'subset=NULL'
        mmfull <- matrix(NA, nfull, ncol(mm), dimnames=list(NULL, colnames(mm)))
        mmfull[ok,] <- mm
        mm <- mmfull
      } else {
        stop(paste("Internal error: model matrix has wrong number of rows:",
                   nr, "should be", nfull, "or", nok),
             call.=FALSE)
      }
    }
  } else {
    mm <- model.matrix(FIT, ...)
    if(nrow(mm) != nok)
      stop("Internal error: model matrix has wrong number of rows", call.=FALSE)
    df <- df[ok, , drop=FALSE]
  }

  ## remove NaN
  if(any(isnan <- is.nan(mm)))
     mm[isnan] <- 0
  
  ## get standard attributes of model matrix
  ctr <- attr(mm, "contrasts")
  if(is.null(ctr) && !is.null(ctr <- FIT[["contrasts"]]))
    attr(mm, "contrasts") <- ctr
  ass <- attr(mm, "assign")
  if(is.null(ass)) {
    ass <- FIT[["assign"]]
    if(length(ass) == ncol(mm)) {
      attr(mm, "assign") <- ass
    } else {
      ass <- NULL
      warning("Unable to determine 'assign' in model.matrix.mppm", call.=FALSE)
    }
  }

  ## split if required
  if(separate) {
    id <- df$id
    mm <- split.data.frame(mm, id)  # see help(split)
    if(!is.null(ass))
      mm <- lapply(mm, "attr<-", which="assign", value=ass)
    if(!is.null(ctr))
      mm <- lapply(mm, "attr<-", which="contrasts", value=ctr)
  }

  return(mm)
}

