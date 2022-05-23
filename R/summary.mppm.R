#
# summary.mppm.R
#
# $Revision: 1.19 $  $Date: 2022/01/04 05:30:06 $
#


summary.mppm <- function(object, ..., brief=FALSE) {
  # y will be the summary 
  y <- object[c("Call", "Info", "Inter", "trend", "iformula",
                "random",
                "npat", "maxlogpl")]
  y$brief <- brief

  Info  <- object$Info
  Inter <- object$Inter
  FIT   <- object$Fit$FIT
  moadf <- object$Fit$moadf

  y$Fit <- object$Fit[c("fitter", "use.gam", "fmla",
                        "Vnamelist", "Isoffsetlist")]
  y$Fit$FIT <- summary(FIT)
  y$Fit$moadf <- list(nrow=nrow(moadf), colnames=colnames(moadf))
  
  ninteract    <- Inter$ninteract
  interaction  <- Inter$interaction
  iused        <- Inter$iused
  itags        <- Inter$itags
  processnames <- Inter$processes
  constant     <- Inter$constant
  trivial      <- Inter$trivial

  npat      <- y$npat
  iformula  <- y$iformula
  random    <- y$random
  Vnamelist <- y$Fit$Vnamelist
  allVnames <- unlist(Vnamelist)
  Isoffsetlist <- y$Fit$Isoffsetlist
  poistags  <- itags[trivial]
  intertags <- c(allVnames, poistags)

  ## does the model depend on covariates?
  y$depends.covar <- Info$has.covar && (length(Info$used.cov.names) > 0)

#  rownames  <- y$Info$rownames
  
  switch(y$Fit$fitter,
         glmmPQL={
           y$coef <- co <- fixed.effects(FIT)
           y$coef.rand <- random.effects(FIT)
         },
         gam=,
         glm={
           y$coef <- co <- coef(FIT)
         })

  ## identify model terms which involve interpoint interaction
  md <- model.depends(FIT)
  is.interpoint <- colnames(md) %in% intertags
  involves.interpoint <- apply(md[ , is.interpoint, drop=FALSE], 1, any)
  y$coef.inter <- co[involves.interpoint]
  ## identify trend and design coefficients 
  systematic <- !involves.interpoint
  y$coef.syst <- co[systematic]

  # random effects
  y$ranef <- if(Info$has.random) summary(FIT$modelStruct) else NULL

  ### Interpoint interactions 
  # model is Poisson ?
  y$poisson <- ispois <- all(trivial[iused])
  # Determine how complicated the interactions are:
  # (0) are there random effects involving the interactions
  randominteractions <-
    !is.null(random) && any(variablesinformula(random) %in% itags)
  # (1) is the interaction formula of the form ~ tag + tag + ... + tag
  isimple  <- identical(sort(variablesinformula(iformula)),
                        sort(termsinformula(iformula)))
  # (2) is it of the form ~tag 
  trivialformula <- (isimple && ninteract == 1)
  # (3) is it of the form ~tag where the interaction is the same in each row
  fixedinteraction <- (trivialformula && constant && !randominteractions)
  
  ### Determine printing of interactions, accordingly ###
  iprint <- list()
  if(randominteractions) {
    toohard <- TRUE
    printeachrow <- FALSE
  } else 
  if(fixedinteraction || ispois) {    
    # exactly the same interaction for all patterns
    interaction <- interaction[1L,1L,drop=TRUE]
    fi.all <- fii(interaction, co, Vnamelist[[1L]], Isoffsetlist[[1L]]) 
    iprint <- list("Interaction for all patterns"=fi.all)
    printeachrow <- FALSE
    toohard      <- FALSE
  } else if(trivialformula) {
    ## same interaction structure for all patterns;
    ## irregular parameters may be different on each row;
    ## regular parameters of interaction do not depend on design
    pname <-  unlist(processnames)[iused]
    iprint <- list("Interaction for each pattern" = pname)
    printeachrow <- TRUE
    toohard      <- FALSE
  } else if(sum(iused) == 1) {
    ## same interaction structure for all patterns;
    ## irregular parameters may be different on each row;
    ## regular parameters of interaction may depend on design
    pname <-  unlist(processnames)[iused]
    iprint <- list("Interaction for each pattern" = pname)
    printeachrow <- TRUE
    toohard      <- FALSE
    ## look for design : interaction terms
    mm <- md[involves.interpoint, !is.interpoint, drop=FALSE]
    tangled <- apply(mm, 2, any)
    if(any(tangled)) {
      tanglednames <- colnames(mm)[tangled]
      textra <- list(commasep(sQuote(tanglednames)))
      names(textra) <- paste("Interaction depends on design",
                             ngettext(length(tanglednames),
                                      "covariate", "covariates"))
      iprint <- append(iprint, textra)
    }
  } else if(isimple && all(constant)) {
    # several interactions involved, each of which is the same for all patterns
    iprint <- list("Interaction formula"=iformula,
                   "Interactions defined for each pattern"=NULL)
    for(j in (1:ninteract)[iused]) {
      name.j <- paste("Interaction", sQuote(itags[j]))
      int.j <- Inter$interaction[1L,j,drop=TRUE]
      Vnames.j <- Vnamelist[[j]]
      Isoffset.j <- Isoffsetlist[[j]]
      fii.j <- fii(int.j, co, Vnames.j, Isoffset.j)
      extra.j <- list(fii.j)
      names(extra.j) <- name.j
      iprint <- append(iprint, extra.j)
    }
    printeachrow <- FALSE
    toohard      <- FALSE
  } else {
    # general case
    # determine which interaction(s) are active on each row
    active <- active.interactions(object)
    if(ninteract > 1 || !all(active)) 
      iprint <- list("Active interactions"=active)
    printeachrow <- TRUE
    toohard <- any(rowSums(active) > 1)
  }

  y$ikind <- list(
                  randominteractions=randominteractions,
                  isimple           =isimple,
                  trivialformula    =trivialformula,
                  fixedinteraction  =fixedinteraction,
                  toohard           =toohard,
                  printeachrow      =printeachrow)

  
  y$depends.on.row <- ("id" %in% variablesinformula(y$trend)) || !fixedinteraction
    
  if(toohard)
    iprint <- append(iprint,
                     list("(Sorry, cannot interpret fitted interactions)"))
  else if(printeachrow) {
    subs <- subfits(object, what="interactions")
    um <- uniquemap(subs)
    uniq <- (um == seq_along(um))
    if(mean(uniq) <= 0.5) {
      icode <- cumsum(uniq)[um]
      inames <- if(max(icode) <= 26) LETTERS[icode] else as.character(icode)
      itable <- data.frame(row=seq_along(um), interaction=inames)
      usubs <- subs[um[uniq]]
      names(usubs) <- inames[uniq]
      iprint <- append(iprint,
                       list("Summary table of interactions"=itable,
                            "key to interaction table"=usubs,
                            "=== Interactions on each row ===="=NULL))
    }
    names(subs) <- paste("Interaction on row", 1:npat)
    iprint <- append(iprint, subs)
  }

  y$iprint <- iprint

  class(y) <- c("summary.mppm", class(list))
  return(y)
}


print.summary.mppm <- function(x, ..., brief=x$brief) {
  # NB: x is an object of class "summary.mppm"
  npat <- x$npat
#  Inter <- x$Inter
#  ninteract   <- Inter$ninteract
#  interaction   <- Inter$interaction
#  iused     <- Inter$iused
#  constant <- Inter$constant
#  iformula <- x$iformula
#  processnames   <- Inter$processes
#  itags   <- Inter$itags
#  trivial  <- Inter$trivial
#  random   <- x$random

  FIT <- x$Fit$FIT
#  Vnamelist <- x$Fit$Vnamelist
  
#  allVnames <- unlist(Vnamelist)
#  poistags <- itags[trivial]

  terselevel <- spatstat.options("terse")
#  rownames <- x$Info$rownames

  splat("Point process model fitted to", npat, "point patterns")
  if(waxlyrical('gory', terselevel))
    splat("Call:", x$Call$callstring)
  splat("Log trend formula:", pasteFormula(x$trend))
  switch(x$Fit$fitter,
         glmmPQL={
           cat("\nFixed effects:\n")
           print(x$coef.syst)
           cat("Random effects:\n")
           print(x$coef.rand)
           co <- fixed.effects(FIT)
         },
         gam=,
         glm={
           cat("\nFitted trend coefficients:\n")
           print(x$coef.syst)
           co <- coef(FIT)
         })

  if(length(x$coef.inter)) {
    cat("\nFitted coefficients of interpoint interaction:\n")
    print(x$coef.inter)
  }
  
  if(!brief && waxlyrical('extras', terselevel)) {
    cat("All fitted coefficients:\n")
    print(co)
  }
    
  parbreak(terselevel)

  if(!is.null(x$ranef)) {
    splat("Random effects summary:")
    print(x$ranef)
    parbreak(terselevel)
  }

  ### Print interaction information ###
  if(waxlyrical('extras', terselevel)) {
    iprint <- x$iprint 
    nama <- names(iprint) %orifnull% rep("", length(iprint))
    for(i in seq_along(iprint)) {
      nami <- nama[i]
      vali <- iprint[[i]]
      if(brief && is.matrix(vali))
        vali <- paren(paste(nrow(vali), "x", ncol(vali), "matrix"))
      if(nami != "") {
        inline <- inherits(vali, "formula") ||
                  is.character(vali) ||
                  (brief && inherits(vali, "fii"))
        if(inline) cat(paste0(nami, ":\t")) else splat(paste0(nami, ":"))
      }
      if(!is.null(vali)) {
        if(inherits(vali, "fii")) {
          print(vali, tiny=brief)
        } else if(is.character(vali)) {
          splat(vali)
        } else {
          print(vali)
        } 
      }
      parbreak(terselevel)
    }
  }

  if(!brief && waxlyrical('gory', terselevel)) {
    splat("--- Gory details: ---")
    splat("Combined data frame has", x$Fit$moadf$nrow, "rows")
    print(FIT)
  }
  invisible(NULL)
}

