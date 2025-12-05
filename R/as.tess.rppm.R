#' as.tess.rppm.R
#'
#' Create tessellation from rppm object
#'
#' $Revision: 1.6 $ $Date: 2025/12/05 06:12:30 $

as.tess.rppm <- function(X) {
  verifyclass(X, "rppm")
  parseRPPM(X, weed=TRUE, result="tess")
}

parseRPPM <- function(X, ..., weed=TRUE,
                      result=c("tess", "expression", "plain"),
                      verbose=FALSE) {
  result <- match.arg(result)
  ## extract ppm fit
  pfit <- X$pfit
  W    <- Window(pfit)
  ## extract the recursive partition
  rp <- X$rp
  ## get name of variable used in each split ("<leaf>" = no split)
  frame <- rp$frame
  nnodes <- nrow(frame)
  varname <- frame$var
  is.leaf <- (varname == "<leaf>")
  yhat    <- frame$yval
  ## detect constant fit
  if(all(is.leaf)) {
    z <- switch(result,
           expression = expression(),
           plain      = "",
           tess       = setmarks(as.tess(W), yhat[1L]))
    return(z)
  }
  ## extract the original covariates (spatial objects)
  spatialcovariates <- pfit$covariates
  ## identify whether each variable is numeric, logical or categorical
  ##     Note: attr(rp, "xlevels") would identify factors & give their levels,
  ##           but would not identify logical covariates
  gd <- getglmdata(pfit)
  isfac <- sapply(gd, is.factor)
  islog <- sapply(gd, is.logical)
  gdtype <- ifelse(isfac, "factor", ifelse(islog, "logical", "numeric"))
  faclev <- lapply(gd[isfac], levels)
  ## check all named covariates in the rpart object are available
  needed <- unique(varname[!is.leaf])
  if("marks" %in% needed)
    stop("Not yet implemented for models that involve marks", call.=FALSE)
  if(!all(found <- needed %in% c(names(spatialcovariates), "x", "y"))) {
    missed <- needed[!found]
    nmissed <- length(missed)
    stop(paste(ngettext(nmissed, "Covariate", "Covariates"),
               commasep(sQuote(missed)),
               ngettext(nmissed, "was", "were"),
               "not found"),
         call.=FALSE)
  }
  ## Labels for left and right split at each node (from rpart package)
  ##     e.g. ">= 1.2", "< 1.2" for numerical variable
  ##          "ac", "bde" for categorical variable with 5 levels
  labe <- labels(rp, collapse=FALSE)
  ## Each split is a division of space into 2 tiles.
  ## Create these tessellations and text labels
  splitTess <- rep(list(NAobject("tess")), nnodes)
  splitLabels <- unname(labe)
  splitValues <- vector(mode="list", length=nnodes)
  for(irow in which(!is.leaf)) {
    var.i <- varname[irow]
    lab.i <- labe[irow, ]
    if(var.i %in% c("x", "y")) {
      ## split on cartesian coordinate
      splitLabels[irow, ] <- splab <- paste(var.i, lab.i)
      critval <- gsub(">", "", gsub("<", "", gsub("=", "", lab.i[1L])))
      critval <- as.numeric(critval)
      whichbelow <- grep("<", lab.i)
      R <- Frame(W)
      til <- switch(var.i,
             x = list(owin(c(R$xrange[1], critval), R$yrange),
                      owin(c(critval, R$xrange[2]), R$yrange)),
             y = list(owin(R$xrange, c(R$yrange[1], critval)),
                      owin(R$xrange, c(critval, R$yrange[2])))
             )
      if(whichbelow == 2)
        til <- rev(til)
      names(til) <- splab
      splitTess[[irow]] <- tess(tiles=til, window=W)
      splitValues[[irow]] <-
        if(whichbelow == 1) list(c(-Inf, critval), c(critval, Inf)) else
                            list(c(critval, Inf), c(-Inf, critval))
    } else {
      ## split on a supplied covariate
      cov.i <- spatialcovariates[[var.i]]
      switch(gdtype[[var.i]],
             numeric = {
               ## ---------------- NUMERIC VALUED ------------------
               splitLabels[irow, ] <- splab <- paste(var.i, lab.i)
               critval <- gsub(">", "", gsub("<", "", gsub("=", "", lab.i[1L])))
               critval <- as.numeric(critval)
               whichbelow <- grep("<", lab.i)
               splitValues[[irow]] <- if(whichbelow == 1) {
                                        list(c(-Inf, critval), c(critval, Inf))
                                      } else {
                                        list(c(critval, Inf), c(-Inf, critval))
                                      }
               if(is.function(cov.i) &&
                  !inherits(cov.i, c("distfun", "tessfun"))) {
                    ## convert function(x,y) to pixel image
                 cov.i <- as.im(cov.i, W=W)
               }
               if(is.im(cov.i)) {
                 isbelow <- eval(parse(text=paste("eval.im(factor(cov.i",
                                                  lab.i[1L],
                                                  ", levels=c(TRUE, FALSE)))")))
                 levels(isbelow) <- splab
                 splitTess[[irow]] <- tess(image=isbelow)
               } else if(inherits(cov.i, "distfun")) {
                 env.i <- environment(cov.i)
                 X <- get("X", envir=env.i)
                 if(is.owin(X) && isTRUE(get("invert", envir=env.i))) {
                   ## distance to complement of a window
                   above <- erosion(X, critval)
                   below <- setminus.owin(W, above)
                 } else {
                   ## distance to some object
                   below <- dilation(X, critval)
                   above <- setminus.owin(W, below)
                 }
                 til <- list(below, above)
                 if(whichbelow == 2) til <- rev(til)
                 names(til) <- splab
                 splitTess[[irow]] <- tess(tiles=til, window=W)
               } else if(inherits(cov.i, "tessfun")) {
                 env.i <- environment(cov.i)
                 Tess <- get("V", envir=env.i)
                 W <- Window(Tess)
                 funvalues <- tessfunvalues(cov.i)
                 leftsplit <- eval(parse(text=paste("funvalues", lab.i[1L])))
                 if(all(leftsplit)) {
                   leftset <- W
                   rightset <- emptywindow(Frame(W))
                 } else if(!any(leftsplit)) {
                   leftset <- emptywindow(Frame(W))
                   rightset <- W
                 } else {
                   til <- tiles(Tess)
                   leftset <- do.call(union.owin, til[leftsplit])
                   rightset <- setminus.owin(W, leftset)
                 }
                 til <- list(leftset, rightset)
                 if(whichbelow == 2) til <- rev(til)
                 names(til) <- splab
                 splitTess[[irow]] <- tess(tiles=til, window=W)
               } else if(is.owin(cov.i)) {
                 inside <- cov.i
                 outside <- setminus.owin(W, cov.i)
                 til <- list(inside, outside)
                 splab <- paste0(var.i, "==", c("TRUE", "FALSE"))
                 if(whichbelow == 2) {
                   til <- rev(til)
                   splab <- rev(splab)
                 }
                 splitLabels[irow, ] <- splab
                 names(til) <- splab
                 splitTess[[irow]] <- tess(tiles=til, window=W)
                 splitValues[[irow]] <- if(whichbelow == 1) {
                                        list(FALSE, TRUE)
                                      } else {
                                        list(TRUE, FALSE)
                                      }
               } else {
                 stop(paste("Unable to interpret covariate",
                            sQuote(var.i),
                            "(expecting numeric image)"),
                      call.=FALSE)
               }
             },
             factor = {
               ## ---------------- FACTOR VALUED ------------------
               lev.i <- levels(cov.i)
               ## lab.i contains levels mapped to a, b, c,...
               levcode <- strsplit(unname(lab.i), "")
               levindex <- lapply(levcode, match, table=letters)
               levactual <- lapply(levindex, function(j) { lev.i[j] })
               splitValues[[irow]] <- levactual
               nlev <- lengths(levindex)
               levstring <- sapply(levactual, paste, collapse=",")
               levQstring <- sapply(levactual,
                         function(x) paste(dQuote(x, FALSE), collapse=","))
               levstring <- ifelse(nlev == 1,
                                   paste("==", levstring),
                                   paste0("%in% c(", levstring, ")"))
               levQstring <- ifelse(nlev == 1,
                                    paste("==", levQstring),
                                    paste0("%in% c(", levQstring, ")"))
               splab <- paste(var.i, levstring)
               if(is.function(cov.i) && !inherits(cov.i, "tessfun")) 
                 cov.i <- as.im(cov.i, W=W)
               if(inherits(cov.i, "tessfun") && as.tess(cov.i)$type == "image")
                 cov.i <- as.im(cov.i, W=W)
               if(is.im(cov.i)) {
                 ## factor valued image.
                 ## Make an image tessellation
                 isleft <- eval(parse(text=paste("eval.im(factor(cov.i",
                                                 levQstring[1L],
                                                 ", levels=c(TRUE, FALSE)))")))
                 levels(isleft) <- splab
                 splitTess[[irow]] <- tess(image=isleft)
               } else if(is.tess(cov.i)) {
                 til <- tiles(cov.i)
                 if(nlev[1L] == 1) {
                   leftlev <- levactual[[1L]]
                   left <- til[[leftlev]]
                   right <- setminus.owin(W, left)
                 } else if(nlev[2L] == 1) {
                   rightlev <- levactual[[2L]]
                   right <- til[[rightlev]]
                   left <- setminus.owin(W, right)
                 } else {
                   leftlevs <- levactual[[1L]]
                   rightlevs <- levactual[[2L]]
                   left <- do.call(union.owin, til[leftlevs])
                   right <- do.call(union.owin, til[rightlevs])
                 }
                 splittiles <- list(left, right)
                 names(splittiles) <- splab
                 splitTess[[irow]] <- tess(tiles=splittiles, window=W)
               } else if(inherits(cov.i, "tessfun")) {
                 ## original covariate object was a logical valued function
                 ## that is constant on each tlle of a tessellation
                 Tess <- get("V", envir=environment(cov.i))
                 val <- tessfunvalues(cov.i)
                 til <- tiles(Tess)
                 if(nlev[1L] == 1) {
                   leftlev <- levactual[[1L]]
                   left <- til[[leftlev]]
                   right <- setminus.owin(W, left)
                 } else if(nlev[2L] == 1) {
                   rightlev <- levactual[[2L]]
                   right <- til[[rightlev]]
                   left <- setminus.owin(W, right)
                 } else {
                   leftlevs <- levactual[[1L]]
                   rightlevs <- levactual[[2L]]
                   left <- do.call(union.owin, til[leftlevs])
                   right <- do.call(union.owin, til[rightlevs])
                 }
                 splittiles <- list(left, right)
                 names(splittiles) <- splab
                 splitTess[[irow]] <- tess(tiles=splittiles, window=W)
               } else {
                 stop(paste("Unable to interpret covariate", sQuote(var.i),
                            "(expecting factor)"),
                      call.=FALSE)
               }
               splitLabels[irow, ] <- splab
             },
             logical = {
               ## ---------------- LOGICAL VALUED ------------------
               whichbelow <- grep("<", lab.i)
               logique <- c("FALSE", "TRUE")
               lab.i <- if(whichbelow == 1) logique else rev(logique)
               splitLabels[irow, ] <- splab <- paste(var.i, "==", lab.i)
               splitValues[[irow]] <- if(whichbelow == 1) {
                                        list(FALSE, TRUE)
                                      } else {
                                        list(TRUE, FALSE)
                                      }
               if(is.function(cov.i) && !inherits(cov.i, "tessfun")) 
                 cov.i <- as.im(cov.i, W=W)
               if(inherits(cov.i, "tessfun") && as.tess(cov.i)$type == "image")
                 cov.i <- as.im(cov.i, W=W)
               if(is.im(cov.i)) {
                 istrue <- eval.im(factor(cov.i, levels=c(TRUE, FALSE)))
                 levels(istrue) <- splab
                 splitTess[[irow]] <- tess(image=istrue)
               } else if(is.owin(cov.i)) {
                 inside <- cov.i
                 outside <- setminus.owin(W, cov.i)
                 til <- list(inside, outside)
                 splab <- paste0(var.i, "==", c("TRUE", "FALSE"))
                 if(whichbelow == 2) {
                   til <- rev(til)
                   splab <- rev(splab)
                 }
                 splitLabels[irow, ] <- splab
                 names(til) <- splab
                 splitTess[[irow]] <- tess(tiles=til, window=W)
               } else if(inherits(cov.i, "tessfun")) {
                 ## original covariate object was a logical valued function
                 ## that is constant on each tlle of a tessellation
                 Tess <- get("V", envir=environment(cov.i))
                 truetiles <- as.logical(tessfunvalues(cov.i))
                 if(!any(truetiles)) {
                   trueset <- emptywindow(W)
                   falseset <- W
                 } else if(all(truetiles)) {
                   trueset <- W
                   falseset <- emptywindow(W)
                 } else {
                   trueset <- do.call(union.owin, tiles(Tess)[truetiles])
                   falseset <- setminus.owin(W, trueset)
                 }
                 splittiles <- list(falseset, trueset)
                 if(whichbelow == 2) splittiles <- rev(splittiles)
                 names(splittiles) <- splab
                 splitTess[[irow]] <- tess(tiles=splittiles, window=W)
               } else {
                 stop(paste("Unable to interpret covariate",
                            sQuote(var.i),
                            "(expecting logical image)"),
                      call.=FALSE)
               }
             })
    }
  }
  ## determine nesting depth of each node
  nodecode <- as.integer(row.names(frame))
  depth <- floor(log2(nodecode) + 1e-07)
  depth <- depth - min(depth)
  ## traverse tree depth-first, collecting conditions for each leaf
  leaves <- list()
  condition <- character(0)
  conditionsets <- leafsets <- list()
  conditionvariable <- character(0)
  conditionvalues <- list()
  leafyhat <- numeric(0)
  if(verbose)
    print(splitLabels)
  for(irow in seq_len(nnodes)) {
    newdepth <- depth[irow]
    prevdepth <- if(irow == 1L) -1L else depth[irow - 1L]
    if(verbose) {
      cat(paste("row", irow, "depth", newdepth, "previous", prevdepth,
                "\n condition:"))
      print(condition)
    }
    if(newdepth > prevdepth) {
      ## current node is left branch of previous node
      if(is.leaf[irow]) {
        if(verbose) cat("save condition as a leaf\n")
        ## save this condition as one of the leaves
        if(weed) {
          if(verbose) cat("weeding...\n")
          ok <- !redundantConditions(conditionvariable, conditionvalues,
                                     verbose=verbose)
        } else {
          ok <- rep(TRUE, length(conditionvariable))
        }
        okcond <- condition[ok]
        oksets <- conditionsets[ok]
        leaves <- append(leaves, list(okcond))
        if(length(okcond) == 1) {
          newleaf <- oksets
        } else {
          newleaf <- list(do.call(intersect.owin, oksets))
          names(newleaf) <- paste(paren(names(oksets)), collapse = "&")
        }
        leafsets <- append(leafsets, newleaf)
        leafyhat <- c(leafyhat, yhat[irow])
      } else {
        ## continue descending, following another left branch
        condition <- c(condition, splitLabels[irow, 1L])
        conditionvariable <- c(conditionvariable, varname[irow])
        conditionvalues <- append(conditionvalues,
                                  list(splitValues[[irow]][[1L]]))
        conditionsets <- append(conditionsets,
                                tiles(splitTess[[irow]])[1L])
        if(verbose) {
          cat("appending, new condition:\n")
          print(condition)
        }
      }
    } else {
      ## current node is the right branch of a previous node
      ## cut condition back to before this depth
      parentdepth <- newdepth - 1L
      backto <- seq_len(parentdepth)
      condition <- condition[backto]
      conditionsets <- conditionsets[backto]
      conditionvariable <- conditionvariable[backto]
      conditionvalues <- conditionvalues[backto]
      if(verbose) {
        cat(paste("Scraping back to depth", parentdepth, "\n"))
        cat("Truncated condition:\n")
        print(condition)
      }
      ## find which node
      prevrows <- seq_len(irow - 1L)
      last <- max(prevrows[depth[prevrows] == newdepth-1L])
      if(verbose) {
        cat(paste("Referring back to previous node", last,
                  ", depth", depth[last], "\n",
                  "which is a split on", varname[last], "\n Labels:"))
        print(splitLabels[last,])
      }
      ## use the right branch of this node
      condition <- c(condition, splitLabels[last, 2L])
      conditionvariable <- c(conditionvariable, varname[last])
      conditionvalues <- append(conditionvalues,
                                list(splitValues[[last]][[2L]]))
      conditionsets <- append(conditionsets,
                              tiles(splitTess[[last]])[2L])
      if(verbose) {
        cat("new condition:\n")
        print(condition)
        print(conditionvariable)
        
      }
      if(is.leaf[irow]) {
        ## save this condition as one of the leaves
        if(verbose) cat("save condition as a leaf\n")
        ## save this condition as one of the leaves
        if(weed) {
          ok <- !redundantConditions(conditionvariable, conditionvalues,
                                     verbose=verbose)
        } else {
          ok <- rep(TRUE, length(conditionvariable))
        }
        okcond <- condition[ok]
        oksets <- conditionsets[ok]
        leaves <- append(leaves, list(okcond))
        if(length(okcond) == 1) {
          newleaf <- oksets
        } else {
          newleaf <- list(do.call(intersect.owin, oksets))
          names(newleaf) <- paste(paren(names(oksets)), collapse = "&")
        }
        leafsets <- append(leafsets, newleaf)
        leafyhat <- c(leafyhat, yhat[irow])
      } else {
        ## continue descending, following another left branch
        condition <- c(condition, splitLabels[irow, 1L])
        conditionvariable <- c(conditionvariable, varname[irow])
        conditionvalues <- append(conditionvalues,
                                  list(splitValues[[irow]][[1L]]))
        conditionsets <- append(conditionsets,
                                tiles(splitTess[[irow]])[1L])
        if(verbose) {
          cat("appending, new condition:\n")
          print(condition)
        }
      }
    }
  }
  z <- switch(result,
              plain      = leaves,
              expression = lapply(leaves, str2expression),
              tess       = tess(tiles=leafsets, window=W, marks=leafyhat))
  return(z)
}
  
  
redundantConditions <- function(splitvariable, splitvalues,
                                verbose=TRUE) {
  nc <- length(splitvariable)
  if(verbose) {
    cat(paste("Entering redundantConditions with", nc, "conditions and",
              length(splitvalues), "sets of split values\n"))
    print(splitvariable)
  }
  redundant <- logical(nc) # all FALSE
  if(!anyDuplicated(splitvariable))
    return(redundant)
  for(j in seq_len(nc)) {
    varj <- splitvariable[j]
    valj <- splitvalues[[j]]
    others <- setdiff(which(splitvariable == varj), j)
    if(length(others)) {
      ## compare with other cs for the same variable
      if(is.numeric(valj)) {
        ## numerical variable: half-infinite intervals
        if(valj[1L] == -Inf) {
          ## left-infinite interval
          crit <- valj[2L]
          ## redundant if a more restrictive interval exists
          v1 <- sapply(splitvalues[others], "[", i=1L)
          v2 <- sapply(splitvalues[others], "[", i=2L)
          redundant[j] <- any(v1 == -Inf & v2 < crit)
        } else {
          ## right-infinite interval
          crit <- valj[1L]
          ## redundant if a more restrictive interval exists
          v1 <- sapply(splitvalues[others], "[", i=1L)
          v2 <- sapply(splitvalues[others], "[", i=2L)
          redundant[j] <- any(v2 == Inf & v1 > crit)
        }
      } else if(is.character(valj)) {
        ## categorical variable
        ## redundant if a more restrictive set of levels is given
        redundant[j] <- any(sapply(splitvalues[others],
                                   function(x) all(x %in% valj)))
      }
    }
  }
  if(verbose) {
    if(!any(redundant)) {
      cat("no redundant conditions\n")
    } else {
      nr <- sum(redundant)
      splat(paste(ngettext(nr, "condition", "conditions"),
                  commasep(which(redundant)),
                  ngettext(nr, "is", "are"), "redundant"))
    }
  }
  return(redundant)
}
