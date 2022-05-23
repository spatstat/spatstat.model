#'
#'     newformula.R
#'
#'    $Revision: 1.4 $ $Date: 2022/01/19 08:50:37 $
#' 
#'   Update formula and expand polynomial

newformula <- function(old, change, eold, enew,
                       expandpoly=spatstat.options("expand.polynom"),
                       dotvars=character(0)) {
  old <- if(is.null(old)) ~1 else eval(old, eold)
  change <- if(is.null(change)) ~1 else eval(change, enew)
  old <- as.formula(old, env=eold)
  change <- as.formula(change, env=enew)
  if(expandpoly) {
    old <- expand.polynom(old)
    change <- expand.polynom(change)
  }
  old <- expandDot(old, dotvars)
  answer <- update.formula(old, change)
  return(answer)
}

expandDot <- local({

  hasDot <- function(x) { "." %in% all.names(x) }

  expandDot <- function(f, dotvars) {
    if(length(dotvars) == 0) return(f)
    dotsum <- paren(paste(dotvars, collapse=" + "))
    dotexpr <- rhs.of.formula(as.formula(paste("~", dotsum)), tilde=FALSE)
    g <- fuddle(f, dotexpr)
    environment(g) <- environment(f)
    return(g)
  }

  fuddle <- function(f, dotexpr) {
    print(f)
    if(!hasDot(f)) return(f)
    if(identical(f, as.name('.')))
      return(dotexpr)
    if(length(f) == 1) return(f)
    if(identical(f[[1]], as.name('I'))) {
      ## expressions enclosed in I() are protected
      return(f)
    } 
    tbd <- unlist(lapply(f, hasDot))
    if(any(tbd)) {
      ## descend recursively
      for(i in which(tbd)) 
        f[[i]] <- fuddle(f[[i]], dotexpr)
    }
    return(f)
  }

  expandDot
})

