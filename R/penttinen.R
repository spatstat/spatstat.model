#
#
#    penttinen.R
#
#    $Revision: 1.4 $	$Date: 2026/01/08 10:04:34 $
#
#    Penttinen pairwise interaction
#
#
# -------------------------------------------------------------------
#	

Penttinen <- local({

  # .......... auxiliary functions ................
  pentTerms <- function(X, Y, idX, idY, radius) {
    stopifnot(is.numeric(radius))
    # sort in increasing order of x coordinate
    oX <- fave.order(X$x)
    oY <- fave.order(Y$x)
    Xsort <- X[oX]
    Ysort <- Y[oY]
    idXsort <- idX[oX]
    idYsort <- idY[oY]
    nX <- npoints(X)
    nY <- npoints(Y)
    # call C routine
    out <- .C(SM_Epent,
            nnsource = as.integer(nX),
            xsource  = as.double(Xsort$x),
            ysource  = as.double(Xsort$y),
            idsource = as.integer(idXsort),
            nntarget = as.integer(nY),
            xtarget  = as.double(Ysort$x),
            ytarget  = as.double(Ysort$y),
            idtarget = as.integer(idYsort),
            radius   = as.double(radius),
            values   = as.double(double(nX)),
            PACKAGE="spatstat.model")
    answer <- numeric(nX)
    answer[oX] <- out$values
    return(answer)
  }
  
  # create blank template object without family and pars

  BlankAntti <-
  list(
       name     = "Penttinen process",
       creator  = "Penttinen",
       family    = "pairwise.family", # evaluated later
       pot      = function(d, par) {
         ans <- numeric(length(d))
         dim(ans) <- dim(d)
         zz <- d/(2 * par$r)
         ok <- (zz < 1)
         z <- zz[ok]
         ans[ok] <- (2/pi) * (acos(z) - z * sqrt(1-z^2))
         return(ans)
       },
       par      = list(r = NULL), # to be filled in
       parnames = "circle radius",
       hasInf = FALSE,
       init     = function(self) {
         r <- self$par$r
         if(!is.numeric(r) || length(r) != 1 || r <= 0)
           stop("interaction distance r must be a positive number")
       },
       update = NULL,  # default OK
       print = NULL,    # default OK
       interpret =  function(coeffs, self) {
         theta <- as.numeric(coeffs[1])
         gamma <- exp(theta)
         return(list(param=list(gamma=gamma),
                     inames="interaction parameter gamma",
                     printable=dround(gamma)))
       },
       valid = function(coeffs, self) {
         theta <- as.numeric(coeffs[1])
         return(is.finite(theta) && (theta <= 0))
       },
       project = function(coeffs, self) {
         if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$r
         if(anyNA(coeffs))
           return(2 * r)
         theta <- coeffs[1]
         if(abs(theta) <= epsilon)
           return(0)
         else
           return(2 * r)
       },
       version=NULL, # to be filled in 
       ## fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                         splitInf=FALSE, ...) {
         dont.complain.about(splitInf)
         ## fast evaluator for Penttinen interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Penttinen")
         radius <- potpars$r
         idX <- seq_len(npoints(X))
         idU <- rep.int(-1L, npoints(U))
         idU[EqualPairs[,2L]] <- EqualPairs[,1L]
         values <- pentTerms(U, X, idU, idX, radius)
         result <- matrix(values, ncol=1L)
         return(result)
       }
    )
  class(BlankAntti) <- "interact"


  # Finally define main function
  
  Penttinen <- function(r) {
    instantiate.interact(BlankAntti, list(r=r))
  }

  Penttinen <- intermaker(Penttinen, BlankAntti)
  
  Penttinen
})

