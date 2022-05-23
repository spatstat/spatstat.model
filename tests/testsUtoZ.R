#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.model
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.model)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#
#  tests/undoc.R
#
#   $Revision: 1.16 $   $Date: 2020/11/02 07:06:49 $
#
#  Test undocumented hacks, experimental code, etc


local({
  if(FULLTEST) {
    ## cases of 'pickoption'
    aliases <- c(Lenin="Ulyanov", Stalin="Djugashvili", Trotsky="Bronstein")
    surname <- "Trot"
    pickoption("leader", surname,  aliases)
    pickoption("leader", surname,  aliases, exact=TRUE, die=FALSE)
  }
  if(ALWAYS) {
    ## pixellate.ppp accepts a data frame of weights
    pixellate(cells, weights=data.frame(a=1:42, b=42:1))
    ## test parts of 'rmhsnoop' that don't require interaction with user
    rmhSnoopEnv(cells, Window(cells), 0.1)
  }
  if(FULLTEST) {
    ## Berman-Turner frame
    A <- bt.frame(quadscheme(cells), ~x, Strauss(0.07), rbord=0.07)
    print(A)
    ## digestCovariates
    D <- distfun(cells)
    Z <- distmap(cells)
    U <- dirichlet(cells)
    stopifnot(is.scov(D))
    stopifnot(is.scov(Z))
    stopifnot(is.scov(U))
    stopifnot(is.scov("x"))
    dg <- digestCovariates(D=D,Z=Z,U=U,"x",list(A="x", B=D))
    ##
    a <- getfields(dg, c("A", "D", "niets"), fatal=FALSE)
    ## util.R
    gg <- pointgrid(owin(), 7)
    checkbigmatrix(1000000L, 1000000L, FALSE, TRUE)
    spatstatDiagnostic("whatever")
    M <- list(list(a=2, b=FALSE),
              list(a=2, b=TRUE))
    stopifnot(!allElementsIdentical(M))
    stopifnot(allElementsIdentical(M, "a"))
    ##
    A <- Strauss(0.1)
    A <- reincarnate.interact(A)
    ##
    ## special lists
    B <- solist(a=cells, b=redwood, c=japanesepines)
    BB <- as.ppplist(B)
    BL <- as.layered(B)
    DB <- as.imlist(lapply(B, density))
    is.solist(B)
    is.ppplist(B)
    is.imlist(DB)
    ## case of density.ppplist 
    DEB <- density(BB, se=TRUE)
  }

  if(ALWAYS) {
    ## fft
    z <- matrix(1:16, 4, 4)
    a <- fft2D(z, west=FALSE)
    if(fftwAvailable())
      b <- fft2D(z, west=TRUE)
  }

  if(ALWAYS) {
    ## experimental interactions
    pot <- function(d, par) { d <= 0.1 }
    A <- Saturated(pot)
    print(A)
    A <- update(A, name="something")
    ppm(amacrine ~ x, A, rbord=0.1)
  }

  if(ALWAYS) { # platform dependent
    #' version-checking
    now <- Sys.Date()
    versioncurrency.spatstat(now + 80, FALSE)
    versioncurrency.spatstat(now + 140, FALSE)
    versioncurrency.spatstat(now + 400, FALSE)
    versioncurrency.spatstat(now + 1000)
  }

  if(FULLTEST) {
    #' general Ord interaction
    gradual <- function(d, pars) {
      y <- pmax(0, 0.005 - d)/0.005
      if(is.matrix(d)) y <- matrix(y, nrow(d), ncol(d))
      return(y)
    }
    B <- Ord(gradual, "gradual Ord process")
  }
})
  

##
##  tests/updateppm.R
##
##  Check validity of update.ppm
##
##  $Revision: 1.7 $ $Date: 2020/11/02 07:07:42 $

local({
  
})
#
#  tests/vcovppm.R
#
#  Check validity of vcov.ppm algorithms
#
#  Thanks to Ege Rubak
#
#  $Revision: 1.12 $  $Date: 2020/05/02 01:32:58 $
#

local({

  set.seed(42)
  X <- rStrauss(200, .5, .05)
  model <- ppm(X, inter = Strauss(.05))

  if(ALWAYS) {
    b  <- vcov(model, generic = TRUE, algorithm = "basic")
    v  <- vcov(model, generic = TRUE, algorithm = "vector")
    vc <- vcov(model, generic = TRUE, algorithm = "vectorclip")
    vn <- vcov(model, generic = FALSE)

    disagree <- function(x, y, tol=1e-7) { max(abs(x-y)) > tol }
    asymmetric <- function(x) { disagree(x, t(x)) }

    if(asymmetric(b))
      stop("Non-symmetric matrix produced by vcov.ppm 'basic' algorithm")
    if(asymmetric(v))
      stop("Non-symmetric matrix produced by vcov.ppm 'vector' algorithm")
    if(asymmetric(vc))
      stop("Non-symmetric matrix produced by vcov.ppm 'vectorclip' algorithm")
    if(asymmetric(vn))
      stop("Non-symmetric matrix produced by vcov.ppm Strauss algorithm")
    
    if(disagree(v, b))
      stop("Disagreement between vcov.ppm algorithms 'vector' and 'basic' ")
    if(disagree(v, vc))
      stop("Disagreement between vcov.ppm algorithms 'vector' and 'vectorclip' ")
    if(disagree(vn, vc))
      stop("Disagreement between vcov.ppm generic and Strauss algorithms")
  }

  if(ALWAYS) { # C code
    ## Geyer code
    xx <- c(0.7375956, 0.6851697, 0.6399788, 0.6188382)
    yy <- c(0.5816040, 0.6456319, 0.5150633, 0.6191592)
    Y <- ppp(xx, yy, window=square(1))
    modelY <- ppm(Y ~1, Geyer(0.1, 1))
    
    b  <- vcov(modelY, generic = TRUE, algorithm = "basic")
    v  <- vcov(modelY, generic = TRUE, algorithm = "vector")
    vc <- vcov(modelY, generic = TRUE, algorithm = "vectorclip")
    
    if(asymmetric(b))
      stop("Non-symmetric matrix produced by vcov.ppm 'basic' algorithm for Geyer model")
    if(asymmetric(v))
      stop("Non-symmetric matrix produced by vcov.ppm 'vector' algorithm for Geyer model")
    if(asymmetric(vc))
      stop("Non-symmetric matrix produced by vcov.ppm 'vectorclip' algorithm for Geyer model")
  
    if(disagree(v, b))
      stop("Disagreement between vcov.ppm algorithms 'vector' and 'basic' for Geyer model")
    if(disagree(v, vc))
      stop("Disagreement between vcov.ppm algorithms 'vector' and 'vectorclip' for Geyer model")
  }

  if(ALWAYS) { # C code
    ## tests of 'deltasuffstat' code
    ##     Handling of offset terms
    modelH <- ppm(cells ~x, Hardcore(0.05))
    a <- vcov(modelH, generic=TRUE) ## may fall over
    b <- vcov(modelH, generic=FALSE)
    if(disagree(a, b))
      stop("Disagreement between vcov.ppm algorithms for Hardcore model")
  
    ##     Correctness of pairwise.family$delta2
    modelZ <- ppm(amacrine ~1, MultiStrauss(radii=matrix(0.1, 2, 2)))
    b <- vcov(modelZ, generic=FALSE)
    g <- vcov(modelZ, generic=TRUE)
    if(disagree(b, g))
      stop("Disagreement between vcov.ppm algorithms for MultiStrauss model")

    ## Test that 'deltasuffstat' works for Hybrids
    modelHyb <- ppm(japanesepines ~ 1, Hybrid(Strauss(0.05), Strauss(0.1)))
    vHyb <- vcov(modelHyb)
  }

  if(FULLTEST) {
    ## Code blocks for other choices of 'what'
    model <- ppm(X ~ 1, Strauss(.05))
    cG <- vcov(model, what="corr")
    cP <- vcov(update(model, Poisson()), what="corr")
    ## outdated usage
    cX <- vcov(model, A1dummy=TRUE)

    ## Model with zero-length coefficient vector
    lam <- intensity(X)
    f <- function(x,y) { rep(lam, length(x)) }
    model0 <- ppm(X ~ offset(log(f)) - 1)
    dd <- vcov(model0)
    cc <- vcov(model0, what="corr")
  
    ## Model with NA coefficients
    fit <- ppm(X ~ log(f))
    vcov(fit)
    fitE <- emend(fit, trace=TRUE)
  
    ## Other weird stuff
    su <- suffloc(ppm(X ~ x))
  }
})
