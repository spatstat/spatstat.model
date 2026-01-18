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
#  tests/envelopes.R
#
#  Test validity of envelope data
#
#  $Revision: 1.29 $  $Date: 2024/01/10 13:45:29 $
#

local({
  

checktheo <- function(fit) {
  fitname <- deparse(substitute(fit))
  en <- envelope(fit, nsim=4, verbose=FALSE, nrep=1e3)
  nama <- names(en)
  expecttheo <- is.poisson(fit) && is.stationary(fit)
  context <- paste("Envelope of", fitname)
  if(expecttheo) {
    if(!("theo" %in% nama))
      stop(paste(context, "did not contain", sQuote("theo")))
    if("mmean" %in% nama)
      stop(paste(context, "unexpectedly contained", sQuote("mmean")))
  } else {
    if("theo" %in% nama)
      stop(paste(context, "unexpectedly contained", sQuote("theo")))
    if(!("mmean" %in% nama))
      stop(paste(context, "did not contain", sQuote("mmean")))
  }
  cat(paste(context, "has correct format\n"))
}

if(ALWAYS) {
  checktheo(ppm(cells ~x))
}
if(FULLTEST) {
  checktheo(ppm(cells))
  checktheo(ppm(cells ~1, Strauss(0.1)))
}


#' check savefuns/savepatterns with global 
fit <- ppm(cells~x)
if(ALWAYS) Ef <- envelope(fit, Kest, nsim=4, savefuns=TRUE, global=TRUE)
if(FULLTEST) Ep <- envelope(fit, Kest, nsim=4, savepatterns=TRUE, global=TRUE)
#' check handling of 'dangerous' cases
if(FULLTEST) {
  fut <- ppm(redwood ~ x)
  Ek <- envelope(fut, Kinhom, update=FALSE, nsim=4)
  kfut <- kppm(redwood3 ~ x)
  Ekk <- envelope(kfut, Kinhom, lambda=density(redwood3), nsim=7)
}


if(ALWAYS) { # invokes C code
  fit <- ppm(japanesepines ~ 1, Strauss(0.04))
  e6 <- envelope(fit, Kest, nsim=4, fix.n=TRUE)
  fit2 <- ppm(amacrine ~ 1, Strauss(0.03))
  e7 <- envelope(fit2, Gcross, nsim=4, fix.marks=TRUE)
}


if(FULLTEST) {
    fit <- ppm(cells ~ 1, Strauss(0.07))
    U <- envelope(fit, nsim=3, simulate=expression(runifpoint(20)))
    kfit <- kppm(redwood3 ~ x)
    UU <- envelope(kfit, nsim=7, simulate=expression(simulate(kfit, drop=TRUE)))
    VV <- envelope(kfit, nsim=7, weights=1:7)
    MM <- envelope(kfit, nsim=7, Kinhom, lambda=density(redwood3))
}

if(FULLTEST) {
  ## from Marcelino de la Cruz - scoping in update.ppm
  X <- cells
  Z <- density(X)
  pfit <- ppm(X ~ Z)
  cat("Fitted ppm(X~Z)", fill=TRUE)
  penv <- envelope(pfit, Kinhom, lambda=pfit, nsim=3)
  RX <- rotate(X, pi/3, centre="centroid")
  RZ <- density(RX)
  Rpfit <- ppm(RX ~ RZ)
  cat("Fitted ppm(RX~RZ)", fill=TRUE)
  Rpenv <- envelope(Rpfit, Kinhom, lambda=Rpfit, nsim=3)
}

if(FULLTEST) {
  #' envelope computations in other functions
  P <- lurking(cells, expression(x), envelope=TRUE, nsim=9)
  print(P)
  #' re-using envelope objects in other functions
  A <- envelope(cells, nsim=9, savepatterns=TRUE, savefuns=TRUE)
  S <- lurking(cells, expression(x), envelope=A, nsim=9)
  #' envelope.envelope
  B <- envelope(cells, nsim=5, savepatterns=TRUE, savefuns=FALSE)
  envelope(B)
}



## close 'local'
})
#'  tests/enveltest.R
#'     Envelope tests (dclf.test, mad.test)
#'     and two-stage tests (bits.test, dg.test, bits.envelope, dg.envelope)
#' 
#'     $Revision: 1.3 $  $Date: 2020/04/28 12:58:26 $ 
#'
if(FULLTEST) {
local({
  #' handling of NA function values (due to empty point patterns)
  set.seed(1234)
  X <- rThomas(5, 0.05, 10) 
  fit <- kppm(X ~ 1, "Thomas")
  set.seed(100000)
  dclf.test(fit)
  set.seed(909)
  dg.test(fit, nsim=9)
  #' other code blocks
  dclf.test(fit, rinterval=c(0, 3), nsim=9)
  envelopeTest(X, exponent=3, clamp=TRUE, nsim=9)
})
}
#'
#'    tests/fasteval.R
#'
#'    Checks validity of fast C implementations of Gibbs interaction terms
#'
#'    $Revision: 1.8 $  $Date: 2026/01/18 10:27:23 $
#'
if(FULLTEST) {  # depends on hardware
local({
  X <- redwood
  Q <- quadscheme(X)
  U <- union.quad(Q)
  E <- equalpairs.quad(Q)
  checkit <- function(A, toler=sqrt(.Machine$double.eps)) {
    ## 'A' is an interaction object
    Aname <- deparse(substitute(A))
    a <- A$family$eval(X, U, E, A$pot, A$par, "border")
    b <-    A$fasteval(X, U, E, A$pot, A$par, "border")
    attr(a, "POT") <- NULL
    fa <- is.finite(a)
    fb <- is.finite(b)
    if(!all(fa & fb)) {
      if(any(fa != fb))
        stop(paste0(Aname, "$family$eval and ",
                    Aname, "$fasteval do not agree on Infinite values"))
      a <- a[fa]
      b <- b[fb]
    }
    if(max(abs(a-b)) > toler) {
      stop(paste0(Aname, "$family$eval and ",
                  Aname, "$fasteval do not agree"))
    }
    range(a-b)
  }
  R <- 0.11
  #' R=0.11 is chosen to avoid hardware numerical effects (gcc bug 323).
  #' It avoids being close any value of pairdist(redwood).
  #' The nearest such values are 0.1077.. and 0.1131..
  #' By contrast if r = 0.1 there are values differing from 0.1 by 3e-17
  #' ................. Strauss interaction ...........................
  checkit(Strauss(R))
  #' ................. Strauss-Hard core interaction ...........................
  checkit(StraussHard(R, 0.015))
  #' ................. Geyer interaction ...........................
  checkit(Geyer(R, 2))
  #' and again for a non-integer value of 'sat'
  #' (spotted by Thordis Linda Thorarinsdottir)  
  checkit(Geyer(R, 2.5))
  #' and again for sat < 1
  #' (spotted by Rolf)  
  checkit(Geyer(R, 0.5))
  #' ................. Penttinen interaction ...........................
  checkit(Penttinen(R/2))
})
}

#' tests/formuli.R
#'
#'  Test machinery for manipulating formulae
#' 
#' $Revision: 1.7 $  $Date: 2020/04/28 12:58:26 $

local({

  ff <- function(A, deletevar, B) {
    D <- reduceformula(A, deletevar)
    if(!spatstat.utils::identical.formulae(D, B)) {
      AD <- as.expression(substitute(reduceformula(A,d),
                                     list(A=A, d=deletevar)))
      stop(paste(AD, "\n\tyields ", spatstat.utils::pasteFormula(D),
                 " instead of ", spatstat.utils::pasteFormula(B)),
           call.=FALSE)
    }
    invisible(NULL)
  }

  ff(~ x + z, "x", ~z)

  ff(y ~ x + z, "x", y~z)

  ff(~ I(x^2) + z, "x",  ~z)

  ff(y ~ poly(x,2) + poly(z,3), "x", y ~poly(z,3))

  ff(y ~ x + z, "g", y ~ x + z)

  reduceformula(y ~ x+z, "g", verbose=TRUE)
  reduceformula(y ~ sin(x-z), "z", verbose=TRUE)
  
  illegal.iformula(~str*g, itags="str", dfvarnames=c("marks", "g", "x", "y"))
})



##  
##     tests/funnymarks.R
##
## tests involving strange mark values
## $Revision: 1.7 $ $Date: 2020/04/28 12:58:26 $

if(ALWAYS) { # depends on locale
local({
  ## ppm() where mark levels contain illegal characters
  hyphenated <- c("a", "not-a")
  spaced <- c("U", "non U")
  suffixed <- c("a+", "a*")
  charred <- c("+", "*")

  irad <- matrix(0.1, 2,2)
  hrad <- matrix(0.005, 2, 2)

  tryit <- function(types, X, irad, hrad) { 
    levels(marks(X)) <- types
    fit <- ppm(X ~marks + polynom(x,y,2),
               MultiStraussHard(types=types,iradii=irad,hradii=hrad))
    print(fit)
    print(coef(fit))
    val <- fitted(fit)
    pred <- predict(fit)
    return(invisible(NULL))
  }

  tryit(hyphenated, amacrine, irad, hrad)
  tryit(spaced, amacrine, irad, hrad)
  tryit(suffixed, amacrine, irad, hrad)
  tryit(charred, amacrine, irad, hrad)

  ## marks which are dates
  X <- cells
  n <- npoints(X)
  endoftime <- rep(ISOdate(2001,1,1), n)
  eotDate   <- rep(as.Date("2001-01-01"), n)
  markformat(endoftime)
  markformat(eotDate)
  marks(X) <- endoftime
  print(X)
  Y <- X %mark% data.frame(id=1:42, date=endoftime, dd=eotDate)
  print(Y)
  md <- markformat(endoftime)
  
  ## mark formats
  Z <- Y
  marks(Z) <- marks(Z)[1,,drop=FALSE]
  ms <- markformat(solist(cells, redwood))
  marks(Z) <- factor(1:npoints(Z))
  marks(Z)[12] <- NA
  mz <- is.multitype(Z)
  cZ <- coerce.marks.numeric(Z)
  marks(Z) <- data.frame(n=1:npoints(Z),
                         a=factor(sample(letters, npoints(Z), replace=TRUE)))
  cZ <- coerce.marks.numeric(Z)
  stopifnot(is.multitype(cells %mark% data.frame(a=factor(1:npoints(cells)))))

  a <- numeric.columns(finpines)
  b1 <- numeric.columns(amacrine)
  b2 <- coerce.marks.numeric(amacrine)
  d <- numeric.columns(cells)
  f <- numeric.columns(longleaf)
  ff <- data.frame(a=factor(letters[1:10]), y=factor(sample(letters, 10)))
  numeric.columns(ff)

  ## mark operations
  df <- data.frame(x=1:2, y=sample(letters, 2))
  h <- hyperframe(z=1:2, p=solist(cells, cells))
  a <- NULL %mrep% 3
  a <- 1:4 %mrep% 3
  a <- df %mrep% 3
  a <- h %mrep% 3
  b <- markcbind(df, h)
  b <- markcbind(h, df)
})
}
