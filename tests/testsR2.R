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
# tests/rmhExpand.R
#
# test decisions about expansion of simulation window
#
#  $Revision: 1.8 $  $Date: 2022/01/05 02:06:32 $
#

local({
  if(FULLTEST) {
  }
})



#
# tests/rmhTrend.R
#
#  Problems with trend images (rmhmodel.ppm or rmhEngine)
#

if(ALWAYS) {
local({
  set.seed(42)

  # Bug folder 37 of 8 feb 2011
  # rmhmodel.ppm -> predict.ppm
  # + rmhResolveTypes -> is.subset.owin

  Z <- rescale(demopat, 7000)
  X <- unmark(Z)
  X1 <- split(Z)[[1]]
  Int  <- density(X,dimyx=200)
  Lint <- eval.im(log(npoints(X1)*Int/npoints(X)))
  M    <- as.owin(Int)
  MR   <- intersect.owin(M,scalardilate(M,0.5,origin="midpoint"))
  X1 <- X1[MR]
  Fut  <- ppm(X1~offset(Lint),covariates=list(Lint=Lint),
              inter=BadGey(r=c(0.03,0.05),sat=3))
  Y   <- rmh(Fut,control=list(expand=M,nrep=1e3), verbose=FALSE)

})
}
#
#      tests/rmhmodel.ppm.R
#
#    $Revision: 1.10 $  $Date: 2020/05/01 05:29:42 $
#
# Case-by-case tests of rmhmodel.ppm
#

if(FULLTEST) {
local({
f <- ppm(cells)
m <- rmhmodel(f)

f <- ppm(cells ~x)
m <- rmhmodel(f)

f <- ppm(cells ~1, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells ~1, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)
print(m)

f <- ppm(cells ~1, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells ~1, DiggleGratton(0.05,0.1))
m <- rmhmodel(f)

f <- ppm(cells ~1, Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

f <- ppm(cells ~1, Geyer(0.07,2))
m <- rmhmodel(f)

f <- ppm(cells ~1, BadGey(c(0.07,0.1,0.13),2))
m <- rmhmodel(f)

f <- ppm(cells ~1, PairPiece(r = c(0.05, 0.1, 0.2)))
m <- rmhmodel(f)

f <- ppm(cells ~1, AreaInter(r=0.06))
m <- rmhmodel(f)
print(m)

# multitype

r <- matrix(0.07, 2, 2)
f <- ppm(amacrine ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)
print(m)

h <- matrix(min(nndist(amacrine))/2, 2, 2)
f <- ppm(amacrine ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

diag(r) <- NA
diag(h) <- NA
f <- ppm(amacrine ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

f <- ppm(amacrine ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

# multitype data, interaction not dependent on type

f <- ppm(amacrine ~marks, Strauss(0.05))
m <- rmhmodel(f)
print(m)

# trends

f <- ppm(cells ~x, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells ~y, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells ~x+y, Hardcore(0.07))
m <- rmhmodel(f)
print(m)

f <- ppm(cells ~polynom(x,y,2), Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

# covariates

Z <- as.im(function(x,y){ x^2+y^2 }, as.owin(cells))
f <- ppm(cells ~z, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))
print(m)

Zim <- as.im(Z, as.owin(cells))
f <- ppm(cells ~z, covariates=list(z=Zim))
m <- rmhmodel(f)

Z <- as.im(function(x,y){ x^2+y }, as.owin(amacrine))
f <- ppm(amacrine ~z + marks, covariates=list(z=Z))
m <- rmhmodel(f)
print(m)
m <- rmhmodel(f, control=list(p=1))
m <- rmhmodel(f, control=list(p=1,fixall=TRUE))
print(m)

Zim <- as.im(Z, as.owin(amacrine))
f <- ppm(amacrine ~z + marks, covariates=list(z=Zim))
m <- rmhmodel(f)
print(m)

})
}
#
#    tests/rmhmodelHybrids.R
#
#  Test that rmhmodel.ppm and rmhmodel.default
#  work on Hybrid interaction models
#
#   $Revision: 1.5 $  $Date: 2020/05/01 05:29:42 $
#

if(ALWAYS) { # involves C code
local({
    


})
}

#
#  tests/rmh.ppm.R
#
#  $Revision: 1.5 $ $Date: 2020/05/01 05:29:42 $
#
#  Examples removed from rmh.ppm.Rd
#  stripped down to minimal tests of validity
#

local({
   op <- spatstat.options()
   spatstat.options(rmh.nrep=10, npixel=10, ndummy.min=10)
   spatstat.options(project.fast=TRUE)
   Nrep <- 10

   X <- swedishpines
   if(FULLTEST) {
     ## Poisson process
     fit <- ppm(X ~1, Poisson())
     Xsim <- rmh(fit)
   }
   if(ALWAYS) { # Gibbs model => C code
     ## Strauss process   
     fit <- ppm(X ~1, Strauss(r=7))
     Xsim <- rmh(fit)

     ## Strauss process simulated on a larger window
     ## then clipped to original window
     Xsim <- rmh(fit, control=list(nrep=Nrep, expand=1.1, periodic=TRUE))

     ## Extension of model to another window (thanks to Tuomas Rajala)
     Xsim <- rmh(fit, w=square(2))
     Xsim <- simulate(fit, w=square(2))
   
     ## Strauss - hard core process
     ##   fit <- ppm(X ~1, StraussHard(r=7,hc=2))
     ##   Xsim <- rmh(fit, start=list(n.start=X$n))

     ## Geyer saturation process
     ##   fit <- ppm(X ~1, Geyer(r=7,sat=2))
     ##   Xsim <- rmh(fit, start=list(n.start=X$n))

     ## Area-interaction process
     fit <- ppm(X ~1, AreaInter(r=7))
     Xsim <- rmh(fit, start=list(n.start=X$n))
  
     ## Penttinen process
     fit <- ppm(X ~1, Penttinen(r=7))
     Xsim <- rmh(fit, start=list(n.start=X$n))
  
     ## soft core interaction process
     ##     X <- quadscheme(X, nd=50)
     ##     fit <- ppm(X ~1, Softcore(kappa=0.1), correction="isotropic")
     ##     Xsim <- rmh(fit, start=list(n.start=X$n))

     ## Diggle-Gratton pairwise interaction model
     ##     fit <- ppm(cells ~1, DiggleGratton(0.05, 0.1))
     ##     Xsim <- rmh(fit, start=list(n.start=cells$n))
     ##     plot(Xsim, main="simulation from fitted Diggle-Gratton model")
   

     ## piecewise-constant pairwise interaction function
     X <- rSSI(0.05, 100)
     fit <- ppm(X ~1, PairPiece(seq(0.02, 0.1, by=0.01)))
     Xsim <- rmh(fit)
   }
   
   ## marked point pattern
   Y <- amacrine

   if(FULLTEST) {
     #' marked Poisson models
     fit <- ppm(Y)
     Ysim <- rmh(fit)

     fit <- ppm(Y~marks)
     Ysim <- rmh(fit)

     fit <- ppm(Y~x)
     Ysim <- rmh(fit)
     
     fit <- ppm(Y~marks+x)
     Ysim <- rmh(fit)
   }

   if(ALWAYS) {
     #' multitype Strauss
     typ <- levels(Y$marks)
     MS <- MultiStrauss(types = typ,
                        radii=matrix(0.07, ncol=2, nrow=2))

     fit <- ppm(Y~marks*x, MS)
     Ysim <- rmh(fit)

     #' multitype Hardcore
     h0 <- minnndist(unmark(Y)) * 0.95
     MH <- MultiHard(types = typ,
                     hradii=matrix(h0, ncol=2, nrow=2))
     fit <- ppm(Y ~ marks+x, MH)
     Ysim <- rmh(fit)
     #' other code blocks
     Ysim <- rmh(fit, control=list(periodic=TRUE, expand=1))
     Ysim <- rmh(fit, control=list(periodic=FALSE, expand=1))
     #' multihard core with invalid initial state
     Ydouble <- superimpose(Y, rjitter(Y, h0/10))
     Ysim <- rmh(fit, start=list(x.start=Ydouble))

     #' Lennard-Jones
     fut <- ppm(unmark(longleaf) ~ 1, LennardJones(), rbord=1)
     Ysim <- rmh(fut)
     Ysim <- rmh(fut, control=list(periodic=TRUE, expand=1))
   }
   
   spatstat.options(op)
 })


reset.spatstat.options()
#'
#'   tests/rmhsnoopy.R
#'
#'   Test the rmh interactive debugger
#' 
#'   $Revision: 1.10 $  $Date: 2020/05/01 05:29:42 $

if(ALWAYS) { # may depend on platform
local({
  R <- 0.1
  ## define a model and prepare to simulate
  W <- Window(amacrine)
  t1 <- as.im(function(x,y){exp(8.2+0.22*x)}, W)
  t2 <- as.im(function(x,y){exp(8.3+0.22*x)}, W)
  model <- rmhmodel(cif="strauss",
                    trend=solist(off=t1, on=t2),
                    par=list(gamma=0.47, r=R, beta=c(off=1, on=1)))
  siminfo <- rmh(model, preponly=TRUE)
  Wsim <- siminfo$control$internal$w.sim
  Wclip <- siminfo$control$internal$w.clip
  if(is.null(Wclip)) Wclip <- Window(cells)

  ## determine debugger interface panel geometry
  Xinit <- runifpoint(ex=amacrine)[1:40]
  P <- rmhsnoop(Wsim=Wsim, Wclip=Wclip, R=R,
                xcoords=Xinit$x,
                ycoords=Xinit$y,
                mlevels=levels(marks(Xinit)),
                mcodes=as.integer(marks(Xinit)) - 1L,
                irep=3L, itype=1L,
                proptype=1, proplocn=c(0.5, 0.5), propmark=0, propindx=0,
                numerator=42, denominator=24,
                panel.only=TRUE)
  boxes <- P$boxes
  clicknames <- names(P$clicks)
  boxcentres <- do.call(concatxy, lapply(boxes, centroid.owin))

  ## design a sequence of clicks
  actionsequence <- c("Up", "Down", "Left", "Right",
                      "At Proposal", "Zoom Out", "Zoom In", "Reset",
                      "Accept", "Reject", "Print Info",
                      "Next Iteration", "Next Shift", "Next Death",
                      "Skip 10", "Skip 100", "Skip 1000", "Skip 10,000",
                      "Skip 100,000", "Exit Debugger")
  actionsequence <- match(actionsequence, clicknames)
  actionsequence <- actionsequence[!is.na(actionsequence)]
  xy <- lapply(boxcentres, "[", actionsequence)

  ## queue the click sequence
  spatstat.utils::queueSpatstatLocator(xy$x,xy$y)

  ## go
  rmh(model, snoop=TRUE)
})
}
