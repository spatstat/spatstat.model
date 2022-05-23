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
#   tests/testaddvar.R
#
# test addvar options
#
#   $Revision: 1.3 $  $Date: 2020/05/02 01:32:58 $

if(FULLTEST) {
local({
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  model <- ppm(X ~y)
  addvar(model, "x", crosscheck=TRUE)
  addvar(model, "x", bw.input="quad")
  w <- square(0.5)
  addvar(model, "x", subregion=w)
  addvar(model, "x", subregion=w, bw.input="points")
  Z <- as.im(function(x,y) { x }, Window(X))
  addvar(model, Z)
})
}
#
#   tests/testparres.R
#
# additional test of parres
#
#  $Revision: 1.8 $  $Date: 2022/05/20 06:59:59 $
#

if(FULLTEST) {
local({
X <-  rpoispp(function(x,y){exp(3+x+2*x^2)})
model <- ppm(X ~x+y)

# options in parres (and code blocks in print.parres)
parres(model, "x")
parres(model, "x", smooth.effect=TRUE)
parres(model, "x", bw.input="quad")
w <- square(0.5)
parres(model, "x", subregion=w)
parres(model, "x", subregion=w, bw.input="quad")
f <- function(x,y) { x + y }
parres(model, f)

# check whether 'update.ppm' has messed up internals
mod2 <- update(model, ~x)
parres(mod2, "x")

#' other kinds of covariates
mod3 <- ppm(X ~ x + offset(y))
parres(mod3, "offset(y)")
Z <- distmap(runifpoint(3))
parres(mod3, Z)
mod4 <- ppm(X ~ sin(x), data=solist(B=Z))
parres(mod4, "sin(x)")
parres(mod4, "B")

#' models with interaction
mod5 <- ppm(cells ~ x, AreaInter(0.06))
parres(mod5, "x")
dlin <- distfun(copper$SouthLines)
copfit <- ppm(copper$SouthPoints ~ dlin, Geyer(1,1))
parres(copfit, "dlin")

#' covariate need not be specified if there is only one.
parres(mod5)
parres(copfit)

#' infrastructure
ltuae <- evaluateCovariate(42, cells)
LTUAE <- evaluateCovariate(ltuae, cells)

fit <- ppm(amacrine ~ x * marks, nd=16)
dmat <- model.depends(fit)
check.separable(dmat, "x", c(x=FALSE, marks=FALSE), FALSE)
check.separable(dmat, "x", c(FALSE, FALSE), FALSE)
check.separable(dmat, "x", c(x=FALSE, marks=TRUE), FALSE)
})
}
#'
#'     tests/threedee.R
#'
#'     Tests of 3D code 
#'
#'      $Revision: 1.8 $ $Date: 2020/05/02 01:32:58 $
#'

local({
  X <- runifpoint3(30)
  Y <- runifpoint3(20)
  if(FULLTEST) {
    A <- runifpoint3(10, nsim=2)
    Z <- ppsubset(X, 2:4)
  }
  ##
  if(ALWAYS) { # includes C code
    d <- pairdist(X, periodic=TRUE, squared=TRUE)
    d <- crossdist(X, Y, squared=TRUE)
    d <- crossdist(X, Y, squared=TRUE, periodic=TRUE)
    #' 
    h <- has.close(X, 0.2)
    h <- has.close(X, 0.2, periodic=TRUE)
    h <- has.close(X, 0.2, Y=Y)
    h <- has.close(X, 0.2, Y=Y, periodic=TRUE)
    #' code blocks not otherwise reached
    rmax <- 0.6 * max(nndist(X))
    g <- G3est(X, rmax=rmax, correction="rs")
    g <- G3est(X, rmax=rmax, correction="km")
    g <- G3est(X, rmax=rmax, correction="Hanisch")
    g <- G3est(X, rmax=rmax, sphere="ideal")
    g <- G3est(X, rmax=rmax, sphere="digital")
    v <- sphere.volume()
    v <- digital.volume()
    #' older code
    co <- coords(X)
    xx <- co$x
    yy <- co$y
    zz <- co$z
    gg1 <- g3engine(xx, yy, zz, correction="Hanisch G3")
    gg2 <- g3engine(xx, yy, zz, correction="minus sampling")
    ff1 <- f3engine(xx, yy, zz, correction="no")
    ff2 <- f3engine(xx, yy, zz, correction="minus sampling")
  }
  ##
  if(ALWAYS) {
    #'class support
    X <- runifpoint3(10)
    print(X)
    print(X %mark% runif(10))
    print(X %mark% factor(letters[c(1:5,5:1)]))
    print(X %mark% data.frame(a=1:10, b=runif(10)))
    da <- as.Date(paste0("2020-01-0", c(1:5,5:1)))
    print(X %mark% da)
    print(X %mark% data.frame(a=1:10, b=da))
  }
})
#
# tests/triplets.R
#
# test code for triplet interaction 
#
# $Revision: 1.9 $ $Date: 2022/05/22 08:45:38 $
#
if(ALWAYS) { # C code, platform dependence
local({
  #' valid model
  fit <- ppm(cells ~1, Triplets(0.1))
  fit
  suffstat(fit)
  #' invalid model 
  fitR <- ppm(redwood ~1, Triplets(0.1))
  fitR
  suffstat(fitR)
  #' hard core (zero triangles, coefficient is NA)
  fit0 <- ppm(cells ~1, Triplets(0.05))
  fit0
  suffstat(fit0)
  #' bug case (1 triangle in data)
  fit1 <- ppm(cells ~1, Triplets(0.15))
  fit1
  suffstat(fit1)
  #' simulation
  X <- simulate(fit)
  mod <- list(cif="triplets",par=list(beta=50,gamma=0.2,r=0.07), w=square(1))
  Xm <- rmh(model=mod,start=list(n.start=5), control=list(nrep=1e5))
  #' hard core
  mod$par$gamma <- 0
  XmHard <- rmh(model=mod,start=list(n.start=5), control=list(nrep=1e5))
})
}
