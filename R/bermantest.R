#
# bermantest.R
#
# Test statistics from Berman (1986)
#
#  $Revision: 1.24 $  $Date: 2022/05/22 02:23:18 $
#
#

## Code for generic berman.test and berman.test.ppp
## is moved to spatstat.explore


berman.test.ppm <- function(model, covariate,
                           which=c("Z1", "Z2"),
                           alternative=c("two.sided", "less", "greater"),
                           ...) {
  modelname <- short.deparse(substitute(model))
  covname <- short.deparse(substitute(covariate))
  force(model)
  force(covariate)
  if(is.character(covariate)) covname <- covariate
  verifyclass(model, "ppm")
  which <- match.arg(which)
  alternative <- match.arg(alternative)
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call(bermantestEngine,
          resolve.defaults(list(quote(model), 
				quote(covariate), 
				which, alternative),
                           list(...),
                           list(modelname=modelname,
                                covname=covname,
                                dataname=model$Qname)))
}


## Code for generic berman.test and berman.test.ppp
## is moved to spatstat.explore
