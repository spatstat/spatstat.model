#
#  cdftest.R
#
#  $Revision: 2.27 $  $Date: 2022/05/22 11:01:37 $
#
#

## generic cdf.test() and cdf.test.ppp() are moved to spatstat.explore


cdf.test.ppm <- 
  function(model, covariate, test=c("ks", "cvm", "ad"), ...,
           interpolate=TRUE, jitter=TRUE, nsim=99, verbose=TRUE) {
  modelname <- short.deparse(substitute(model))
  covname <- singlestring(short.deparse(substitute(covariate)))
  test <- match.arg(test)
  verifyclass(model, "ppm")
  if(is.character(covariate)) covname <- covariate
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call(spatialCDFtest,
          resolve.defaults(list(model=quote(model), 
				covariate=quote(covariate), 
				test=test),
                           list(interpolate=interpolate, jitter=jitter,
                                nsim=nsim, verbose=verbose),
                           list(...),
                           list(modelname=modelname,
                                covname=covname)))
}





cdf.test.slrm <- function(model, covariate,
                          test=c("ks", "cvm", "ad"), ...,
                          modelname=NULL, covname=NULL) {
  # get names
  if(is.null(modelname))
    modelname <- short.deparse(substitute(model))
  if(is.null(covname))
    covname <- short.deparse(substitute(covariate))
  dataname <- model$CallInfo$responsename
  test <- match.arg(test)
  #
  stopifnot(is.slrm(model))
  stopifnot(is.im(covariate))
  # extract data
  prob <- fitted(model)
  covim <- as.im(covariate, W=as.owin(prob))
  probvalu <- as.matrix(prob)
  covvalu  <- as.matrix(covim)
  ok <- !is.na(probvalu) & !is.na(covvalu)
  probvalu <- as.vector(probvalu[ok])
  covvalu <- as.vector(covvalu[ok])
  # compile weighted cdf's
  FZ <- ewcdf(covvalu, probvalu/sum(probvalu))
  X <- model$Data$response
  ZX <- safelookup(covim, X)
  # Ensure support of cdf includes the range of the data
  xxx <- knots(FZ)
  yyy <- FZ(xxx)
  if(min(xxx) > min(ZX)) {
    xxx <- c(min(ZX), xxx)
    yyy <- c(0, yyy)
  }
  if(max(xxx) < max(ZX)) {
    xxx <- c(xxx, max(ZX))
    yyy <- c(yyy, 1)
  }
  if(length(xxx) > 1) {
    #' non-degenerate cdf
    ## replace by piecewise linear approximation
    FZ <- approxfun(xxx, yyy, rule=2)
  }
  # now apply cdf
  U <- FZ(ZX)
  # Test uniformity of transformed values
  result <- switch(test,
                   ks  = ks.test(U, "punif", ...),
                   cvm = cvm.test(U, "punif", ...),
                   ad = ad.test(U, "punif", ...))
  testname <- switch(test,
                     ks="Kolmogorov-Smirnov",
                     cvm="Cramer-Von Mises",
                     ad="Anderson-Darling")

  # modify the 'htest' entries
  result$method <- paste("Spatial", testname, "test of",
                         "inhomogeneous Poisson process",
                         "in two dimensions")
  result$data.name <-
    paste("covariate", sQuote(paste(covname, collapse="")),
          "evaluated at points of", sQuote(dataname),
          "\n     and transformed to uniform distribution under",
          sQuote(modelname))
  # additional class 'cdftest'
  class(result) <- c("cdftest", class(result))
  attr(result, "prep") <-
    list(Zvalues=covvalu, ZX=ZX, FZ=FZ, FZX=ecdf(ZX), U=U)
  attr(result, "info") <- list(modelname=modelname, covname=covname,
                               dataname=dataname, csr=FALSE)
  return(result)        
}

## Internal code is moved to spatstat.explore

