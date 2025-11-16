#
#   quadrattest.R
#
#   $Revision: 1.71 $  $Date: 2025/11/16 09:34:27 $
#

## Code for generic quadrat.test() and quadrat.test.ppp()
## is moved to spatstat.explore


quadrat.test.slrm <- quadrat.test.ppm <-
  function(X, nx=5, ny=nx,
           alternative = c("two.sided", "regular", "clustered"),      
           method=c("Chisq", "MonteCarlo"),
           conditional=TRUE, CR=1, df.est=NULL, ...,
           xbreaks=NULL, ybreaks=NULL,
           tess=NULL, nsim=1999)
{
   fitname <- short.deparse(substitute(X))
   dataname <- paste("data from", fitname)
   method <- match.arg(method)
   alternative <- match.arg(alternative)
   if(!is.poisson(X))
    stop("Test is only defined for Poisson point process models")
   if(is.marked(X))
    stop("Sorry, not yet implemented for marked point process models")
   Xdata <- response(X)
   dont.complain.about(Xdata)
   do.call(quadrat.testEngine,
          resolve.defaults(list(quote(Xdata), nx=nx, ny=ny,
                                alternative=alternative,
                                method=method,
                                conditional=conditional, CR=CR,
                                xbreaks=xbreaks, ybreaks=ybreaks,
                                tess=tess,
                                nsim=nsim, 
                                fit=X,
                                df.est=df.est),
                           list(...),
                           list(Xname=dataname, fitname=fitname)))
}



## code for quadrat.test.quadratcount is moved to spatstat.explore
## Infrastructure for quadrat.test is moved to spatstat.explore
