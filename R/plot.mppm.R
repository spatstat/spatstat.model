#
# plot.mppm.R
#
#   $Revision: 1.8 $  $Date: 2023/02/28 04:04:51 $
#
#

plot.mppm <- function(x, ..., trend=TRUE, cif=FALSE, se=FALSE,
                      how=c("image", "contour", "persp"), main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  how <- match.arg(how)
  subs <- subfits(x)
  if(trend) 
    plot.anylist(x=subs, how=how, main=main, ...,
                 trend=TRUE, cif=FALSE, se=FALSE)
  if(cif) 
    plot.anylist(x=subs, how=how, main=main, ...,
                 trend=FALSE, cif=TRUE, se=FALSE)
  if(se) 
    plot.anylist(x=subs, how=how, main=main, ...,
                 trend=FALSE, cif=FALSE, se=TRUE)
  invisible(NULL)
}

