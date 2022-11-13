#'
#'    traj.R
#'
#'  Trajectories of function evaluations
#'
#'  An object of class 'traj' is a data frame containing the sequence of
#'  function evaluations (input and output) performed by 'optim'.
#'
#'  If kppm is called with pspace$save=TRUE,
#'  the resulting kppm object has an attribute 'h' containing the trajectory.
#'  This is extracted by attr(, "h") or traj()
#'
#'  There are methods for print, plot and lines
#' 
#'  $Revision: 1.3 $  $Date: 2022/11/10 07:03:33 $
#' 
#'  Copyright (c) Adrian Baddeley 2022
#'  GNU Public Licence >= 2.0
#'

traj <- function(object) {
  h <- attr(object, "h")
  if(inherits(h, "traj")) return(h)
  return(NULL)
}

print.traj <- function(x, ...) {
  cat("Trajectory of function evaluations\n")
  cat(paste0("Variables: ", commasep(sQuote(colnames(x))), "\n"))
  invisible(NULL)
}

plot.traj <- function(x, ..., show.ends=TRUE, add=FALSE, xlab=NULL, ylab=NULL) {
  if(add) {
    ## add to existing plot (handles 'type')
    points(x[,1], x[,2], ...)
  } else {
    ## new plot
    nama <- colnames(x)
    plot(x[,1], x[,2],
         xlab=xlab %orifnull% nama[1],
         ylab=ylab %orifnull% nama[2], ...)
  }
  ## indicate first and last states
  if(show.ends) {
    n <- nrow(x)
    points(x[1,1], x[1,2], pch=1, col="blue", cex=3)
    points(x[n,1], x[n,2], pch=3, col="red", cex=2)
  }
  return(invisible(NULL))
}

lines.traj <- function(x, ..., directed=FALSE) {
  xx <- x[,1]
  yy <- x[,2]
  nn <- length(xx)
  if(directed) {
    arrows(xx[-nn], yy[-nn], xx[-1], yy[-1], ...)
  } else {
    lines(x=xx, y=yy, ...)
  }
  return(invisible(NULL))
}

