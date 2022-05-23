##  spatstat.model/R/First.R

.onLoad <- function(...) reset.spatstat.options()

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.model"),
                 fields="Version")
  vs <- as.character(vs)
  putSpatstatVariable("SpatstatModelVersion", vs)
  packageStartupMessage(paste("spatstat.model", vs))
  return(invisible(NULL))
}

  
