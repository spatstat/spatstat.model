#
#
#    poisson.S
#
#    $Revision: 1.9 $	$Date: 2022/03/07 03:58:22 $
#
#    The Poisson process
#
#    Poisson()    create an object of class 'interact' describing
#                 the (null) interpoint interaction structure
#                 of the Poisson process.
#	
#
# -------------------------------------------------------------------
#	

Poisson <- local({

  BlankPoisson <- list(
    name     = "Poisson process",
    creator  = "Poisson",
    family   = NULL,
    order    = 1,
    pot      = NULL,
    par      = NULL,
    parnames = NULL,
    init     = function(...) { },
    update   = function(...) { },
    print    = function(self) {
      cat("Poisson process\n")
      invisible()
    },
    valid = function(...) { TRUE },
    project = function(...) NULL, 
    irange = function(...) { 0 },
    version=NULL
    )
  
  class(BlankPoisson) <- "interact"

  Poisson <- function() { BlankPoisson }

  Poisson <- intermaker(Poisson, BlankPoisson)

  Poisson
})
                 
