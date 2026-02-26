#
#  Kmodel.R
#
# Kmodel and pcfmodel
#
#  $Revision: 1.1 $  $Date: 2011/05/30 14:02:21 $
#

Kmodel <- function(model, ...) {
  UseMethod("Kmodel")
}

Lmodel <- function(model, ...) {
  UseMethod("Lmodel")
}

pcfmodel <- function(model, ...) {
  UseMethod("pcfmodel")
}
