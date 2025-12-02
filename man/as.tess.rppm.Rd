\name{as.tess.rppm}
\alias{as.tess.rppm}
\title{
  Convert Recursively Partitioned Point Process Model to a Tessellation
}
\description{
  Given a recursively partitioned point process model,
  create a spatial tessellation representing the partition.
}
\usage{
\method{as.tess}{rppm}(X)
}
\arguments{
  \item{X}{
    A recursively partitioned point process model
    (object of class \code{"rppm"} created by \code{\link{rppm}}).
  }
}
\details{
  The function \code{\link[spatstat.geom]{as.tess}} is generic. The function
  \code{as.tess.rppm} is the method for class \code{"rppm"}.
  
  In a recursively partitioned point process model,
  the intensity function 
  of the point process can be evaluated at any spatial location
  by following a decision tree. Each fork in the tree is defined by a split
  on the value of one of the spatial covariates. The tree is followed
  until reaching a leaf (terminal node) and this determines the value
  of the intensity.
  The tree itself can be visualised by \code{\link{plot.rppm}}.

  This function \code{as.tess.rppm} produces a spatial representation
  of the model. It interprets each fork in the tree as a division
  of space into two complementary subsets. Each leaf (terminal node)
  is an intersection of several such sets. The model corresponds to a
  division of space into disjoint subsets, that is, a tessellation.

  The result is a tessellation in which each tile
  corresponds to a leaf of the decision tree.
  The \code{tilenames} of the tessellation are expressions
  representing the logical conditions which define each tile
  (redundant conditions are omitted).
  The \code{marks} of the tessellation are the intensity values.
}
\value{
  A tessellation (object of class \code{"tess"})
  with numeric marks.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rppm}}
}
\examples{
  fit <- rppm(nztrees ~ x + y)
  V <- as.tess(fit)
  plot(V, do.col=TRUE)
  tilenames(V)
  marks(V)
}
\keyword{spatial}
\keyword{models}
