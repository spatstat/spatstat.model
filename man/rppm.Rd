\name{rppm}
\alias{rppm}
\title{
  Recursively Partitioned Point Process Model
}
\description{
  Fits a recursive partition model to point pattern data.
}
\usage{
  rppm(\dots, rpargs=list())
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link{ppm}}
    specifying the point pattern data and the explanatory covariates.
  }
  \item{rpargs}{
    Optional list of arguments passed to \code{\link[rpart]{rpart}}
    controlling the recursive partitioning procedure.
  }
}
\details{
  This function attempts to find a simple rule for predicting
  the intensity of a point process using explanatory covariates.
  
  The arguments \code{\dots} specify the point pattern data
  and explanatory covariates in the same way as they would be
  in the function \code{\link{ppm}}. 

  The recursive partitioning algorithm \code{\link[rpart]{rpart}}
  is then used to find a partitioning rule.
  
  In the resulting model, a recursively partitioned point process model,
  the intensity function 
  of the point process can be evaluated at any spatial location
  by following a decision tree. Each fork in the tree is defined by a split
  on the value of one of the spatial covariates at the query location.
  For example, for a numerical covariate \eqn{Z}, the left branch of the fork
  is followed if \eqn{Z < c} and the right fork is followed if
  \eqn{Z \ge c}{Z >= c}, where \eqn{c} is the critical value of the
  split.
  For a categorical covariate, the left branch of the fork is followed
  if the value of \eqn{Z} is one of a specified set of levels of
  \eqn{Z}, and the right branch is followed otherwise.
  The tree is followed until reaching a leaf (terminal node)
  and this determines the value of the intensity.
  
  The decision tree can be visualised as a tree by
  \code{\link{plot.rppm}}.
  The split rules are printed in text form by \code{print.rppm}.
  The predicted intensity at all spatial locations
  is computed by \code{\link{predict.rppm}}. The division of space
  into subsets with different predicted intensities is
  computed by \code{\link{as.tess.rppm}}.
}
\value{
  An object of class \code{"rppm"}. There are methods
  for \code{print}, \code{plot}, \code{fitted}, \code{predict} and
  \code{prune} for this class.
}
\references{
  Breiman, L., Friedman, J. H., Olshen, R. A., and Stone, C. J. (1984)
  \emph{Classification and Regression Trees}. Wadsworth.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{plot.rppm}}, 
  \code{\link{predict.rppm}}, 
  \code{\link{update.rppm}}, 
  \code{\link{prune.rppm}},
  \code{\link{as.tess.rppm}}.
}
\examples{
    # New Zealand trees data: trees planted along border
    # Use covariates 'x', 'y'
    nzfit <- rppm(nztrees ~ x + y)
    nzfit
    nzfit2 <- prune(nzfit, cp=0.035)
    nzfit2
    plot(nzfit2)
    # Murchison gold data: numeric and logical covariates
    mur <- solapply(murchison, rescale, s=1000, unitname="km")
    mur$dfault <- distfun(mur$faults)
    # 
    mfit <- rppm(gold ~ dfault + greenstone, data=mur)
    mfit
    # Gorillas data: factor covariates
    #          (symbol '.' indicates 'all variables')
    gfit <- rppm(unmark(gorillas) ~ . , data=gorillas.extra)
    gfit
}
\keyword{spatial}
\keyword{models}
