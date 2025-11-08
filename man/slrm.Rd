\name{slrm}
\alias{slrm}
\title{Spatial Logistic Regression}
\description{
  Fits a spatial logistic regression model
  to a spatial point pattern. 
}
\usage{
slrm(formula, ..., data = NULL, offset = TRUE, link = c("logit", "cloglog"),
                   ruleAtPoints=c("s", "r"), dataAtPoints=NULL, splitby=NULL)
}
\arguments{
  \item{formula}{The model formula. See Details.}
  \item{\dots}{
    Optional arguments passed to \code{\link[spatstat.geom]{as.mask}}
    determining the pixel resolution for the discretisation
    of the point pattern.
  }
  \item{data}{
    Optional. A list containing data required in the formula.
    The names of entries in the list should correspond to variable
    names in the formula. The entries should be point patterns,
    pixel images or windows.
  }
  \item{offset}{
    Logical flag indicating whether the model formula
    should be augmented by an offset equal to the logarithm of the
    pixel area.
  }
  \item{link}{The link function for the regression model.
    A character string (partially matched), specifying the link function
    for binary regression.
  }
  \item{ruleAtPoints}{
    Rule for obtaining covariate data at the observed point locations.
    See Details.
    Ignored if \code{dataAtPoints} is given.
  }
  \item{dataAtPoints}{Optional.
    Exact values of the covariates at the data points.
    A data frame, with column names corresponding to
    variables in the \code{formula}, with one row for each
    point in the point pattern dataset.
  }
  \item{splitby}{
    Optional. Character string identifying a window. The window will be used
    to split pixels into sub-pixels. 
  }
}
\details{
  This function fits a Spatial Logistic Regression model
  (Tukey, 1972; Agterberg, 1974) to a spatial point pattern dataset. 
  The logistic function may be replaced by another link function.

  The \code{formula} specifies the form of the model to be fitted,
  and the data to which it should be fitted. The \code{formula}
  must be an \R formula with a left and right hand
  side.

  The left hand side of the \code{formula} is the name of the
  point pattern dataset, an object of class \code{"ppp"}. 

  The right hand side of the \code{formula} is an expression,
  in the usual \R formula syntax, representing the functional form of
  the linear predictor for the model.

  Each variable name that appears in the formula may be 
  \itemize{
    \item
    one of the reserved names \code{x} and \code{y},
    referring to the Cartesian coordinates;
    \item
    the name of an entry in the list \code{data}, if this argument is given;
    \item
    the name of an object in the
    parent environment, that is, in the environment where the call
    to \code{slrm} was issued.
  }
  Each object appearing on the right hand side of the formula may be
  \itemize{
    \item a pixel image (object of class \code{"im"})
    containing the values of a covariate;
    \item a window (object of class \code{"owin"}), which will be
    interpreted as a logical covariate which is \code{TRUE} inside the
    window and \code{FALSE} outside it;
    \item a \code{function} in the \R language, with arguments
    \code{x,y}, which can be evaluated at any location to
    obtain the values of a covariate;
    \item a tessellation (object of class \code{"tess"}), which will be
    interpreted as a factor-valued covariate, with each level of the
    factor associated with a tile of the tessellation.
  }
  See the Examples below.

  The fitting algorithm discretises the point pattern onto a pixel grid. The
  value in each pixel is 1 if there are any points of the point pattern
  in the pixel, and 0 if there are no points in the pixel.
  The dimensions of the pixel grid will be determined as follows:
  \itemize{
    \item
    The pixel grid will be determined by the extra
    arguments \code{\dots} if they are specified (for example the argument
    \code{dimyx} can be used to specify the number of pixels).
    \item
    Otherwise, if the right hand side of the \code{formula} includes
    the names of any pixel images containing covariate values,
    these images will determine the pixel grid for the discretisation.
    The covariate image with the finest grid (the smallest pixels) will
    be used.
    \item
    Otherwise, the default pixel grid size is given by
    \code{spatstat.options("npixel")}.
  }

  The covariates mentioned in the right hand side of the model formula
  are evaluated at each pixel, as described in the section on
  \emph{Discretisation rule}.
  
  If \code{link="logit"} (the default), the algorithm fits a Spatial Logistic
  Regression model. This model states that the probability
  \eqn{p} that a given pixel contains a data point, is related to the
  covariates through
  \deqn{\log\frac{p}{1-p} = \eta}{log(p/(1-p)) = eta}
  where \eqn{\eta}{eta} is the linear predictor of the model
  (a linear combination of the covariates,
  whose form is specified by the \code{formula}).

  If \code{link="cloglog"} then the algorithm fits a model stating that
  \deqn{\log(-\log(1-p)) = \eta}{log(-log(1-p)) = eta}.

  If \code{offset=TRUE} (the default), the model formula will be
  augmented by adding an offset term equal to the logarithm of the pixel
  area. This ensures that the fitted parameters are
  approximately independent of pixel size.
  If \code{offset=FALSE}, the offset is not included, and the
  traditional form of Spatial Logistic Regression is fitted.
}
\section{Discretisation rule}{
  Results depend on the rule used to convert continuous vector data
  to pixel data.
  
  The model-fitting procedure requires the values of the covariates
  mentioned on the right hand side of the model formula,
  evaluated for each pixel.

  Covariates are initially evaluated at the centre of each pixel.

  For pixels which contain a point of the spatial point pattern,
  the covariate values will be updated if 
  one of the arguments \code{dataAtPoints} or \code{ruleAtPoints}
  is given:
  \itemize{
    \item 
    If \code{dataAtPoints} is given, then 
    the covariate values for pixels which contain a point of the
    spatial point pattern will be determined by the entries of
    \code{dataAtPoints}
    (and the spatial coordinates \code{x} and \code{y}
    are overwritten by the exact spatial
    coordinates of the data points).
    \item
    If \code{dataAtPoints} is not given, but \code{ruleAtPoints="r"} is
    specified, then the covariate values for pixels which contain a
    point of the spatial point pattern are determined by evaluating the
    covariate at the exact spatial coordinates of the data points.
    This is called rule R in Baddeley et al (2025).
    \item
    If \code{dataAtPoints} is not given and \code{ruleAtPoints="s"}
    (the defaults), covariate values for all pixels are evaluated
    at the centre of each pixel, and there is no special treatment
    of pixels which contain points of the spatial point pattern.
    This is called rule S in Baddeley et al (2025) and is the
    standard practice in many software packages.
  }
  Rule S and rule R are different if one of the covariates is
  either a \code{function(x,y)},
  a polygonal window (object of class \code{"owin"}
  representing a polygon or rectangle), or a polygonal tessellation
  (object of class \code{"tess"} whose tiles are polygons or
  rectangles). In any of these cases, 
  the exact spatial location of a data point
  and the centre of the pixel that contains the data point
  may yield different values of the covariate.
  For a covariate which is a polygonal window or polygonal tessellation,
  the discrepancy is enough to cause instability to changes in
  pixel size when rule S is used, and we recommend rule R in these
  cases (Baddeley et al 2025).

  The best possible stability against
  changes in pixel size is achieved by setting \code{ruleAtPoints="r"},
  \code{link="cloglog"} and \code{offset=TRUE}.
}
\value{
  An object of class \code{"slrm"} representing the fitted model.

  There are many methods for this class, including methods for
  \code{print}, \code{fitted}, \code{predict},
  \code{anova}, \code{coef}, \code{logLik}, \code{terms},
  \code{update}, \code{formula} and \code{vcov}.
  Automated stepwise model selection is possible using
  \code{\link{step}}. Confidence intervals for the parameters can be
  computed using \code{\link[stats]{confint}}. 
}
\seealso{
  \code{\link{anova.slrm}},
  \code{\link{coef.slrm}},
  \code{\link{fitted.slrm}},
  \code{\link{logLik.slrm}},
  \code{\link{plot.slrm}},
  \code{\link{predict.slrm}},
  \code{\link{vcov.slrm}}
}
\references{
  Agterberg, F.P. (1974)
  Automatic contouring of geological maps to detect target areas for
  mineral exploration.
  \emph{Journal of the International Association for Mathematical Geology}
  \bold{6}, 373--395.

  Baddeley, A., Brown, W., Nair, G.M., Milne, R.K., Rakshit, S. and Fu, S.C.
  (2025)
  Mineral prospectivity analysis is unstable to changes in pixel size.
  \emph{Computers and Geosciences}, 105965.
  \code{DOI: 10.1016/j.cageo.2025.105965}.

  Baddeley, A., Berman, M., Fisher, N.I., Hardegen, A., Milne, R.K.,
  Schuhmacher, D., Shah, R. and Turner, R. (2010)
  Spatial logistic regression and change-of-support
  for spatial Poisson point processes.
  \emph{Electronic Journal of Statistics}
  \bold{4}, 1151--1201.
  \code{DOI: 10.1214/10-EJS581}

  Tukey, J.W. (1972)
  Discussion of paper by F.P. Agterberg and S.C. Robinson.
  \emph{Bulletin of the International Statistical Institute}
  \bold{44} (1) p. 596.
  Proceedings, 38th Congress, International Statistical Institute.
}
\examples{
if(offline <- !interactive()) op <- spatstat.options(npixel=32)

     X <- copper$SouthPoints
     slrm(X ~ 1)
     slrm(X ~ x+y)

     slrm(X ~ x+y, link="cloglog")
     # specify a grid of 2-km-square pixels
     slrm(X ~ 1, eps=2)

     Y <- copper$SouthLines
     Z <- distmap(Y)
     slrm(X ~ Z)
     slrm(X ~ Z, dataAtPoints=list(Z=nncross(X,Y,what="dist")))

     mur <- solapply(murchison, rescale, s=1000, unitname="km")
     mur$dfault <- distfun(mur$faults)
     slrm(gold ~ dfault, data=mur)
     slrm(gold ~ dfault, data=mur, ruleAtPoints="r")
     slrm(gold ~ dfault + greenstone, data=mur)
     slrm(gold ~ dfault, data=mur, splitby="greenstone")

if(offline) spatstat.options(op)
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{models}
