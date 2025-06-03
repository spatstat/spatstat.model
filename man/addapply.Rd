\name{addapply}
\alias{addapply}
\title{
  Significance Tests or Effect Size for Single Term Additions to a Model
}
\description{
  Given a fitted point process model,
  consider adding new explanatory variables,
  and apply a significance test (or effect size calculation)
  for the effect of each new variable.
}
\usage{
addapply(object, 
         action = c("berman.test", "cdf.test", "rhohat", "roc"),
         scope, 
         \dots, high = TRUE)
}
\arguments{
  \item{object}{
    A fitted point process model (object of class \code{"ppm"},
    \code{"kppm"}, \code{"dppm"}, \code{"slrm"} or \code{"lppm"})
    specifying the model to be extended. 
  }
  \item{action}{
    Character string (partially matched) specifying
    the hypothesis test to be performed, or other calculation.
  }
  \item{scope}{
    A formula or a character vector specifying the variable or variables that
    are to be considered for addition, or a fitted point process model
    containing all these variables.
  }
  \item{high}{
    Argument passed to \code{\link[spatstat.explore]{roc}} to specify whether the
    ROC curves should be based on high or low values of the covariates,
    when \code{action="roc"}.
    Either a logical value, or a logical vector of the same length
    as \code{scope}.
  }
  \item{\dots}{
    Other arguments passed to the relevant function
    \code{\link[spatstat.explore]{roc}}, 
    \code{\link[spatstat.explore]{rhohat}},
    \code{\link[spatstat.explore]{berman.test}} or
    \code{\link[spatstat.explore]{cdf.test}}.
    This includes arguments passed to
    \code{\link[spatstat.geom]{as.mask}} to control the
    pixel resolution used for calculation.
  }
}
\details{
  This function is like \code{\link[stats]{add1}}
  in that it considers adding new terms to the model \code{object}
  and measures the change in the model.
  In this case the change is measured by performing the \code{action}.
  Options are:
  \describe{
    \item{\code{action="roc"}:}{
      the ROC curve for the added covariate
      is computed using the original model \code{object} as a baseline.
    }
    \item{\code{action="berman.test"}:}{
      One of Berman's tests is applied
      (see \code{\link[spatstat.explore]{berman.test}}), using the
      original model \code{object}
      as the null hypothesis, and the extended model
      as the alternative. 
    }
    \item{\code{action="cdf.test"}:}{
      One of the CDF tests is applied
      (see \code{\link[spatstat.explore]{cdf.test}}), using the original
      model \code{object} as the null hypothesis, and the extended model 
      as the alternative.
    }
    \item{\code{action="rhohat"}:}{
      taking the original model \code{object} as a baseline, the true intensity
      (ratio of true intensity to baseline intensity)
      is estimated as a function of the added explanatory variable,
      using the function \code{\link[spatstat.explore]{rhohat}}.
    }
  }
  Note that \code{addapply(object, "roc", scope)}
  is equivalent to \code{\link[spatstat.model]{addROC}(object, scope)}.

  Either \code{object} or \code{scope} should be a fitted point process
  model, and the other argument may be a fitted point process model or a
  formula. If \code{object} is a fitted model then \code{scope} may be a
  character vector of the names of variables to be added.
}
\value{
  A named list containing the results for each added variable.

  The list belongs to the class \code{"anylist"}
  so it can be printed and plotted in its entirety.

  If \code{action="roc"} the individual entries are
  ROC curves belonging to class \code{"fv"}.
  If \code{action="rhohat"} the individual entries are
  curves belonging to class \code{"fv"} and class \code{"rhohat"}.
  If \code{action="berman.test"} the individual entries
  are hypothesis tests of class \code{"htest"} and \code{"bermantest"}.
  If \code{action="cdf.test"} the individual entries
  are hypothesis tests of class \code{"htest"} and \code{"cdftest"}.
}
\author{
  \rocketAuthors.
}
\references{
  \rocpaper
}
\seealso{
  \code{\link[spatstat.model]{addROC}},
  \code{\link[spatstat.model]{dropROC}},
  \code{\link[spatstat.model]{roc.ppm}},
  \code{\link[spatstat.explore]{rhohat}},
  \code{\link[spatstat.explore]{berman.test}},
  \code{\link[spatstat.explore]{cdf.test}}
}
\examples{
  dimyx <- if(interactive()) NULL else 32
  fit0 <- ppm(bei ~ 1, data=bei.extra)
  z <- addapply(fit0, "ber", . ~ grad+elev, dimyx=dimyx)
  z
  plot(z, mar.panel=5)

  ## how to extract p-values from each test
  sapply(z, getElement, name="p.value")
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{htest}
