\name{dropply}
\alias{dropply}
\title{
  Significance Tests for all Single Term Deletions from a Model
}
\description{
  Given a fitted point process model,
  consider dropping each possible term in the model,
  and apply a significance test for the dropped covariate.
}
\usage{
dropply(object,
        action = c("berman.test", "cdf.test", "rhohat", "roc"),
        scope = NULL, \dots, high = TRUE)
}
\arguments{
  \item{object}{
    A fitted point process model (object of class \code{"ppm"},
    \code{"kppm"}, \code{"dppm"}, \code{"slrm"} or \code{"lppm"}).
  }
  \item{action}{
    Character string (partially matched) specifying
    the hypothesis test to be performed, or other calculation.
  }
  \item{scope}{
    A formula or a character vector specifying the model term or terms that
    are to be considered for deletion, or a fitted model obtained by
    deleting some terms from \code{object}.
    The default is to consider all possible terms.
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
  This function is like \code{\link[stats]{drop1}}
  in that it considers each possible term in the model \code{object}
  (or only the terms listed in the \code{scope} argument),
  deletes each such term from the model, and measures the change in the model.
  In this case the change is measured by performing the \code{action}.
  Options are:
  \describe{
    \item{\code{action="roc"}:}{
      the ROC curve for the deleted covariate
      is computed using the updated model as a baseline.
    }
    \item{\code{action="berman.test"}:}{
      One of Berman's tests is applied
      (see \code{\link[spatstat.explore]{berman.test}}), using the updated model
      as the null hypothesis, and the original model 
      as the alternative. 
    }
    \item{\code{action="cdf.test"}:}{
      One of the CDF tests is applied
      (see \code{\link[spatstat.explore]{cdf.test}}), using the updated model
      as the null hypothesis, and the original model 
      as the alternative.
    }
    \item{\code{action="rhohat"}:}{
      taking the updated model as a baseline, the true intensity
      (ratio of true intensity to baseline intensity)
      is estimated as a function of the deleted covariate,
      using the function \code{\link[spatstat.explore]{rhohat}}.
    }
  }
  Note that \code{dropply(model, "roc", \dots)} is equivalent to
  \code{\link[spatstat.model]{dropROC}(model, \dots)}.
}
\value{
  A named list containing the results for each possible deletion.

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
  \code{\link[spatstat.model]{dropROC}},
  \code{\link[spatstat.model]{roc.ppm}},
  \code{\link[spatstat.explore]{rhohat}},
  \code{\link[spatstat.explore]{berman.test}},
  \code{\link[spatstat.explore]{cdf.test}}
}
\examples{
  dimyx <- if(interactive()) NULL else 32
  fut <- ppm(bei ~ grad + elev, data=bei.extra)
  z <- dropply(fut, "b", dimyx=dimyx)
  z
  plot(z, mar.panel=5)

  ## how to extract p-values from each test
  sapply(z, getElement, name="p.value")
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{htest}
