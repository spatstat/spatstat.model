\name{ppm.ppp}
\alias{ppm.ppp}
\alias{ppm.quad}
\concept{point process model}
\concept{Poisson point process}
\concept{Gibbs point process}
\title{
  Fit Point Process Model to Point Pattern Data
}
\description{
  Fits a point process model to an observed point pattern.
}
\usage{
   \method{ppm}{ppp}(Q, trend=~1, interaction=Poisson(),
       \dots,
       covariates=data,
       data=NULL,
       covfunargs = list(),
       subset,
       clipwin,
       correction="border",
       rbord=reach(interaction),
       use.gam=FALSE,
       method=c("mpl", "logi", "VBlogi"),
       forcefit=FALSE,
       improve.type = c("none", "ho", "enet"),
       improve.args=list(),
       emend=project,
       project=FALSE,
       prior.mean = NULL,
       prior.var = NULL,
       nd = NULL,
       eps = NULL,
       quad.args=list(),
       gcontrol=list(),
       nsim=100, nrmh=1e5, start=NULL, control=list(nrep=nrmh),
       verb=TRUE,
       callstring=NULL)

   \method{ppm}{quad}(Q, trend=~1, interaction=Poisson(),
       \dots,
       covariates=data,
       data=NULL,
       covfunargs = list(),
       subset,
       clipwin,
       correction="border",
       rbord=reach(interaction),
       use.gam=FALSE,
       method=c("mpl", "logi", "VBlogi"),
       forcefit=FALSE,
       improve.type = c("none", "ho", "enet"),
       improve.args=list(),
       emend=project,
       project=FALSE,
       prior.mean = NULL,
       prior.var = NULL,
       nd = NULL,
       eps = NULL,
       quad.args=list(),
       gcontrol=list(),
       nsim=100, nrmh=1e5, start=NULL, control=list(nrep=nrmh),
       verb=TRUE,
       callstring=NULL)
}
\arguments{
  \item{Q}{
    A data point pattern (of class \code{"ppp"})
    to which the model will be fitted,
    or a quadrature scheme (of class \code{"quad"})
    containing this pattern.
  }
  \item{trend}{
  An \R formula object specifying the spatial trend to be fitted. 
  The default formula, \code{~1}, indicates the model is stationary
  and no trend is to be fitted. 
  }
  \item{interaction}{
    An object of class \code{"interact"}
    describing the point process interaction
    structure, or a function that makes such an object,
    or \code{NULL} indicating that a Poisson process (stationary
    or nonstationary) should be fitted.
  }
  \item{\dots}{Ignored.}
  \item{data,covariates}{
    The values of any spatial covariates (other than the Cartesian
    coordinates) required by the model.
    Either a data frame, or a list whose entries are images,
    functions, windows, tessellations or single numbers. See Details.
  }
  \item{subset}{
    Optional.
    An expression (which may involve the names of the
    Cartesian coordinates \code{x} and \code{y}
    and the names of entries in \code{data})
    defining a subset of the spatial domain,
    to which the likelihood or pseudolikelihood should be restricted.
    See Details.
    The result of evaluating the expression should be either a logical
    vector, or a window (object of class \code{"owin"})
    or a logical-valued pixel image (object of class \code{"im"}).
  }
  \item{clipwin}{
    Optional. A spatial window (object of class \code{"owin"})
    to which data will be restricted, before model-fitting
    is performed. See Details.
  }
  \item{covfunargs}{
    A named list containing the values of any additional arguments
    required by covariate functions.
  }
  \item{correction}{
    The name of the edge correction to be used. The default 
    is \code{"border"} indicating the border correction.
    Other possibilities may include \code{"Ripley"}, \code{"isotropic"},
    \code{"periodic"}, \code{"translate"} and \code{"none"}, depending on the 
    \code{interaction}.
  }
  \item{rbord}{
    If \code{correction = "border"}
    this argument specifies the distance by which
    the window should be eroded for the border correction.
  }
  \item{use.gam}{
    Logical flag; if \code{TRUE} then computations are performed
    using \code{gam} instead of \code{\link{glm}}.
  }
  \item{method}{
    String (partially matched) specifying the method used to fit the
    model. Options are 
    \code{"mpl"} for the method of Maximum PseudoLikelihood (the default),
    \code{"logi"} for the Logistic Likelihood method and
    \code{"VBlogi"} for the Variational Bayes Logistic Likelihood method.
  }
  \item{forcefit}{
    Logical flag for internal use.
    If \code{forcefit=FALSE}, some trivial models will be
    fitted by a shortcut. If \code{forcefit=TRUE},
    the generic fitting method will always be used. 
  }
  \item{improve.type}{
    String (partially matched) specifying a method for improving the
    initial fit. 
    If \code{improve.type = "none"} (the default), no improvement is
    performed. If \code{improve.type="ho"}, the Huang-Ogata approximate
    maximum likelihood method is used. If \code{improve.type="enet"},
    the model coefficients are re-estimated using a regularized version
    of the composite likelihood.
  }
  \item{improve.args}{
    Arguments used to control the algorithm
    for improving the initial fit. See Details.
  }
  \item{emend,project}{
    (These are equivalent: \code{project} is an older name for
    \code{emend}.)
    Logical value. Setting \code{emend=TRUE} will ensure that the
    fitted model is always a valid point process by
    applying \code{\link{emend.ppm}}.
  }
  \item{prior.mean}{
    Optional vector of prior means for canonical parameters (for
       \code{method="VBlogi"}). See Details.
  }
  \item{prior.var}{
    Optional prior variance covariance matrix for canonical parameters (for \code{method="VBlogi"}). See Details.
  }
  \item{nd}{
    Optional. Integer or pair of integers.
    The dimension of the grid of dummy points (\code{nd * nd}
    or \code{nd[1] * nd[2]})
    used to evaluate the integral in the pseudolikelihood.
    Incompatible with \code{eps}.
  }
  \item{eps}{
    Optional. 
    A positive number, or a vector of two positive numbers, giving the
    horizontal and vertical spacing, respectively, of the grid of
    dummy points. Incompatible with \code{nd}.
  }
  \item{quad.args}{
    Arguments controlling the construction of the quadrature scheme,
    when \code{Q} is a point pattern. A list of arguments that will be
    passed to \code{\link[spatstat.geom]{quadscheme}}
    or (if \code{method="logi"}) to \code{\link[spatstat.geom]{quadscheme.logi}}.
  }
  \item{gcontrol}{
    Optional. List of parameters passed to \code{\link{glm.control}}
    (or passed to \code{\link[mgcv]{gam.control}} if \code{use.gam=TRUE})
    controlling the model-fitting algorithm. 
  }
  \item{nsim}{
    Number of simulated realisations
    to generate (for \code{improve.type="ho"})
  }
  \item{nrmh}{
    Number of Metropolis-Hastings iterations
    for each simulated realisation (for \code{improve.type="ho"})
  }
  \item{start,control}{
    Arguments passed to \code{\link[spatstat.random]{rmh}} controlling the behaviour
    of the Metropolis-Hastings algorithm (for \code{improve.type="ho"})
  }
  \item{verb}{
    Logical flag indicating whether to print progress reports
    (for \code{improve.type="ho"})
  }
  \item{callstring}{
    Internal use only.
  }
}
\value{
  An object of class \code{"ppm"} describing a fitted point process
  model.
 
  See \code{\link{ppm.object}} for details of the format of this object
  and methods available for manipulating it.
}
\details{
  \bold{NOTE:} This help page describes the \bold{old syntax} of the
  function \code{ppm}, described in many older documents.
  This old syntax is still supported. However, if you are learning about
  \code{ppm} for the first time, we recommend you use the
  \bold{new syntax} described in the help file for \code{\link{ppm}}.
  
  This function fits a point process model
  to an observed point pattern.
  The model may include
  spatial trend, interpoint interaction, and dependence on covariates.

  \describe{
    \item{basic use:}{
      In basic use, \code{Q} is a point pattern dataset
      (an object of class \code{"ppp"}) to which we wish to fit a model.

      The syntax of \code{ppm()} is closely analogous to the \R functions
      \code{\link{glm}} and \code{\link[mgcv]{gam}}.
      The analogy is:
      \tabular{ll}{
	\bold{glm} \tab \bold{ppm} \cr
	\code{formula} \tab \code{trend} \cr
	\code{family} \tab \code{interaction}
      }
      The point process model to be fitted is specified by the 
      arguments \code{trend} and \code{interaction}
      which are respectively analogous to
      the \code{formula} and \code{family} arguments of glm(). 
 
      Systematic effects (spatial trend and/or dependence on 
      spatial covariates) are specified by the argument
      \code{trend}. This is an \R formula object, which may be expressed
      in terms of the Cartesian coordinates \code{x}, \code{y},
      the marks \code{marks},
      or the variables in \code{covariates} (if supplied), or both.
      It specifies the \bold{logarithm} of the first order potential
      of the process.
      The formula should not use any names beginning with \code{.mpl}
      as these are reserved for internal use.
      If \code{trend} is absent or equal to the default, \code{~1}, then
      the model to be fitted is stationary (or at least, its first order 
      potential is constant). 
 
      The symbol \code{.} in the trend expression stands for
      all the covariates supplied in the argument \code{data}.
      For example the formula \code{~ .} indicates an additive
      model with a main effect for each covariate in \code{data}.
  
      Stochastic interactions between random points of the point process
      are defined by the argument \code{interaction}. This is an object of
      class \code{"interact"} which is initialised in a very similar way to the
      usage of family objects in \code{\link{glm}} and \code{\link[mgcv]{gam}}.
      The models currently available are:
      \GibbsInteractionsList.
      See the examples below.
      It is also possible to combine several interactions
      using \code{\link{Hybrid}}.
 
      If \code{interaction} is missing or \code{NULL},
      then the model to be fitted
      has no interpoint interactions, that is, it is a Poisson process
      (stationary or nonstationary according to \code{trend}). In this case
      the methods of maximum pseudolikelihood and maximum logistic likelihood
      coincide with maximum likelihood. 

      The fitted point process model returned by this function can be printed 
      (by the print method \code{\link{print.ppm}})
      to inspect the fitted parameter values.
      If a nonparametric spatial trend was fitted, this can be extracted using
      the predict method \code{\link{predict.ppm}}.
    }
    \item{Models with covariates:}{
      To fit a model involving spatial covariates
      other than the Cartesian coordinates \eqn{x} and \eqn{y},
      the values of the covariates should be supplied in the
      argument \code{covariates}. 
      Note that it is not sufficient to have observed
      the covariate only at the points of the data point pattern; 
      the covariate must also have been observed at other 
      locations in the window.

      Typically the argument \code{covariates} is a list,
      with names corresponding to variables in the \code{trend} formula.
      Each entry in the list is either
      \describe{
	\item{a pixel image,}{
	  giving the values of a spatial covariate at 
	  a fine grid of locations. It should be an object of
	  class \code{"im"}, see \code{\link[spatstat.geom]{im.object}}.
	}
	\item{a function,}{
	  which can be evaluated
	  at any location \code{(x,y)} to obtain the value of the spatial
	  covariate. It should be a \code{function(x, y)}
	  or \code{function(x, y, ...)} in the \R language.
	  For marked point pattern data, the covariate can be a
	  \code{function(x, y, marks)} or \code{function(x, y, marks, ...)}.
	  The first two arguments of the function should be the
	  Cartesian coordinates \eqn{x} and \eqn{y}. The function may have
	  additional arguments; if the function does not have default
	  values for these additional arguments, then the user must
	  supply values for them, in \code{covfunargs}.
	  See the Examples.
	}
	\item{a window,}{
	  interpreted as a logical variable
	  which is \code{TRUE} inside the window and \code{FALSE} outside
	  it. This should be an object of class \code{"owin"}.
	}
	\item{a tessellation,}{
	  interpreted as a factor covariate.
	  For each spatial location, the factor value indicates
	  which tile of the tessellation it belongs to.
	  This should be an object of class \code{"tess"}.
	}
	\item{a single number,}{indicating a covariate that is
	  constant in this dataset.
	}
      }
      The software will look up
      the values of each covariate at the required locations
      (quadrature points).

      Note that, for covariate functions, only the \emph{name} of the
      function appears in the trend formula. A covariate function is
      treated as if it were a single variable. The function arguments do not
      appear in the trend formula. See the Examples.

      If \code{covariates} is a list,
      the list entries should have names corresponding to
      the names of covariates in the model formula \code{trend}.
      The variable names \code{x}, \code{y} and \code{marks}
      are reserved for the Cartesian 
      coordinates and the mark values,
      and these should not be used for variables in \code{covariates}.

      If \code{covariates} is a data frame, \code{Q} must be a
      quadrature scheme (see under Quadrature Schemes below).
      Then \code{covariates} must have
      as many rows as there are points in \code{Q}.
      The \eqn{i}th row of \code{covariates} should contain the values of
      spatial variables which have been observed
      at the \eqn{i}th point of \code{Q}. 
    }
    \item{Quadrature schemes:}{
      In advanced use, \code{Q} may be a `quadrature scheme'.
      This was originally just a technicality but it has turned out
      to have practical uses, as we explain below.

      Quadrature schemes are required for our implementation of
      the method of maximum pseudolikelihood.
      The definition of the pseudolikelihood involves an integral over
      the spatial window containing the data. In practice this integral
      must be approximated by a finite sum over a set of quadrature points.
      We use the technique of Baddeley and Turner (2000), a generalisation
      of the Berman-Turner (1992) device. In this technique the quadrature
      points for the numerical approximation include all the data points
      (points of the observed point pattern) as well as
      additional `dummy' points. 

      Quadrature schemes are also required for 
      the method of maximum logistic likelihood, which
      combines the data points with additional `dummy' points.
      
      A quadrature scheme is an object of class \code{"quad"}
      (see \code{\link[spatstat.geom]{quad.object}})
      which specifies both the data point pattern and the dummy points
      for the quadrature scheme, as well as the quadrature weights
      associated with these points.
      If \code{Q} is simply a point pattern
      (of class \code{"ppp"}, see \code{\link[spatstat.geom]{ppp.object}})
      then it is interpreted as specifying the
      data points only; a set of dummy points specified
      by \code{\link[spatstat.geom]{default.dummy}()} is added,
      and the default weighting rule is
      invoked to compute the quadrature weights.
 
      Finer quadrature schemes (i.e. those with more dummy
      points) generally yield a better approximation, at the
      expense of higher computational load. 

      An easy way to fit models using a finer quadrature scheme
      is to let \code{Q} be the original point pattern data,
      and use the argument \code{nd}
      to determine the number of dummy points in the quadrature scheme.
      
      Complete control over the quadrature scheme is possible.
      See \code{\link[spatstat.geom]{quadscheme}} for an overview.
      Use \code{quadscheme(X, D, method="dirichlet")} to compute
      quadrature weights based on the Dirichlet tessellation,
      or \code{quadscheme(X, D, method="grid")} to compute
      quadrature weights by counting points in grid squares,
      where \code{X} and \code{D} are the patterns of data points
      and dummy points respectively.
      Alternatively use \code{\link[spatstat.geom]{pixelquad}} to make a quadrature
      scheme with a dummy point at every pixel in a pixel image.

      The argument \code{quad.args} can be used to control the
      construction of the quadrature scheme. For example
      \code{quad.args=list(quasi=TRUE, method="dirichlet", eps=0.1)} would
      create dummy points according to a quasirandom pattern, with a
      typical spacing of 0.1 units between dummy points,
      and compute quadrature weights based on the Dirichlet tessellation.

      A practical advantage of quadrature schemes arises when we want to fit
      a model involving covariates (e.g. soil pH). Suppose we have only been
      able to observe the covariates at a small number of locations.
      Suppose \code{cov.dat} is a data frame containing the values of
      the covariates at the data points (i.e.\ \code{cov.dat[i,]}
      contains the observations for the \code{i}th data point)
      and \code{cov.dum} is another data frame (with the same columns as
      \code{cov.dat}) containing the covariate values at another
      set of points whose locations are given by the point pattern \code{Y}.
      Then setting \code{Q = quadscheme(X,Y)} combines the data points
      and dummy points into a quadrature scheme, and 
      \code{covariates = rbind(cov.dat, cov.dum)} combines the covariate
      data frames. We can then fit the model by calling
      \code{ppm(Q, ..., covariates)}.
    }
    \item{Model-fitting technique:}{
      There are several choices for the technique used
      to fit the model.
      \describe{
	\item{method="mpl"}{
	  (the default):
	  the model will be fitted by maximising the 
	  pseudolikelihood (Besag, 1975) using the
	  Berman-Turner computational approximation
	  (Berman and Turner, 1992; Baddeley and Turner, 2000).
	  Maximum pseudolikelihood is equivalent to maximum likelihood
	  if the model is a Poisson process. 
	  Maximum pseudolikelihood is biased if the
	  interpoint interaction is very strong, unless there
	  is a large number of dummy points.
	  The default settings for \code{method='mpl'}
	  specify a moderately large number of dummy points,
	  striking a compromise between speed and accuracy.
	}
	\item{method="logi":}{
	  the model will be fitted by maximising the 
	  logistic likelihood (Baddeley et al, 2014).
	  This technique is roughly equivalent in speed to
	  maximum pseudolikelihood, but is 
	  believed to be less biased. Because it is less biased,
	  the default settings for \code{method='logi'}
	  specify a relatively small number of dummy points,
	  so that this method is the fastest, in practice.
	}
	\item{method="VBlogi":}{	  
	  the model will be fitted in a Bayesian setup by maximising the
	  posterior probability density for the canonical model
	  parameters. This uses the variational Bayes approximation to
	  the posterior derived from the logistic likelihood as described
	  in Rajala (2014). The prior is assumed to be multivariate
	  Gaussian with mean vector \code{prior.mean} and variance-covariance
	  matrix \code{prior.var}.	  
	}
      }
      Note that \code{method='logi'} and \code{method='VBlogi'}
      involve randomisation, so that the results are
      subject to random variation.

      After this initial fit, there are several ways to improve the fit:
      \describe{
	\item{improve.type="none":}{
	  No further improvement is performed.
	}
	\item{improve.type="ho":}{
	  the model will be re-fitted
	  by applying the approximate maximum likelihood
	  method of Huang and Ogata (1999). See below.
	  The Huang-Ogata method is slower than the other options,
	  but has better statistical properties.
	  This method involves randomisation, so the results are subject
	  to random variation.
	}
	\item{improve.type="enet":}{
	  The model will be re-fitted using a regularized version of the
	  composite likelihood. See below.
	}
      }
    }
    \item{Huang-Ogata method:}{
      If \code{improve.type="ho"} then the model will be fitted using
      the Huang-Ogata (1999) approximate maximum likelihood method.
      First the model is fitted by maximum pseudolikelihood as
      described above, yielding an initial estimate of the parameter
      vector \eqn{\theta_0}{theta0}.
      From this initial model, \code{nsim} simulated
      realisations are generated. The score and Fisher information of
      the model at \eqn{\theta=\theta_0}{theta=theta0}
      are estimated from the simulated realisations. Then one step
      of the Fisher scoring algorithm is taken, yielding an updated
      estimate \eqn{\theta_1}{theta1}. The corresponding model is
      returned.

      Simulated realisations are generated using \code{\link[spatstat.random]{rmh}}.
      The iterative behaviour of the Metropolis-Hastings algorithm
      is controlled by the arguments \code{start} and \code{control}
      which are passed to \code{\link[spatstat.random]{rmh}}.

      As a shortcut, the argument
      \code{nrmh} determines the number of Metropolis-Hastings
      iterations run to produce one simulated realisation (if
      \code{control} is absent). Also
      if \code{start} is absent or equal to \code{NULL}, it defaults to
      \code{list(n.start=N)} where \code{N} is the number of points
      in the data point pattern.
    }
    \item{Regularization:}{
      This requires the package \pkg{glmnet}. 
      \bold{Details to be written.}
   }
    \item{Edge correction}{
      Edge correction should be applied to the sufficient statistics
      of the model, to reduce bias.
      The argument \code{correction} is the name of an edge correction
      method.
      The default \code{correction="border"} specifies the border correction,
      in which the quadrature window (the domain of integration of the 
      pseudolikelihood) is obtained by trimming off a margin of width
      \code{rbord} from the observation window of the data pattern.
      Not all edge corrections are implemented (or implementable)
      for arbitrary windows.
      Other options depend on the argument \code{interaction}, but these
      generally include \code{correction="periodic"} (the periodic or toroidal edge
      correction in which opposite edges of a rectangular window are
      identified) and \code{correction="translate"} (the translation correction,
      see Baddeley 1998 and Baddeley and Turner 2000).
      For pairwise interaction models
      there is also Ripley's isotropic correction,
      identified by \code{correction="isotropic"} or \code{"Ripley"}.
    }
    \item{Subsetting}{
      The arguments \code{subset} and \code{clipwin} specify that the
      model should be fitted to a restricted subset of the available
      data. These arguments are equivalent for Poisson point process models,
      but different for Gibbs models.
      If \code{clipwin} is specified, then all the available data will
      be restricted to this spatial region, and data outside this region
      will be discarded, before the model is fitted.
      If \code{subset} is specified, then no data are deleted, but
      the domain of integration of the likelihood or pseudolikelihood
      is restricted to the \code{subset}.
      For Poisson models, these two arguments have the same effect;
      but for a Gibbs model, 
      interactions between points inside and outside the \code{subset}
      are taken into account, while
      interactions between points inside and outside the \code{clipwin}
      are ignored.
    }
  }
}
\section{Interaction parameters}{
  Apart from the Poisson model, every point process model fitted by
  \code{ppm} has parameters that determine the strength and
  range of \sQuote{interaction} or dependence between points.
  These parameters are of two types:
  \describe{
    \item{regular parameters:}{
      A parameter \eqn{\phi}{phi} is called \emph{regular}
      if the log likelihood is a linear function of \eqn{\theta}{theta} where 
      \eqn{\theta = \theta(\psi)}{theta = theta(psi)} is some transformation of 
      \eqn{\psi}{psi}. [Then \eqn{\theta}{theta} is called the canonical
      parameter.]
    }
    \item{irregular parameters}{
      Other parameters are called \emph{irregular}. 
    }
  }
  Typically, regular parameters determine the \sQuote{strength}
  of the interaction, while irregular parameters determine the
  \sQuote{range} of the interaction. For example, the Strauss process
  has a regular parameter \eqn{\gamma}{gamma} controlling the strength
  of interpoint inhibition, and an irregular parameter \eqn{r}
  determining the range of interaction.

  The \code{ppm} command is only designed to estimate regular
  parameters of the interaction.
  It requires the values of any irregular parameters of the interaction
  to be fixed. For example, to fit a Strauss process model to the \code{cells}
  dataset, you could type \code{ppm(cells, ~1, Strauss(r=0.07))}.
  Note that the value of the irregular parameter \code{r} must be given.
  The result of this command will be a fitted model in which the
  regular parameter \eqn{\gamma}{gamma} has been estimated.

  To determine the irregular parameters, there are several
  practical techniques, but no general statistical theory available.
  Useful techniques include maximum profile pseudolikelihood, which
  is implemented in the command \code{\link{profilepl}},
  and Newton-Raphson maximisation, implemented in the
  experimental command \code{\link{ippm}}. 

  Some irregular parameters can be estimated directly from data:
  the hard-core radius in the model \code{\link{Hardcore}}
  and the matrix of hard-core radii in \code{\link{MultiHard}} can be
  estimated easily from data. In these cases, \code{ppm} allows the user
  to specify the interaction without giving
  the value of the irregular parameter. The user can give the
  hard core interaction as \code{interaction=Hardcore()}
  or even \code{interaction=Hardcore}, and 
  the hard core radius will then be estimated from the data.
}
\references{
  Baddeley, A., Coeurjolly, J.-F., Rubak, E. and Waagepetersen, R. (2014)
  Logistic regression for spatial Gibbs point processes.
  \emph{Biometrika} \bold{101} (2) 377--392.

  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Berman, M. and Turner, T.R. 
  Approximating point process likelihoods with GLIM.
  \emph{Applied Statistics} \bold{41} (1992) 31--38.
 
  Besag, J.
  Statistical analysis of non-lattice data.
  \emph{The Statistician} \bold{24} (1975) 179-195.
 
  Diggle, P.J., Fiksel, T., Grabarnik, P., Ogata, Y., Stoyan, D. and
  Tanemura, M.
  On parameter estimation for pairwise interaction processes.
  \emph{International Statistical Review} \bold{62} (1994) 99-117.

  Huang, F. and Ogata, Y.
  Improvements of the maximum pseudo-likelihood estimators
  in various spatial statistical models.
  \emph{Journal of Computational and Graphical Statistics}
  \bold{8} (1999) 510-530.
  
  Jensen, J.L. and Moeller, M.
  Pseudolikelihood for exponential family models of spatial point processes.
  \emph{Annals of Applied Probability} \bold{1} (1991) 445--461.
 
  Jensen, J.L. and Kuensch, H.R. 
  On asymptotic normality of pseudo likelihood
  estimates for pairwise interaction processes,
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{46} (1994) 475-486.

  Rajala T. (2014)
  \emph{A note on Bayesian logistic regression for spatial exponential family
  Gibbs point processes},
  Preprint on ArXiv.org. \url{https://arxiv.org/abs/1411.0539}

}

\seealso{
  \code{\link{ppm.object}} for details of how to
  print, plot and manipulate a fitted model.

  \code{\link[spatstat.geom]{ppp}} and \code{\link[spatstat.geom]{quadscheme}}
  for constructing data.
  
  Interactions:
  \GibbsInteractionsList.

  See \code{\link{profilepl}} for advice on
  fitting nuisance parameters in the interaction,
  and \code{\link{ippm}} for irregular parameters in the trend.

  See \code{\link{valid.ppm}} and \code{\link{emend.ppm}} for
  ensuring the fitted model is a valid point process.
}
\section{Error and Warning Messages}{
  Some common error messages and warning messages 
  are listed below, with explanations.
  \describe{
    \item{\dQuote{Model is invalid} or \dQuote{Model is not valid}}{
      The fitted model coefficients do not define a valid point
      process. This can occur because some of the fitted coefficients
      are \code{NA} (perhaps because the model formula included redundant
      covariates so that the coefficients cannot be estimated), or
      because the fitted interaction coefficients do not define
      a valid point process (e.g. because a point process
      model which always has inhibition between points was fitted to a
      clustered point pattern). See \code{\link{valid.ppm}} for detailed
      information. 
    }
    \item{\dQuote{System is computationally singular} or
      \dQuote{Fisher information matrix is singular}}{
      The software cannot calculate standard errors or confidence
      intervals for the coefficients of the fitted model.
      This requires the (asymptotic) variance-covariance
      matrix, which is the inverse matrix of the Fisher information
      matrix of the fitted model. The error message states
      that the determinant of the Fisher information matrix is zero,
      or close to zero, so that the matrix cannot be inverted.
      This error is usually reported when the model is printed,
      because the \code{print} method calculates standard errors for the
      fitted parameters. Singularity usually occurs because the spatial
      coordinates in the original data were very large numbers
      (e.g. expressed in metres) so that the fitted coefficients were
      very small numbers. The simple remedy is to
      \bold{rescale the data}, for example, to convert from metres to
      kilometres by \code{X <- \link[spatstat.geom]{rescale}(X, 1000)}, then re-fit the
      model. Singularity can also occur if the covariate values are
      very large numbers, or if the covariates are approximately
      collinear.       
    }
    \item{\dQuote{Covariate values were NA or undefined at X\%
	(M out of N) of the quadrature points}}{
      The covariate data (typically a pixel image) did not provide
      values of the covariate at some of the spatial locations in the 
      observation window of the point pattern. This means that the
      spatial domain of the pixel image does not completely cover the
      observation window of the point pattern. If the percentage
      is small, this warning can be ignored - typically it happens
      because of rounding effects which cause the pixel image
      to be one-pixel-width narrower than the observation window.
      However if more than a few percent of covariate values are
      undefined, it would be prudent to check that the pixel images
      are correct, and are correctly registered in their spatial relation to
      the observation window.
    }
    \item{\dQuote{Some tiles with positive area do not contain any
	quadrature points: relative error = X\%}}{
      A problem has arisen when creating the quadrature scheme
      used to fit the model. In the default rule for computing the
      quadrature weights, space is divided into rectangular tiles,
      and the number of quadrature points (data and dummy points) in
      each tile is counted. It is possible for a tile with non-zero area
      to contain no quadrature points; in this case, the quadrature
      scheme will contribute a bias to the model-fitting procedure.
      \bold{A small relative error (less than 2 percent) is not	important.}
      Relative errors of a few percent can occur because of the shape of
      the window.
      If the relative error is greater than about 5 percent, we
      recommend trying different parameters for the quadrature scheme,
      perhaps setting a larger value of \code{nd} to increase the number
      of dummy points. A relative error greater than 10 percent
      indicates a major problem with the input data: in this case,
      extract the quadrature scheme by applying \code{\link{quad.ppm}}
      to the fitted model, and inspect it. 
      (The most likely cause of this problem is that the spatial coordinates
      of the original data were not handled correctly, for example,
      coordinates of the locations and the window boundary were incompatible.)
    }
    \item{\dQuote{Model is unidentifiable}}{
      It is not possible to estimate all the model parameters
      from this dataset. The error message gives a further explanation,
      such as \dQuote{data pattern is empty}.
      Choose a simpler model, or check the data.
    }
    \item{\dQuote{N data points are illegal (zero conditional intensity)}}{
      In a Gibbs model (i.e. with interaction between
      points), the conditional intensity may be zero at some spatial
      locations, indicating that the model forbids the presence of a 
      point at these locations. However if the conditional intensity is
      zero \emph{at a data point}, this means that the model is
      inconsistent with the data. Modify the interaction parameters so
      that the data point is not illegal (e.g. reduce the value of the
      hard core radius) or choose a different interaction.
    }
  }
}
\section{Warnings}{
  The implementation of the Huang-Ogata method is experimental;
  several bugs were fixed in \pkg{spatstat} 1.19-0.
  
  See the comments above about the possible inefficiency
  and bias of the maximum pseudolikelihood estimator.
 
  The accuracy of the Berman-Turner approximation to
  the pseudolikelihood depends on the number of dummy points used
  in the quadrature scheme. The number of dummy points should 
  at least equal the number of data points.
 
  The parameter values of the fitted model
  do not necessarily determine a valid point process.
  Some of the point process models are only defined when the parameter
  values lie in a certain subset. For example the Strauss process only 
  exists when the interaction parameter \eqn{\gamma}{gamma}
  is less than or equal to \eqn{1},
  corresponding to a value of \code{ppm()$theta[2]}
  less than or equal to \code{0}.

  By default (if \code{emend=FALSE}) the algorithm
  maximises the pseudolikelihood
  without constraining the parameters, and does not apply any checks for
  sanity after fitting the model.
  This is because the fitted parameter value
  could be useful information for data analysis.
  To constrain the parameters to ensure that the model is a valid
  point process, set \code{emend=TRUE}. See also the functions
  \code{\link{valid.ppm}} and \code{\link{emend.ppm}}.
  
  The \code{trend} formula should not use any variable names
  beginning with the prefixes \code{.mpl} or \code{Interaction}
  as these names are reserved
  for internal use. The data frame \code{covariates} should have as many rows
  as there are points in \code{Q}. It should not contain
  variables called \code{x}, \code{y} or \code{marks}
  as these names are reserved for the Cartesian coordinates
  and the marks.
 
  If the model formula involves one of the functions
  \code{poly()}, \code{bs()}
  or \code{ns()}
  (e.g. applied to spatial coordinates \code{x} and \code{y}),
  the fitted coefficients can be misleading.
  The resulting fit is not to the raw spatial variates
  (\code{x}, \code{x^2}, \code{x*y}, etc.) 
  but to a transformation of these variates.  The transformation is implemented
  by \code{poly()} in order to achieve better numerical stability.
  However the
  resulting coefficients are appropriate for use with the transformed
  variates, not with the raw variates.  
  This affects the interpretation of the constant
  term in the fitted model, \code{logbeta}. 
  Conventionally, \eqn{\beta}{beta} is the background intensity, i.e. the  
  value taken by the conditional intensity function when all predictors
  (including spatial or ``trend'' predictors) are set equal to \eqn{0}.
  However the coefficient actually produced is the value that the
  log conditional intensity takes when all the predictors, 
  including the \emph{transformed}
  spatial predictors, are set equal to \code{0}, which is not the same thing.

  Worse still, the result of \code{\link{predict.ppm}} can be
  completely wrong if the trend formula contains one of the
  functions \code{poly()}, \code{bs()}
  or \code{ns()}. This is a weakness of the underlying
  function \code{\link{predict.glm}}. 

  If you wish to fit a polynomial trend, 
  we offer an alternative to \code{\link[stats]{poly}()},
  namely \code{polynom()}, which avoids the
  difficulty induced by transformations.  It is completely analogous
  to \code{poly} except that it does not orthonormalise.
  The resulting coefficient estimates then have
  their natural interpretation and can be predicted correctly. 
  Numerical stability may be compromised.

  Values of the maximised pseudolikelihood are not comparable
  if they have been obtained with different values of \code{rbord}.
}
\examples{
 # fit the stationary Poisson process
 # to point pattern 'nztrees'

 ppm(nztrees)
 ppm(nztrees ~ 1)
 # equivalent.

 Q <- quadscheme(nztrees) 
 ppm(Q) 
 # equivalent.

 fit1 <- ppm(nztrees, ~ x)
 # fit the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(a + bx)
 # where x,y are the Cartesian coordinates
 # and a,b are parameters to be estimated

 # For other examples, see help(ppm)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{models}
