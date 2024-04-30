\name{spatstat.model-package}
\alias{spatstat.model-package} 
\alias{spatstat.model} 
\docType{package}
\title{The spatstat.model Package}
\description{
  The \pkg{spatstat.model} package
  belongs to the \pkg{spatstat} family of packages.
  It contains the core functionality
  for parametric statistical modelling of spatial data.
}
\details{
  \pkg{spatstat} is
  a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  The original \pkg{spatstat} package
  has now been split into several
  sub-packages.

  This sub-package \pkg{spatstat.model} contains all the main
  user-level functions
  that perform parametric statistical modelling of spatial data.

  (The main exception is that functions for linear networks
  are in the separate sub-package \pkg{spatstat.linnet}.)
}
\section{Structure of the spatstat family}{  
  The orginal \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.explore} containing the functionality
    for exploratory analysis and nonparametric modelling of spatial data
    \item \pkg{spatstat.model} containing the main functionality
    for parametric modelling, analysis and inference for spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.

  For an overview of all the functions available in
  the sub-packages of \pkg{spatstat},
  see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package.

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{Overview of Functionality in \pkg{spatstat.model}}{
  The \pkg{spatstat} family of packages
  is designed to support a complete statistical analysis
  of spatial data. It supports
  \itemize{
    \item creation, manipulation and plotting of point patterns;
    \item exploratory data analysis;
    \item spatial random sampling;
    \item simulation of point process models;
    \item parametric model-fitting;
    \item non-parametric smoothing and regression;
    \item formal inference (hypothesis tests, confidence intervals);
    \item model diagnostics.
  }
  For an overview, see the help file for \code{"spatstat-package"}
  in the \pkg{spatstat} package. 

  Following is a list of the functionality provided in the
  \pkg{spatstat.model} package only.

  \bold{To simulate a random point pattern:}
  
   Functions for generating random point patterns
   are now contained in the \pkg{spatstat.random} package.

     \bold{Exploratory analysis}

     Exploratory graphics, 
     smoothing, and exploratory analysis of spatial data
     are now provided in the \pkg{spatstat.explore} package.


 \bold{Model fitting (Cox and cluster models)}
  
  Cluster process models (with homogeneous or inhomogeneous intensity)
  and Cox processes can be fitted by the function \code{\link[spatstat.model]{kppm}}.
  Its result is an object of class \code{"kppm"}.
  The fitted model can be printed, plotted, predicted, simulated
  and updated.

  \tabular{ll}{
    \code{\link[spatstat.model]{kppm}} \tab  Fit model\cr
    \code{\link[spatstat.model]{plot.kppm}} \tab  Plot the fitted model\cr
    \code{\link[spatstat.model]{summary.kppm}} \tab  Summarise the fitted model\cr
    \code{\link[spatstat.model]{fitted.kppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.model]{predict.kppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.model]{update.kppm}} \tab Update the model \cr
    \code{\link[spatstat.model]{improve.kppm}} \tab Refine the estimate of trend \cr
    \code{\link[spatstat.model]{simulate.kppm}} \tab Generate simulated realisations \cr
    \code{\link[spatstat.model]{vcov.kppm}} \tab Variance-covariance matrix of coefficients \cr
    \code{\link[spatstat.model:methods.kppm]{coef.kppm}}
    \tab Extract trend coefficients \cr
    \code{\link[spatstat.model:methods.kppm]{formula.kppm}}
    \tab Extract trend formula \cr
    \code{\link[spatstat.model]{parameters}} \tab Extract all model parameters \cr
    \code{\link[spatstat.model]{clusterfield.kppm}} \tab Compute offspring density \cr
    \code{\link[spatstat.model]{clusterradius.kppm}} \tab Radius of support of offspring density \cr
    \code{\link[spatstat.model]{Kmodel.kppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link[spatstat.model]{pcfmodel.kppm}} \tab Pair correlation of fitted model 
  }
  
  For model selection, you can also use 
  the generic functions \code{\link[stats]{step}}, \code{\link[stats]{drop1}} 
  and \code{\link[stats]{AIC}} on fitted point process models.
  For variable selection, see \code{\link[spatstat.explore]{sdr}}.
  
  The theoretical models can also be simulated,
  for any choice of parameter values,
  using \code{\link[spatstat.random]{rThomas}}, \code{\link[spatstat.random]{rMatClust}},
  \code{\link[spatstat.random]{rCauchy}}, \code{\link[spatstat.random]{rVarGamma}},
  and \code{\link[spatstat.random]{rLGCP}}.
  
  Lower-level fitting functions include:

  \tabular{ll}{
    \code{\link[spatstat.model]{lgcp.estK}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link[spatstat.model]{lgcp.estpcf}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link[spatstat.model]{thomas.estK}} \tab fit the Thomas process model \cr
    \code{\link[spatstat.model]{thomas.estpcf}} \tab fit the Thomas process model \cr
    \code{\link[spatstat.model]{matclust.estK}} \tab fit the \Matern Cluster process model \cr
    \code{\link[spatstat.model]{matclust.estpcf}} \tab fit the \Matern Cluster process model \cr
    \code{\link[spatstat.model]{cauchy.estK}} \tab fit a Neyman-Scott Cauchy cluster process \cr
   \code{\link[spatstat.model]{cauchy.estpcf}} \tab fit a Neyman-Scott Cauchy cluster process\cr
   \code{\link[spatstat.model]{vargamma.estK}} \tab fit a Neyman-Scott Variance Gamma process\cr
   \code{\link[spatstat.model]{vargamma.estpcf}} \tab fit a Neyman-Scott Variance Gamma process\cr
    \code{\link[spatstat.model]{mincontrast}} \tab low-level algorithm for fitting models
    \cr \tab by the method of minimum contrast 
  }

  \bold{Model fitting (Poisson and Gibbs models)}
  
  Poisson point processes are the simplest models for point patterns.
  A Poisson model assumes that the points are stochastically
  independent. It may allow the points to have a non-uniform spatial
  density. The special case of a Poisson process with a uniform
  spatial density is often called Complete Spatial Randomness.
  
  Poisson point processes are included in the more general class of Gibbs point
  process models. In a Gibbs model, there is \emph{interaction}
  or dependence between points. Many different types of interaction
  can be specified.
  
  For a detailed explanation of how to fit Poisson or
  Gibbs point process models to point pattern data using \pkg{spatstat},
  see Baddeley and Turner (2005b) or Baddeley (2008).
  
  \bold{To fit a Poisson or Gibbs point process model:}

  Model fitting in \pkg{spatstat} is performed mainly by the function
  \code{\link[spatstat.model]{ppm}}. Its result is an object of class \code{"ppm"}.
  
  Here are some examples, where \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{ppm(X)} \tab Complete Spatial Randomness \cr
    \code{ppm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{ppm(X ~ x)} \tab Poisson process with \cr
                             \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{ppm(X ~ 1, Strauss(0.1))} \tab Stationary Strauss process \cr
    \code{ppm(X ~ x, Strauss(0.1))} \tab Strauss process with \cr
                              \tab conditional intensity loglinear in \eqn{x}
  }
  It is also possible to fit models that depend on
  other covariates.

  \bold{Manipulating the fitted model:}

  \tabular{ll}{
    \code{\link[spatstat.model]{plot.ppm}} \tab 		Plot the fitted model\cr
    \code{\link[spatstat.model]{predict.ppm}}
    \tab   Compute the spatial trend and conditional intensity\cr
    \tab   of the fitted point process model \cr
    \code{\link[spatstat.model]{coef.ppm}} \tab Extract the fitted model coefficients\cr
    \code{\link[spatstat.model]{parameters}} \tab Extract all model parameters\cr
    \code{\link[spatstat.model]{formula.ppm}} \tab Extract the trend formula\cr
    \code{\link[spatstat.model]{intensity.ppm}} \tab Compute fitted intensity \cr
    \code{\link[spatstat.model]{Kmodel.ppm}} \tab \eqn{K} function of fitted model \cr
    \code{\link[spatstat.model]{pcfmodel.ppm}} \tab pair correlation of fitted model \cr
    \code{\link[spatstat.model]{fitted.ppm}} \tab Compute fitted conditional intensity at quadrature points \cr
    \code{\link[spatstat.model]{residuals.ppm}} \tab Compute point process residuals at quadrature points \cr
    \code{\link[spatstat.model]{update.ppm}} \tab Update the fit \cr
    \code{\link[spatstat.model]{vcov.ppm}} \tab Variance-covariance matrix of estimates\cr
    \code{\link[spatstat.model]{rmh.ppm}} \tab Simulate from fitted model  \cr
    \code{\link[spatstat.model]{simulate.ppm}} \tab Simulate from fitted model  \cr
    \code{\link[spatstat.model]{print.ppm}} \tab Print basic information about a fitted model\cr
    \code{\link[spatstat.model]{summary.ppm}} \tab Summarise a fitted model\cr
    \code{\link[spatstat.model]{effectfun}} \tab Compute the fitted effect of one covariate\cr
    \code{\link[spatstat.model]{logLik.ppm}} \tab log-likelihood or log-pseudolikelihood\cr
    \code{\link[spatstat.model]{anova.ppm}} \tab Analysis of deviance \cr
    \code{\link[spatstat.model]{model.frame.ppm}} \tab Extract data frame used to fit model  \cr
    \code{\link[spatstat.model]{model.images}} \tab Extract spatial data used to fit model  \cr
    \code{\link[spatstat.model]{model.depends}} \tab Identify variables in the model \cr
    \code{\link[spatstat.model]{as.interact}} \tab Interpoint interaction component of model \cr
    \code{\link[spatstat.model]{fitin}} \tab Extract fitted interpoint interaction \cr
    \code{\link[spatstat.model]{is.hybrid}} \tab Determine whether the model is a hybrid \cr
    \code{\link[spatstat.model]{valid.ppm}} \tab Check the model is a valid point process \cr
    \code{\link[spatstat.model]{project.ppm}} \tab Ensure the model is a valid point process 
  }
  For model selection, you can also use 
  the generic functions \code{\link[stats]{step}}, \code{\link[stats]{drop1}} 
  and \code{\link[stats]{AIC}} on fitted point process models.
  For variable selection, see \code{\link[spatstat.explore]{sdr}}.
  
  See \code{\link[spatstat.geom]{spatstat.options}} to control plotting of fitted model.
  
  \bold{To specify a point process model:}
  
  The first order ``trend'' of the model is determined by an \R 
  language formula. The formula specifies the form of the
  \emph{logarithm} of the trend.
  
  \tabular{ll}{
    \code{X ~ 1}  \tab No trend (stationary) \cr
    \code{X ~ x}   \tab Loglinear trend
      \eqn{\lambda(x,y) = \exp(\alpha + \beta x)}{lambda(x,y) =	exp(alpha + beta * x)} \cr
    \tab where \eqn{x,y} are Cartesian coordinates \cr
    \code{X ~ polynom(x,y,3)}  \tab Log-cubic polynomial trend  \cr
    \code{X ~ harmonic(x,y,2)}  \tab Log-harmonic polynomial trend \cr
    \code{X ~ Z}   \tab Loglinear function of covariate \code{Z} \cr
    \tab  \eqn{\lambda(x,y) = \exp(\alpha + \beta Z(x,y))}{lambda(x,y) =	exp(alpha + beta * Z(x,y))} 
  }

  The higher order (``interaction'') components are described by
  an object of class \code{"interact"}. Such objects are created by:
  \tabular{ll}{
    \code{\link[spatstat.model]{Poisson}()} \tab		the Poisson point process\cr
    \code{\link[spatstat.model]{AreaInter}()}	 \tab Area-interaction process\cr
    \code{\link[spatstat.model]{BadGey}()} \tab	multiscale Geyer process\cr
    \code{\link[spatstat.model]{Concom}()} \tab	connected component interaction\cr
    \code{\link[spatstat.model]{DiggleGratton}() } \tab Diggle-Gratton potential \cr
    \code{\link[spatstat.model]{DiggleGatesStibbard}() } \tab Diggle-Gates-Stibbard potential \cr
    \code{\link[spatstat.model]{Fiksel}()}	 \tab Fiksel pairwise interaction process\cr
    \code{\link[spatstat.model]{Geyer}()}	 \tab Geyer's saturation process\cr
    \code{\link[spatstat.model]{Hardcore}()}	 \tab Hard core process\cr
    \code{\link[spatstat.model]{HierHard}()}	 \tab Hierarchical multiype hard core process\cr
    \code{\link[spatstat.model]{HierStrauss}()}	 \tab Hierarchical multiype Strauss process\cr
    \code{\link[spatstat.model]{HierStraussHard}()}	 \tab
    Hierarchical multiype Strauss-hard core process\cr
    \code{\link[spatstat.model]{Hybrid}()}	 \tab Hybrid of several interactions\cr
    \code{\link[spatstat.model]{LennardJones}() } \tab Lennard-Jones potential \cr
    \code{\link[spatstat.model]{MultiHard}()} \tab 		multitype hard core process \cr
    \code{\link[spatstat.model]{MultiStrauss}()} \tab 		multitype Strauss process \cr
    \code{\link[spatstat.model]{MultiStraussHard}()} \tab 	multitype Strauss/hard core process \cr
    \code{\link[spatstat.model]{OrdThresh}()} \tab		Ord process, threshold potential\cr
    \code{\link[spatstat.model]{Ord}()} \tab 		        Ord model, user-supplied potential \cr
    \code{\link[spatstat.model]{PairPiece}()} \tab		pairwise interaction, piecewise constant \cr
    \code{\link[spatstat.model]{Pairwise}()} \tab	pairwise interaction, user-supplied potential\cr
    \code{\link[spatstat.model]{Penttinen}()} \tab	Penttinen pairwise interaction\cr
    \code{\link[spatstat.model]{SatPiece}()} \tab	Saturated pair model, piecewise  constant potential\cr
    \code{\link[spatstat.model]{Saturated}()} \tab	Saturated pair model, user-supplied potential\cr
    \code{\link[spatstat.model]{Softcore}()} \tab pairwise interaction, soft core potential\cr
    \code{\link[spatstat.model]{Strauss}()} \tab  Strauss process \cr
    \code{\link[spatstat.model]{StraussHard}()} \tab Strauss/hard core point process \cr
    \code{\link[spatstat.model]{Triplets}()} \tab Geyer triplets process
  }
  Note that it is also possible to combine several such interactions
  using \code{\link[spatstat.model]{Hybrid}}.
  
  \bold{Simulation and goodness-of-fit for fitted models:}
  
  \tabular{ll}{
    \code{\link[spatstat.model]{rmh.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link[spatstat.model]{simulate.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link[spatstat.model]{envelope}} \tab compute simulation envelopes for a
    fitted model 
  }

  \bold{Model fitting (determinantal point process models)}

  Code for fitting \emph{determinantal point process models} has 
  recently been added to \pkg{spatstat}.

  For information, see the help file for \code{\link[spatstat.model]{dppm}}.

  \bold{Model fitting (spatial logistic regression)}
  
  Pixel-based spatial logistic regression is an alternative
  technique for analysing spatial point patterns
  that is widely used in Geographical Information Systems.
  It is approximately equivalent to fitting a Poisson point process
  model.
  
  In pixel-based logistic regression, the spatial domain is
  divided into small pixels, the presence or absence of a
  data point in each pixel is recorded, and logistic regression
  is used to model the presence/absence indicators as a function
  of any covariates.
  
  Facilities for performing spatial logistic regression are
  provided in \pkg{spatstat} for comparison purposes.
  
  \bold{Fitting a spatial logistic regression}
  
  Spatial logistic regression is performed by the function
  \code{\link[spatstat.model]{slrm}}. Its result is an object of class \code{"slrm"}.
  There are many methods for this class, including methods for
  \code{print}, \code{fitted}, \code{predict}, \code{simulate},
  \code{anova}, \code{coef}, \code{logLik}, \code{terms},
  \code{update}, \code{formula} and \code{vcov}. 
  
  For example, if \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{slrm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{slrm(X ~ x)} \tab Poisson process with \cr
                \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{slrm(X ~ Z)} \tab Poisson process with \cr
                \tab intensity loglinear in covariate \code{Z}
  }

  \bold{Manipulating a fitted spatial logistic regression}
  
  \tabular{ll}{
    \code{\link[spatstat.model]{anova.slrm}} \tab Analysis of deviance \cr
    \code{\link[spatstat.model]{coef.slrm}}  \tab Extract fitted coefficients \cr
    \code{\link[spatstat.model]{vcov.slrm}}  \tab Variance-covariance matrix of fitted coefficients \cr
    \code{\link[spatstat.model]{fitted.slrm}} \tab Compute fitted probabilities or
    intensity \cr
    \code{\link[spatstat.model]{logLik.slrm}}   \tab Evaluate loglikelihood of fitted
    model \cr
    \code{\link[spatstat.model]{plot.slrm}}    \tab Plot fitted probabilities or
    intensity \cr
    \code{\link[spatstat.model]{predict.slrm}} \tab Compute predicted probabilities or
    intensity with new data \cr
    \code{\link[spatstat.model]{simulate.slrm}} \tab Simulate model
  }
  
  There are many other undocumented methods for this class,
  including methods for \code{print}, \code{update}, \code{formula}
  and \code{terms}. Stepwise model selection is
  possible using \code{step} or \code{stepAIC}.
  For variable selection, see \code{\link[spatstat.explore]{sdr}}.

  \bold{Simulation}

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}
  Functions for random generation are now contained in the \pkg{spatstat.random} package.

  See also \code{\link[spatstat.explore]{varblock}} for estimating the variance
  of a summary statistic by block resampling, and
  \code{\link[spatstat.explore]{lohboot}} for another bootstrap technique.
  
  \bold{Fitted point process models:}

   If you have fitted a point process model to a point pattern dataset,
   the fitted model can be simulated.

   Cluster process models 
   are fitted by the function \code{\link[spatstat.model]{kppm}} yielding an
   object of class \code{"kppm"}. To generate one or more simulated
   realisations of this fitted model, use 
   \code{\link[spatstat.model]{simulate.kppm}}.

   Gibbs point process models 
   are fitted by the function \code{\link[spatstat.model]{ppm}} yielding an
   object of class \code{"ppm"}. To generate a simulated
   realisation of this fitted model, use \code{\link[spatstat.model]{rmh.ppm}}.
   To generate one or more simulated realisations of the fitted model,
   use \code{\link[spatstat.model]{simulate.ppm}}.

   \bold{Other random patterns:}
   Functions for random generation are now contained in the
   \pkg{spatstat.random} package.

   \bold{Simulation-based inference}

   Simulation-based inference including simulation envelopes
   and hypothesis tests is now supported by the package \pkg{spatstat.explore}.
   

  \bold{Sensitivity diagnostics:}

  Classical measures of model sensitivity such as leverage and influence
  have been adapted to point process models.
  
    \tabular{ll}{
      \code{\link[spatstat.model]{leverage.ppm}} \tab Leverage for point process model\cr
      \code{\link[spatstat.model]{influence.ppm}} \tab Influence for point process model\cr
      \code{\link[spatstat.model]{dfbetas.ppm}} \tab Parameter influence\cr
      \code{\link[spatstat.model]{dffit.ppm}} \tab Effect change diagnostic
    }
  
  \bold{Diagnostics for covariate effect:}

  Classical diagnostics for covariate effects have been adapted to
  point process models.

  \tabular{ll}{
    \code{\link[spatstat.model]{parres}} \tab Partial residual plot\cr
    \code{\link[spatstat.model]{addvar}} \tab Added variable plot \cr
    \code{\link[spatstat.model]{rhohat.ppm}} \tab Kernel estimate of covariate effect\cr
    \code{\link[spatstat.explore]{rho2hat}} \tab Kernel estimate of covariate effect
    (bivariate) 
  }
  
  \bold{Residual diagnostics:}
  
  Residuals for a fitted point process model, and diagnostic plots
  based on the residuals, were introduced in Baddeley et al (2005) and
  Baddeley, Rubak and \Moller (2011).
  
  Type \code{demo(diagnose)}
  for a demonstration of the diagnostics features.

  \tabular{ll}{
    \code{\link[spatstat.model]{diagnose.ppm}} \tab diagnostic plots for spatial trend\cr
    \code{\link[spatstat.model]{qqplot.ppm}} \tab diagnostic Q-Q plot for interpoint interaction\cr
    \code{\link[spatstat.data]{residualspaper}} \tab examples from Baddeley et al (2005) \cr
    \code{\link[spatstat.model]{Kcom}} \tab model compensator of \eqn{K} function \cr
    \code{\link[spatstat.model]{Gcom}} \tab model compensator of \eqn{G} function \cr
    \code{\link[spatstat.model]{Kres}} \tab score residual of \eqn{K} function \cr
    \code{\link[spatstat.model]{Gres}} \tab score residual of \eqn{G} function \cr
    \code{\link[spatstat.model]{psst}} \tab pseudoscore residual of summary function \cr
    \code{\link[spatstat.model]{psstA}} \tab pseudoscore residual of empty space function \cr
    \code{\link[spatstat.model]{psstG}} \tab pseudoscore residual of \eqn{G} function \cr
    \code{\link[spatstat.model]{compareFit}} \tab compare compensators of several fitted
    models
  }


  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation and resampling
  using the following capabilities:
  
  \tabular{ll}{
    \code{\link[spatstat.random]{quadratresample}} \tab block resampling \cr
    \code{\link[spatstat.random]{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link[spatstat.random]{rthin}} \tab  random thinning  
  }

}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie,
  Tilman Davies,
  Julian Gilbey, 
  Yongtao Guan,
  Ute Hahn,
  Kassel Hingee,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Greg McSwiggan,
  Tuomas Rajala,
  Suman Rakshit,
  Dominic Schuhmacher,
  Rasmus Waagepetersen
  and
  Hangsheng Wang
  made substantial contributions of code.

  For comments, corrections, bug alerts and suggestions, we thank
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia} Cobo Sanchez,
  Jean-Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \ifelse{latex}{\out{Fr\'{e}d\'{e}ric}}{Frederic} Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

