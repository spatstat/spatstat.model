\name{spatstat.model-internal} 
\title{Internal spatstat.model functions}
\alias{accumulateStatus}
\alias{active.interactions}
\alias{adaptcoef}
\alias{affine.msr}
\alias{areadelta2}
\alias{as.isf}
\alias{augment.msr}
\alias{bigvaluerule}
\alias{blankcoefnames}
\alias{bt.frame}
\alias{cannot.update}
\alias{check.separable}
\alias{coef.summary.kppm}
\alias{coef.summary.ppm}
\alias{coef.summary.slrm}
\alias{coef.vblogit} 
\alias{condSimCox}
\alias{damaged.ppm}
\alias{data.mppm}
\alias{deltasuffstat}
\alias{dfbetas.ppmInfluence}
\alias{diagnose.ppm.engine}
\alias{dim.msr}
\alias{dimnames.msr}
\alias{doMultiStraussHard}
\alias{dppDpcf}
\alias{dppmFixAlgorithm}
\alias{dppmFixIntensity}
\alias{DPPSaddle}
\alias{DPPSaddlePairwise}
\alias{enet.engine}
\alias{equalpairs}          
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{evalPairPotential}
\alias{evalPairwiseTerm}
\alias{expandDot}
\alias{extractAIC.slrm}
\alias{fakefii}
\alias{family.hackglmmPQL} 
\alias{family.vblogit} 
\alias{fill.coefs}
\alias{findCovariate}
\alias{fii}
\alias{flipxy.msr}
\alias{forbid.logi}
\alias{formula.hackglmmPQL} 
\alias{getdataname}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{getppmdatasubset}
\alias{getppmOriginalCovariates}
\alias{geyercounts}
\alias{geyerdelta2}
\alias{GLMpredict}
\alias{hackglmmPQL}
\alias{hasglmfit}
\alias{hierarchicalordering}
\alias{hiermat}
\alias{ho.engine}
\alias{illegal.iformula}
\alias{impliedpresence}
\alias{impliedcoefficients}
\alias{influence.ppmInfluence}
\alias{instantiate.interact}
\alias{interactionfamilyname}
\alias{intermaker}
\alias{is.expandable.ppm}
\alias{is.interact}
\alias{is.marked.mppm}    
\alias{is.marked.msr}    
\alias{is.marked.slrm}
\alias{is.mppm}
\alias{is.multitype.mppm}    
\alias{is.multitype.msr}    
\alias{is.multitype.slrm}    
\alias{is.poisson.mppm}
\alias{is.poisson.rppm}
\alias{Kpcf.kppm}
\alias{Kmodel.slrm}
\alias{killinteraction}
\alias{kppmComLik}
\alias{kppmMinCon}
\alias{kppmPalmLik}
\alias{kppmCLadap}
\alias{labels.ppm}
\alias{leverage.ppmInfluence}
\alias{logi.engine}
\alias{logLik.vblogit}
\alias{LurkEngine}
\alias{mapInterVars}
\alias{Mayer}
\alias{model.se.image}
\alias{modelFrameGam}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{mpl.usable}
\alias{newformula}
\alias{newstyle.coeff.handling}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{objsurfEngine}
\alias{optimConverged}
\alias{optimStatus}
\alias{optimNsteps}
\alias{outdated.interact}
\alias{oversize.quad}    
\alias{PairPotentialType}
\alias{partialModelMatrix}
\alias{pcfmodel.slrm}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.addvar}
\alias{plot.lurk}
\alias{plot.minconfit}
\alias{plot.parres}
\alias{plot.qqppm}
\alias{PoisSaddle}
\alias{PoisSaddleArea}
\alias{PoisSaddleGeyer}
\alias{PoisSaddlePairwise}
\alias{PPMmodelmatrix}
\alias{ppm.default}
\alias{ppmCovariates}
\alias{ppmDerivatives}
\alias{ppmInfluenceEngine}
\alias{predict.profilepl}
\alias{predict.vblogit}
\alias{printStatus}
\alias{printStatusList}
\alias{print.addvar}    
\alias{print.bt.frame}
\alias{print.diagppm}
\alias{print.detpointprocfamily}
\alias{print.detpointprocfamilyfun}
\alias{print.hierarchicalordering}
\alias{print.influence.ppm}       
\alias{print.interact}       
\alias{print.intermaker}       
\alias{print.isf}
\alias{print.leverage.ppm}
\alias{print.lurk}
\alias{print.minconfit}
\alias{print.mppm}
\alias{print.msr}
\alias{print.parres}
\alias{print.plotppm}
\alias{print.profilepl}
\alias{print.qqppm}
\alias{print.rppm}
\alias{print.summary.mppm}
\alias{print.summary.slrm}
\alias{print.vblogit}
\alias{quad.mppm}
\alias{quadBlockSizes}
\alias{reach.slrm}
\alias{reduceformula}
\alias{reincarnate.interact}
\alias{rescale.msr}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{rotate.msr}
\alias{SaddleApprox}
\alias{safeFiniteValue}
\alias{safePositiveValue}
\alias{scalardilate.msr}
\alias{shift.influence.ppm}              
\alias{shift.leverage.ppm}              
\alias{shift.msr}              
\alias{signalStatus}
\alias{simulate.profilepl}
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{slrmInfluence}
\alias{spatialCovariateEvidence.ppm}
\alias{spatialCovariateEvidence.slrm}
\alias{spatialCovariateUnderModel}
\alias{spatialCovariateUnderModel.ppm}
\alias{spatialCovariateUnderModel.kppm}
\alias{spatialCovariateUnderModel.dppm}
\alias{spatialCovariateUnderModel.slrm}
\alias{spatstatDPPModelInfo}
\alias{splitHybridInteraction}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{strausscounts}
\alias{stripGLMM}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summary.mppm}
\alias{summary.msr}
\alias{summary.profilepl}
\alias{summary.vblogit}
\alias{tweak.coefs}
\alias{unitname.msr}
\alias{unitname<-.msr}
\alias{update.ippm}
\alias{update.msr}
\alias{varcountEngine}
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{windows.mppm}




%%%%%%%
\description{
  Internal spatstat.model functions.
}
\usage{


accumulateStatus(x, stats)
active.interactions(object)
adaptcoef(new.coef, fitcoef, drop)
\method{affine}{msr}(X, mat, vec, \dots)
areadelta2(X, r, \dots, sparseOK)
as.isf(object)
augment.msr(x, \dots, sigma, recompute)
blankcoefnames(x)
bt.frame(Q, trend, interaction, \dots, covariates,
         correction, rbord, use.gam, allcovar)
bigvaluerule(objfun, objargs, startpar, \dots)
cannot.update(\dots)
check.separable(dmat, covname, isconstant, fatal)
\method{coef}{summary.kppm}(object, \dots)
\method{coef}{summary.ppm}(object, \dots)
\method{coef}{summary.slrm}(object, \dots)
\method{coef}{vblogit}(object, \dots) 
condSimCox(object, nsim, \dots, window, n.cond, w.cond,
           giveup, maxchunk, verbose, drop)
damaged.ppm(object)
data.mppm(x)
deltasuffstat(model, \dots,
              restrict, dataonly, sparseOK, quadsub,
              force, warn.forced, verbose, use.special)
\method{dfbetas}{ppmInfluence}(model, \dots)
diagnose.ppm.engine(object, \dots, type, typename, opt,
                         sigma, rbord, compute.sd, compute.cts,
                         envelope, nsim, nrank,
                         rv, oldstyle, splineargs, verbose)
\method{dim}{msr}(x)
\method{dimnames}{msr}(x)
doMultiStraussHard(iradii, hradii, types)
dppDpcf(model, \dots)
dppmFixIntensity(DPP, lambda, po)
dppmFixAlgorithm(algorithm, changealgorithm, clusters, startpar)
DPPSaddle(beta, fi)
DPPSaddlePairwise(beta, fi)
enet.engine(model, \dots, standardize, lambda, alpha, adaptive) 
equalpairs(U, X, marked=FALSE)
evalInteraction(X,P,E,interaction,correction,splitInf,\dots,
                precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,splitInf,\dots,
                Reach,precomputed,savecomputed)
evalPairPotential(X,P,E,pairpot,potpars,R)
evalPairwiseTerm(fint, d)
expandDot(f, dotvars)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, \dots)
fakefii(model)
\method{family}{hackglmmPQL}(object, \dots)
\method{family}{vblogit}(object, \dots)
fill.coefs(coefs, required)
findCovariate(covname, scope, scopename=NULL)
fii(interaction, coefs, Vnames, IsOffset, vnameprefix)
\method{flipxy}{msr}(X)
forbid.logi(object)
\method{formula}{hackglmmPQL}(x, \dots)
getdataname(defaultvalue, \dots, dataname)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
getppmdatasubset(object)
getppmOriginalCovariates(object)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
geyerdelta2(X,r,sat,\dots,sparseOK, correction)
GLMpredict(fit, data, coefs, changecoef, type)
hackglmmPQL(fixed, random, family, data, correlation, weights,
            control, niter, verbose, subset, \dots, reltol)
hasglmfit(object)
hierarchicalordering(i, s)
hiermat(x, h)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
illegal.iformula(ifmla, itags, dfvarnames)
impliedpresence(tags, formula, df, extranames=character(0))
impliedcoefficients(object, tag, new.coef)
\method{influence}{ppmInfluence}(model, \dots)
instantiate.interact(x, par)
interactionfamilyname(object)
intermaker(f, blank)
\method{is.expandable}{ppm}(x)  
is.interact(x)
\method{is.marked}{mppm}(X, \dots)
\method{is.marked}{msr}(X, \dots)
\method{is.marked}{rppm}(X, \dots)
\method{is.marked}{slrm}(X, \dots)
is.mppm(x)
\method{is.multitype}{mppm}(X, \dots)
\method{is.multitype}{msr}(X, \dots)
\method{is.multitype}{rppm}(X, \dots)
\method{is.multitype}{slrm}(X, \dots)
\method{is.poisson}{mppm}(x)
\method{is.poisson}{rppm}(x)
Kpcf.kppm(model, what)
\method{Kmodel}{slrm}(model, \dots)
killinteraction(model)
kppmComLik(X, Xname, po, clusters, control, stabilize, weightfun, rmax,
           algorithm, DPP, \dots, pspace)
kppmMinCon(X, Xname, po, clusters, control, stabilize, statistic, statargs,
           algorithm, DPP, \dots, pspace)
kppmPalmLik(X, Xname, po, clusters, control, stabilize, weightfun, rmax, 
           algorithm, DPP, \dots, pspace)
kppmCLadap(X, Xname, po, clusters, control, weightfun, 
           rmax, epsilon, DPP, algorithm, \dots,
           startpar, globStrat)
\method{labels}{ppm}(object, \dots)
\method{leverage}{ppmInfluence}(model, \dots)
logi.engine(Q, trend, interaction, \dots, 
            covariates, subsetexpr, clipwin,
            correction, rbord, covfunargs, allcovar, 
            vnamebase, vnameprefix, justQ, savecomputed, precomputed,
            VB)
\method{logLik}{vblogit}(object, \dots)
LurkEngine(object, type, cumulative, plot.sd, 
           quadpoints, wts, Z, subQset, 
           covvalues, resvalues, 
           clip, clipwindow, cov.is.im, covrange, 
           typename, covname,
           cl, clenv,
           oldstyle, check,
           verbose, nx, splineargs,
           envelope, nsim, nrank, Xsim,
           internal, checklength)
mapInterVars(object, subs, mom)
Mayer(fi, exponent)
model.se.image(fit, W, \dots, what, new.coef)
modelFrameGam(formula, \dots)
mpl.engine(Q, trend, interaction, \dots,
         covariates, subsetexpr, clipwin, covfunargs, correction,
	 rbord, use.gam, gcontrol,
         GLM, GLMfamily, GLMcontrol, famille,
         forcefit, nd, eps, allcovar, callstring,
         precomputed, savecomputed, preponly,
         rename.intercept, justQ, weightfactor)
mpl.get.covariates(covariates, locations, type, covfunargs, need.deriv)
mpl.prepare(Q, X, P, trend, interaction, covariates, 
            want.trend, want.inter, correction, rbord, Pname,
            callstring, \dots,
            subsetexpr,
            covfunargs, allcovar, precomputed, savecomputed,
            vnamebase, vnameprefix, warn.illegal, warn.unidentifiable,
            weightfactor, skip.border, clip.interaction, splitInf)
mpl.usable(x)
newformula(old, change, eold, enew, expandpoly, dotvars)
newstyle.coeff.handling(object)
nndcumfun(X, \dots, r)
no.trend.ppm(x)
objsurfEngine(objfun, optpar, objargs,
              \dots, dotargs, objname,
              ngrid, xlim, ylim, ratio, verbose)
optimConverged(x)
optimStatus(x, call)
optimNsteps(x)
outdated.interact(object)
oversize.quad(Q, \dots, nU, nX, p)
PairPotentialType(pairpot)
partialModelMatrix(X,D,model,callstring,\dots)
\method{pcfmodel}{slrm}(model, \dots)
ploterodewin(W1, W2, col.edge, col.inside, do.plot, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ, do.plot)
\method{plot}{addvar}(x, \dots, do.points=FALSE)
\method{plot}{lurk}(x, \dots, shade)
\method{plot}{minconfit}(x, \dots)
\method{plot}{parres}(x, \dots)
\method{plot}{qqppm}(x, \dots, limits=TRUE,
           monochrome=spatstat.options('monochrome'),
           limcol=if(monochrome) "black" else "red")
PoisSaddle(beta, fi)
PoisSaddleArea(beta, fi)
PoisSaddleGeyer(beta, fi)
PoisSaddlePairwise(beta, fi)
PPMmodelmatrix(object, data, \dots, subset, Q, keepNA, irregular,
               splitInf)
\method{ppm}{default}(Q, trend, interaction,
       \dots, covariates, data, covfunargs, subset, clipwin,
       correction, rbord, use.gam, method, forcefit,
       improve.type, improve.args, emend, project,
       prior.mean, prior.var,
       nd, eps, gcontrol, nsim, nrmh, start, control,
       verb, callstring)
ppmCovariates(model)
ppmDerivatives(fit, what, Dcovfun, loc, covfunargs)
ppmInfluenceEngine(fit, what, \dots, iScore, iHessian, iArgs,
              drop, method, fine, precomputed, sparseOK,
              fitname, multitypeOK, entrywise, matrix.action,
              dimyx, eps,
              geomsmooth)
\method{predict}{vblogit}(object, newdata, type, se.fit, dispersion,
                          terms, na.action, \dots)
\method{predict}{profilepl}(object, \dots)
printStatus(x, errors.only)
printStatusList(stats)
\method{print}{addvar}(x, \dots)
\method{print}{bt.frame}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{detpointprocfamily}(x, \dots)
\method{print}{detpointprocfamilyfun}(x, \dots)
\method{print}{hierarchicalordering}(x, \dots)
\method{print}{influence.ppm}(x, \dots)
\method{print}{interact}(x, \dots, family, brief, banner)       
\method{print}{intermaker}(x, \dots)
\method{print}{isf}(x, \dots)
\method{print}{leverage.ppm}(x, \dots)
\method{print}{lurk}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{mppm}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{parres}(x, \dots)
\method{print}{plotppm}(x, \dots)
\method{print}{profilepl}(x, \dots)
\method{print}{qqppm}(x, \dots)
\method{print}{rppm}(x, \dots)
\method{print}{summary.mppm}(x, \dots, brief)
\method{print}{summary.slrm}(x, \dots)
\method{print}{vblogit}(x, \dots)
quad.mppm(x)
quadBlockSizes(nX, nD, p, nMAX, announce)
\method{reach}{slrm}(x, \dots)
reduceformula(fmla, deletevar, verbose)
reincarnate.interact(object)
\method{rescale}{msr}(X, s, unitname)
resid4plot(RES, plot.neg, plot.smooth,
           spacing, outer, srange, monochrome, main,
           xlab, ylab, rlab, col.neg, col.smooth, \dots)
resid1plot(RES, opt, plot.neg, plot.smooth,
              srange, monochrome, main,
              add, show.all, do.plot, col.neg, col.smooth, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD,
            xlab,ylab, \dots, do.plot)
\method{rotate}{msr}(X, angle, \dots, centre)
SaddleApprox(beta, fi, approx)
safeFiniteValue(x, default)
safePositiveValue(x, default)
\method{scalardilate}{msr}(X, f, \dots)
\method{shift}{influence.ppm}(X, \dots)
\method{shift}{leverage.ppm}(X, \dots)
\method{shift}{msr}(X, \dots)
signalStatus(x, errors.only)
\method{simulate}{profilepl}(object, \dots)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby, clip)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
slrmInfluence(model, what, \dots)
\method{spatialCovariateEvidence}{ppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, jitterfactor,
          modelname, covname, dataname, subset, clip.predict) 
\method{spatialCovariateEvidence}{slrm}(model, covariate, \dots, lambdatype,
          jitter, jitterfactor,
          modelname, covname, dataname, subset) 
spatialCovariateUnderModel(model, covariate, \dots)
\method{spatialCovariateUnderModel}{ppm}(model, covariate, \dots)
\method{spatialCovariateUnderModel}{kppm}(model, covariate, \dots)
\method{spatialCovariateUnderModel}{dppm}(model, covariate, \dots)
\method{spatialCovariateUnderModel}{slrm}(model, covariate, \dots)
spatstatDPPModelInfo(model)
splitHybridInteraction(coeffs, inte)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
strausscounts(U,X,r,EqualPairs)
stripGLMM(object)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
\method{summary}{mppm}(object, \dots, brief=FALSE)
\method{summary}{msr}(object, \dots)
\method{summary}{profilepl}(object, \dots)
\method{summary}{vblogit}(object, \dots)
tweak.coefs(model, new.coef)
\method{unitname}{msr}(x) 
\method{unitname}{msr}(x) <- value 
\method{update}{ippm}(object, \dots, envir) 
\method{update}{msr}(object, \dots) 
varcountEngine(g, B, lambdaB, f, R, what)
versionstring.interact(object)
versionstring.ppm(object)
windows.mppm(x)



}
\details{
  These internal \pkg{spatstat.model} functions should not be called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.model} to the next.
}
\value{
The return values of these functions are not documented,
and may change without warning.
}
\keyword{internal}
