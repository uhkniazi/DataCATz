# Name: linearRegressionRandomSlope.R
# Auth: uhkniazi
# Date: 15/12/2017
# Desc: linear regression, using fixed and random effects (2 level model)


## following some guideline from https://rpubs.com/kaz_yos/stan-multi-1
## Simulate the data according to guide from the link above
# some adjustments made to data simulation, as Rossella has very good eyes and found a typo
set.seed(201503031)

## Random intercepts
## 5 groups with 5 means
## groups are any form of clusters or grouping in the data, e.g. repeated measurements 
## from the same individual
# we assume that each group has its own intercept which is distributed
# around the population intercept as N(mean=0, sd=1)
gSigmaRan = 1
gInterceptRan = rnorm(n = 5, mean = 0, sd = gSigmaRan)

## At the Population level or the second level of the model
## the intercept is 3
gPopulationIntercept = 3
## Varying intercepts
## now adding the population intercept and the random variation for each group
## we should be able to generate varying intercepts for each group
## we have 6 observations per group
gGroupIntercept  = rep(gPopulationIntercept + gInterceptRan, 6)

## Cluster indicator
## this is a mapping variable - used for tracking the population intercepts
iClusterIndicator = rep(1:5, times = 6)
## population level error, how much
## each data point varies
gSigmaPop = 2
# we have 5 groups X 6 per group = 30 observations
iErrors = rnorm(n = 30, mean = 0, sd = gSigmaPop)

## We generate some random predictor and
## generate some data
## predictor variable
gPredictor = runif(n = 30, min = 0, max = 20)
## Population slope.
## we are not using individual slops for groups, which is also an option
## and we can have correlations between intercepts, slopes or various groups
gSlope = 1.5
## Outcome
gResponse = gGroupIntercept + gSlope*gPredictor + iErrors

dfData = data.frame(resp = gResponse, pred=gPredictor, group=factor(iClusterIndicator), groupIndex=iClusterIndicator)


library(lattice)
xyplot(resp ~ pred | group, data=dfData, type=c('g', 'p', 'r'), 
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7))

xyplot(resp ~ pred, data=dfData, type=c('g', 'p', 'r'), groups= group,
       index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',# layout=c(8,2), 
       par.strip.text=list(cex=0.7))

xyplot(resp ~ pred, data=dfData, type=c('g', 'p', 'r'))

################### first use normal regression and lme4 package
fit.lm = lm(resp ~ pred, data=dfData)
summary(fit.lm)

# residual se = 2.616
# a model with uncorrelated slope
library(lme4)
fit.lme = lmer(resp ~ 1 + pred + (1 | group) #+ (0 + pred | group)
               , data=dfData, REML=F)
summary(fit.lme)

## random effect or group level sd = 1.632
## population level sd = 1.9434
## total = 1.6322+1.9434 = 3.5756

#################### use stan to generate MCMC sample
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='BayesRegression/linearRegressionRandomEffectsRandomSlope.stan')

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters=nlevels(dfData$group), 
                 NgroupMap=dfData$groupIndex, X=dfData$pred,
                 y=dfData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, #pars=c('betas', 'sigmaRan', 'sigmaPop', 'rGroupsJitter'),
                    cores=2)
print(fit.stan, digits=3)


####### now lets write our own model function
library(numDeriv)
library(LearnBayes)
#library(optimx) optional optimiser to use
# op = optimx(start, mylogpost, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), method='Nelder-Mead', data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}


mylogpost = function(theta, data){
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  groupIndex = data$groupIndex  ## mapping variable to map each random effect to its respective response variable
  ## parameters to track/estimate
  sigmaRan = exp(theta['sigmaRan']) # random effect scale/sd
  sigmaPop = exp(theta['sigmaPop']) # population level sd
  betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
  iGroupsJitter = theta[grep('ran', names(theta))]# random effects jitters for the group deflections
  
  ## random effect jitter for the population intercept
  # each group contributes a jitter centered on 0
  # population slope + random jitter
  ivBetaRand = betas[1] + iGroupsJitter
  # create a matrix of betas with the new interceptr/unique intercept for each random effect
  ivIntercept = ivBetaRand[groupIndex] # expand this intercept using the mapping variable
  iFitted = as.matrix(mModMatrix[,2:ncol(mModMatrix)]) %*% betas[2:ncol(mModMatrix)]
  iFitted = ivIntercept + iFitted # using identity link so no transformation required
  ## priors
  # sigma priors
  lhp = dunif(sigmaRan, 0, 2, log=T) + dunif(sigmaPop, 0, 5, log=T)
  # random effects prior
  lran = sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T))
  # priors for the betas
  lp = sum(dcauchy(betas, 0, 10, log=T))
  # write the likelihood function
  lik = sum(dnorm(resp, iFitted, sigmaPop, log=T))
  val = lhp + lran + lp + lik
  return(val)
}

## select starting values for betas using linear regression
betas = coef(lm(resp ~ pred, data=dfData))
start = c('sigmaPop'=log(2), 'sigmaRan'=log(2), 'betas'=as.numeric(betas), 'ran'=rep(0, times=nlevels(dfData$group)))
lData = list('resp'=dfData$resp, 'mModMatrix'=model.matrix(resp ~ pred, data=dfData), 'groupIndex'=dfData$groupIndex)
mylogpost(start, lData)

fit.lap = mylaplace(mylogpost, start, lData)

### lets use the results from laplace and SIR sampling with a t proposal density
### to take a sample
tpar = list(m=fit.lap$mode, var=fit.lap$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 5000, lData)
apply(s, 2, mean)[-c(1:2)]
exp(apply(s, 2, mean)[1:2])
mSir = s


############### you can compare the results from stan and lmer and the logposterior function
# Inference for Stan model: linearRegressionRandomEffects.
# 4 chains, each with iter=5000; warmup=2500; thin=1; 
# post-warmup draws per chain=2500, total post-warmup draws=10000.
# 
#                    mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff  Rhat
# betas[1]           1.775   0.038 1.620  -1.736   0.935   1.841   2.738   4.786  1853 1.000
# betas[2]           1.549   0.001 0.070   1.411   1.503   1.549   1.597   1.685  5649 1.001
# sigmaRan           2.961   0.051 2.471   0.875   1.728   2.414   3.487   8.093  2377 1.001
# sigmaPop           2.094   0.005 0.327   1.563   1.863   2.056   2.280   2.843  4928 1.001
# rGroupsJitter[1]   0.503   0.036 1.623  -2.559  -0.419   0.434   1.326   3.996  2083 1.000
# rGroupsJitter[2]   0.901   0.037 1.635  -2.068  -0.055   0.808   1.733   4.493  1997 1.000
# rGroupsJitter[3]  -0.636   0.036 1.618  -3.838  -1.516  -0.646   0.225   2.813  2036 1.000
# rGroupsJitter[4]   2.140   0.036 1.660  -0.808   1.131   2.013   3.008   5.892  2077 1.000
# rGroupsJitter[5]  -2.448   0.035 1.652  -5.842  -3.378  -2.406  -1.502   0.792  2179 1.000

# > summary(fit.lme)
# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: resp ~ pred + (1 | group)
# Data: dfData
# 
# AIC      BIC   logLik deviance df.resid 
# 141.3    146.9    -66.6    133.3       26 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.5264 -0.7096 -0.2699  0.4328  2.6193 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# group    (Intercept) 2.664    1.632   
# Residual             3.777    1.943   
# Number of obs: 30, groups:  group, 5
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  1.89979    1.00879   1.883
# pred         1.54593    0.06604  23.408
# 
# Correlation of Fixed Effects:
#   (Intr)
# pred -0.594
# > t(ranef(fit.lme)$group)
#              1         2          3       4         5
# (Intercept) 0.3888984 0.7627657 -0.6863616 1.94236 -2.407663

# > fit.lap$mode[-c(1:2)]
# betas1     betas2       ran1       ran2       ran3       ran4       ran5 
# 1.9037216  1.5456933  0.3897490  0.7787912 -0.7075750  1.8900023 -2.3491993 
# > exp(fit.lap$mode)[1:2]
# sigmaPop sigmaRan 
# 1.805638 1.463412

## SIR sampling - these will slighly change each time you run it
# > apply(s, 2, mean)[-c(1:2)]
# betas1     betas2       ran1       ran2       ran3       ran4       ran5 
# 1.8598225  1.5462439  0.3754429  0.7316874 -0.5849880  1.7242024 -2.1096244 
# > exp(apply(s, 2, mean)[c(1:2)])
# sigmaPop sigmaRan 
# 2.082105 1.351142 

######################### calculate model fit parameters
## Log likelihood and AIC
## define the log predictive density function
## this function is a little awkward to define, as
## we need to perform notational acrobatics to create the covariance matrix
## i figured this out after reading the post here
## https://stats.stackexchange.com/questions/271903/understand-marginal-likelihood-of-mixed-effects-models
### calculate model fits
## first write the log predictive density function
lpd = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  group = data$group
  sigmaRan = exp(data$sigmaRan) # random effect scale/sd
  sigmaPop = exp(data$sigmaPop) # population level sd
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  ######## perform notational acrobatics
  # construct the variance covariance matrix
  z = model.matrix(resp ~ 0 + group)
  zt = t(z)
  mRan = diag(sigmaRan^2, nrow=nlevels(group), ncol=nlevels(group))
  mRan = z %*% mRan %*% zt
  mPop = diag(sigmaPop^2, nrow=length(resp), ncol=length(resp))
  mCov = mRan + mPop
  ## likelihood function with posterior theta
  return(dmnorm(resp, iFitted, mCov, log=T))
}

lData$group = dfData$group
lData$sigmaRan = log(1.632)
lData$sigmaPop = log(1.943)
theta = fixef(fit.lme)

lpd(theta, lData)
## compare with lme
logLik(fit.lme)

## use the arm library to extract these directly from lme object
library(arm)
extractAIC(fit.lme)
extractDIC(fit.lme)

######## as the lpd function gives same results, now calculate using our calculations
lData$sigmaRan = apply(mSir, 2, mean)[2]
lData$sigmaPop = apply(mSir, 2, mean)[1]
theta = apply(mSir, 2, mean)[3:4]

lpd(theta, lData)

## calculate AIC
iAIC = (lpd(theta, lData) - 4) * -2

## calculate DIC
eLPD = mean(sapply(1:5000, function(x) {
  lData$sigmaRan = mSir[x, 2]
  lData$sigmaPop = mSir[x, 1]
  theta = mSir[x, 3:4]
  lpd(theta, lData)}))
# calculate lpd(E(theta)) and pDIC
lData$sigmaRan = apply(mSir, 2, mean)[2]
lData$sigmaPop = apply(mSir, 2, mean)[1]
theta = apply(mSir, 2, mean)[3:4]
pDIC = 2 *(lpd(theta, lData) - eLPD)
iDIC = (lpd(theta, lData) - pDIC) * -2



