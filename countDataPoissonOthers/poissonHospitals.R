# Name: poissonHospitals.R
# Auth: umar niazi
# Date: 17/5/2019
# Desc: simple poisson regression models with complete and no pooling

library(LearnBayes)
setwd('countDataPoissonOthers/')
data("hearttransplants")
dfData = hearttransplants
dfData$hospital = factor(1:nrow(dfData))

fit.1 = glm(y ~ 1 , data=dfData, family='poisson', offset=log(e))
summary(fit.1)

## utility functions
# calculates the gamma prior parameters for a poisson sampling distribution
# see page 5 in notes here: https://www.evernote.com/shard/s288/res/659dc700-ccb9-457e-aea6-aa54bc3abbb9
# and for an example see page 154, chapter on Hierarchical Modeling Bayesian Computation with R.
## DESC: using the poisson sampling model, the data vector is used to count values of alpha (shape), beta (rate)
## parameters for the gamma prior
getalphabeta.poisson = function(lambda){
  m = mean(lambda)
  v = var(lambda)
  alpha = (m^2)/v
  beta = alpha/m
  return(c(alpha=alpha, beta=beta))
}

simGammaPost = function(data, exposure, prior){
  alpha = prior['alpha']+sum(data)
  beta = prior['beta']+sum(exposure)
  return(rgamma(1000, shape = alpha, rate = beta))
}

simPostPredictPoisson = function(post, exposure, len, nc=20){
  mDraws = matrix(NA, nrow = len, ncol=nc)
  for (i in 1:nc){
    p = sample(post, size = 1)
    mDraws[,i] = rpois(len, p*exposure)
  }
  return(mDraws)
}

ebayes.prior = getalphabeta.poisson(dfData$y)
fit.2 = simGammaPost(dfData$y, dfData$e, ebayes.prior)

mSim.fit.1 = simulate(fit.1, 20)
mSim.fit.2 = simPostPredictPoisson(fit.2, dfData$e, nrow(dfData), 20)
hist(dfData$y, xlab='', ylab='', main='', prob=T)
temp = apply(mSim.fit.1, 2, function(x) lines(density(x), col='green'))
temp = apply(mSim.fit.2, 2, function(x) lines(density(x), col='blue'))

## complete pooling analysis
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='poissonRegressionNoPooling.stan')

m = model.matrix(y ~ 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'mu'),
                    cores=2)
print(fit.stan, c('betas'), digits=3)

mFitted = extract(fit.stan)$mu
mFitted = mFitted[sample(1:nrow(mFitted), 20, replace = F),]
mSim.fit.stan = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})

temp = apply(mSim.fit.stan, 2, function(x) lines(density(x), col='red'))

###########################################################################################################
########## Second analysis with adding grouping information, includes no pooling and partial pooling
###########################################################################################################
fit.1.2 = glm(y ~ 1 + hospital , data=dfData, family='poisson', offset=log(e))
summary(fit.1.2)
anova(fit.1, fit.1.2)

fit.2.2 = sapply(1:nrow(dfData), function(x){
  return(simGammaPost(dfData$y[x], dfData$e[x], ebayes.prior))
})


simPostPredictPoisson2 = function(post, exposure, len, nc=20){
  mDraws = matrix(NA, nrow = len, ncol=nc)
  for (i in 1:nc){
    p = sample(1:nrow(post), size = 1)
    mDraws[,i] = rpois(len, post[p,] * exposure)
  }
  return(mDraws)
}

mSim.fit.1.2 = simulate(fit.1.2, 20)
mSim.fit.2.2 = simPostPredictPoisson2(fit.2.2, dfData$e, nrow(dfData), 20)
hist(dfData$y, xlab='', ylab='', main='', prob=T)
temp = apply(mSim.fit.1.2, 2, function(x) lines(density(x), col='green'))
temp = apply(mSim.fit.2.2, 2, function(x) lines(density(x), col='blue'))

## no pooling analysis
m = model.matrix(y ~ 1 + hospital, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'mu'),
                    cores=2)
print(fit.stan, c('betas'), digits=3)

mFitted = extract(fit.stan)$mu
dim(mFitted)
mFitted = mFitted[sample(1:nrow(mFitted), 20, replace = F),]
mSim.fit.stan = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})

temp = apply(mSim.fit.stan, 2, function(x) lines(density(x), col='red'))

### partial pooling, multilevel model
stanDso = rstan::stan_model(file='poissonRegressionPartialPooling_1.stan')

m = model.matrix(y ~ hospital - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('populationMean', 'sigmaRan', 'betas', 'mu'),
                    cores=2)
print(fit.stan, c('betas', 'populationMean', 'sigmaRan'), digits=3)

pairs(fit.stan, pars = c("sigmaRan", "populationMean", "lp__"))
# some diagnostics for stan
traceplot(fit.stan, c('sigmaRan'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('populationMean'), ncol=1, inc_warmup=F)

library(lme4)
fit.lme = glmer(y ~ 1 + (1 | hospital), offset=log(dfData$e), data=dfData, family=poisson(link = "log"))
## compare lme4 and stan
summary(fit.lme)

## plot posterior predictive values
mFitted = extract(fit.stan)$mu
dim(mFitted)
mFitted = mFitted[sample(1:nrow(mFitted), 20, replace = F),]
mSim.fit.stan = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})

temp = apply(mSim.fit.stan, 2, function(x) lines(density(x), col='yellow2'))

## can we replicate the figure from the book
## complete pooling
mSim.fit.1 = simulate(fit.1, 1000)
mSim.fit.2 = simPostPredictPoisson(fit.2, dfData$e, nrow(dfData), 1000)

par(mfrow=c(1,3))
hist(as.numeric(mSim.fit.1[94,]))
abline(v = dfData$y[94])
hist(as.numeric(mSim.fit.2[94,]))
abline(v = dfData$y[94])
hist(mSim.fit.stan[94,])
abline(v = dfData$y[94])

## no pooling
mSim.fit.1.2 = simulate(fit.1.2, 1000)
mSim.fit.2.2 = simPostPredictPoisson2(fit.2.2, dfData$e, nrow(dfData), 1000) 

hist(as.numeric(mSim.fit.1.2[94,]))
abline(v = dfData$y[94])
hist(as.numeric(mSim.fit.2.2[94,]))
abline(v = dfData$y[94])
hist(mSim.fit.stan[94,])
abline(v = dfData$y[94])

## partial pooling


