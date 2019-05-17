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

## no-pooling analysis
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='poissonRegressionNoPooling.stan')

m = model.matrix(y ~ 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'mu'),
                    cores=2)
print(fit.stan, c('betas', 'mu'), digits=3)






