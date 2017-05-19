# Name: breadAndPeace.R
# Auth: uhkniazi
# Date: 18/05/2017
# Desc: calculating various predictive quantities to assist in model comparisons

library(LearnBayes)

p.old = par()
dfData = read.csv('BayesModelChecks/breadAndPeace.csv', header=T)
str(dfData)

dim(dfData)

fit.lm = lm(VoteShare ~ IncomeGrowth, data=dfData)
summary(fit.lm)

## calculate using completing the square
## see Gelman [2013] page 356
lData = list(pred=dfData$IncomeGrowth, resp=dfData$VoteShare)
X = model.matrix(lData$resp ~ lData$pred)
bs = solve(t(X) %*% X) %*% t(X) %*% lData$resp
mVbetas = solve(t(X) %*% X)
sig2 = 1/(15-2) * t(lData$resp - X %*% bs) %*% (lData$resp - X %*% bs)
sqrt(mVbetas * as.numeric(sig2)) # compare with SE calculated using optimiser

### first get the estimates using stan and mcmc
library(rstan)
stanDso = rstan::stan_model(file='BayesModelChecks/breadAndPeace.stan')

lStanData = list(Ntotal=nrow(X), Ncol=ncol(X), X=X, y=lData$resp)
fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=4, pars=c('betas', 'sigma'))
print(fit.stan)

mStan = do.call(cbind, extract(fit.stan))[,-4]
colnames(mStan) = c('beta0', 'beta1', 'sigma')
apply(mStan, 2, mean)
apply(mStan, 2, sd)

## write the model with the log posterior function
# use the likelihood with non-informative priors
mylogpost = function(theta, data){
  sigma = exp(theta['sigma'])
  betas = theta[-1]
  x = data$pred
  y = data$resp
  
  # get the predicted or fitted value 
  mModMatrix = model.matrix(y ~ x)
  mCoef = matrix(betas, nrow = length(betas))
  iFitted = mModMatrix %*% mCoef # using identity link
  ## likelihood function
  llik = sum(dnorm(y, iFitted, sigma, log=T))
  ## prior, non informative
  lprior = 1
  lpost = llik + lprior
  return(lpost)
}

## choose some starting values 
start = c('sigma' = log(sd(dfData$VoteShare)), 'beta0'=0, 'beta1'=0)
mylogpost(start, lData)

## lets optimize using an optimizer
fit = laplace(mylogpost, start, lData)
fit

se = sqrt(diag(fit$var))
se
fit$mode+1.96*se
fit$mode-1.96*se

## coef and sd
fit$mode[-1]
exp(fit$mode[1])

## lets take a sample from this
#s = rmnorm(1000, fit$mode, fit$var)
## parameters for the multivariate t density
tpar = list(m=fit$mode, var=fit$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 10000, lData)
#s = rmt(10000, fit$mode, S = fit$var)
sigSample = s[,'sigma']
beta0Sample = s[,'beta0']
beta1Sample = s[,'beta1']

## plot the distributions
par(mfrow=c(2,2))
hist(exp(sigSample), main='Sigma distribution', xlab='', prob=T)
lines(density(exp(mStan[,'sigma'])))
hist(beta0Sample, main='Beta0 distribution', xlab='', prob=T)
lines(density(mStan[,'beta0']))
hist(beta1Sample, main='Beta1 distribution', xlab='', prob=T)
lines(density(mStan[,'beta1']))

# We are interested in prediction accuracy for two reasons: first, to measure the performance
# of a model that we are using; second, to compare models. Our goal in model
# comparison is not necessarily to pick the model with lowest estimated prediction error or
# even to average over candidate models ... but at least to put different models
# on a common scale. Even models with completely different parameterizations can be used
# to predict the same measurements. [Gelman 2013]


############# We will calculate some model checking parameters
## first write the log predictive density function
lpd = function(beta0, beta1, sig){
  sigma = exp(sig)
  x = lData$pred
  y = lData$resp
  # get the predicted or fitted value 
  mModMatrix = model.matrix(y ~ x)
  mCoef = matrix(c(beta0, beta1), nrow = 2)
  iFitted = mModMatrix %*% mCoef # using identity link
  ## likelihood function with posterior theta
  return(sum(dnorm(y, iFitted, sigma, log=T)))
}

## get a distribution of observed log predictive density
lpdSample = sapply(1:10000, function(x) lpd(beta0Sample[x], beta1Sample[x], sigSample[x]))
summary(lpdSample)
max(lpdSample) - mean(lpdSample)
# The mean of the posterior distribution of the log predictive density is −42.0, and the difference
# between the mean and the maximum is 1.7, which is close to the value of 2/3 that would
# be predicted from asymptotic theory, given that 3 parameters are being estimated. [Gelman 2013]


temp = sapply(1:1000, function(x) {
  s = c('sigma'=sigSample[x], beta0Sample[x], beta1Sample[x])
  mylogpost(s, lData)
})

