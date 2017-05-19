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
lData = list(pred=dfData$IncomeGrowth, resp=dfData$VoteShare)

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
sigSample = rnorm(10000, fit$mode['sigma'], se['sigma'])
beta0Sample = rnorm(1000, fit$mode['beta0'], se['beta0'])
beta1Sample = rnorm(1000, fit$mode['beta1'], se['beta1'])

## plot the distributions
par(mfrow=c(2,2))
hist(exp(sigSample), main='Sigma distribution', xlab='')
hist(beta0Sample, main='Beta0 distribution', xlab='')
hist(beta1Sample, main='Beta1 distribution', xlab='')

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
lpdSample = sapply(1:1000, function(x) lpd(beta0Sample[x], beta1Sample[x], sigSample[x]))
summary(lpdSample)


temp = sapply(1:1000, function(x) {
  s = c('sigma'=sigSample[x], beta0Sample[x], beta1Sample[x])
  mylogpost(s, lData)
})

