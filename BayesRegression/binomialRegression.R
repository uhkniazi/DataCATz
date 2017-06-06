# Name: binomialRegression.R
# Auth: uhkniazi
# Date: 06/06/2017
# Desc: Binomial regression, using fixed and random effects (2 level model)


## using the contraception data set
## see page 120 Chapter 6 of Lme4 book.
library(lattice)
library(lme4)
library(mlmRev)
library(car)
library(LearnBayes)

data("Contraception")

str(Contraception)

## we plot the data as the author does in this book
xyplot(ifelse(use == 'Y', 1, 0) ~ age|urban, data = Contraception, groups = livch, lty=1:4, col=1,
       type=c('g', 'smooth'), key = list(text=list(levels(Contraception$livch)), lines=list(lty=1:4),
                                         columns=4, cex=0.8),
       xlab='Age', ylab='Contraception Use', main=list(label='Contraception use Given Urban', cex=0.8))

## Figure suggests following things:
## 1- Contraception use and age have a quadradic trend.
## 2- middle age range women are more likely to use contraception
## 3- Urban women are more likely to use contraception
## 4- women with 0 childeren are less likely to use contraception.
## 5- women who have children are less different from each other

## using point 5, merge the groups 1, 2, and 3 in the variable livch
Contraception$children = factor(Contraception$livch != 0, labels=c('N', 'Y'))

str(Contraception)

xyplot(ifelse(use == 'Y', 1, 0) ~ age|urban, data = Contraception, groups = children, lty=1:2, col=1,
       type=c('g', 'smooth'), key = list(text=list(levels(Contraception$children)), lines=list(lty=1:2),
                                         columns=2, cex=0.8),
       xlab='Age', ylab='Contraception Use', main=list(label='Contraception use Given Urban', cex=0.8))

## we should perhaps also use an interaction term between age and number of children, as the slopes and interceptrs
## in cases with children are different to those without children.

## lets fit a model using the binomial glm without random effects
fit.1 = glm(use ~ age + I(age^2) + urban + children, data=Contraception, family = binomial(link='logit'))
summary(fit.1)


logit.inv = function(p) {exp(p)/(exp(p)+1) }

## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood 
  lp = dcauchy(betas[1], 0, 10, log=T) + sum(dcauchy(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

lData = list(resp=ifelse(Contraception$use == 'N', 0, 1), mModMatrix=model.matrix(use ~ age + I(age^2) + urban + children, data=Contraception))
start = c(rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.2 = laplace(mylogpost, start, lData)
se = sqrt(diag(fit.2$var))

### lets take a sample from this 
## parameters for the multivariate t density
tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 1000, lData)
colnames(s) = colnames(lData$mModMatrix)
apply(s, 2, mean)
apply(s, 2, sd)
pairs(s, pch=20)
fit.2$sir = s

#######################################################################
########### lets fit a more complex random effects model
fit.3 = glmer(use ~ age + I(age^2) + urban + children + (1 | district), data=Contraception, family = binomial)
summary(fit.3)
## ignore the convergence error

############################# test with stan and compare results
library(rstan)
stanDso = rstan::stan_model(file='BayesRegression/binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Nclusters=length(unique(lData$groupIndex)), 
                 NgroupMap=lData$groupIndex, Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas', 'sigmaRan'))
print(fit.stan)
plot(fit.stan)
traceplot(fit.stan, ncol=1, nrow=6, inc_warmup=F)
stan_diag(fit.stan)
## some sample diagnostic plots
library(coda)
oCoda = As.mcmc.list(fit.stan)
xyplot(oCoda[[1]])
autocorr.plot(oCoda[[1]])

# write a new log posterior function with random effects
myloglike.random = function(theta, data){
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  groupIndex = data$groupIndex  ## mapping variable to map each random effect to its respective response variable
  ## parameters to track/estimate
  sigmaRan = exp(theta['sigmaRan']) # random effect scale/sd
  betas = theta[2:(ncol(mModMatrix)+1)] # vector of betas i.e. regression coefficients for population
  iGroupsJitter = theta[-c(1,(ncol(mModMatrix)+1))] # random effects jitters for the group deflections
  
  ## random effect jitter for the population intercept
  # each group contributes a jitter centered on 0
  # population slope + random jitter
  ivBetaRand = betas[1] + iGroupsJitter
  # create a matrix of betas with the new interceptr/unique intercept for each random effect
  ivIntercept = ivBetaRand[groupIndex] # expand this intercept
  iFitted = mModMatrix[,2:ncol(mModMatrix)] %*% betas[2:ncol(mModMatrix)]
  iFitted = ivIntercept + iFitted
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## first level, log priors
  if (is.nan(sigmaRan) | sigmaRan <= 0) return(-Inf)
  lran = sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T))
  lp = dcauchy(betas[1], 0, 10, log=T) + sum(dcauchy(betas[-1], 0, 10, log=T))
  # write the likelihood function
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lp + lran + lik
  return(val)
}

lData$groupIndex = as.numeric(Contraception$district)

start = c(sigmaRan=log(sd(lData$resp)), fit.2$mode, rep(0, times=nlevels(Contraception$district)))

myloglike.random(start, lData)

fit.4 = laplace(myloglike.random, start, lData)
op = optim(start, myloglike.random, gr=NULL,
           control = list(fnscale = -1, maxit=1000), method='SANN', data=lData)
start = op$par
### define a custom laplace function 
library(numDeriv)

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), method='CG', data=data)
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

fit.4 = mylaplace(myloglike.random, start, lData)



start['mu'] = mean(extract(fit.stan)$mu)
start['tau'] = log(mean(extract(fit.stan)$tau))
start

mylaplace(mylogpost, start, lData)



