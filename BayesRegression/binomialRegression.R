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

p.old = par()
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

########################### analytical calculation and meaning of regression coefficient
fit.0 = glm(use ~ children, data=Contraception, family = binomial(link='logit'))
summary(fit.0)

## contingency table
mCont = as.matrix(table(Contraception$use, Contraception$children))
dimnames(mCont) = list(c('UseN', 'UseY'), c('ChildN', 'ChildY'))

iMarginal.use = rowSums(mCont)
iMarginal.child = colSums(mCont)



## we should perhaps also use an interaction term between age and number of children, as the slopes and intercepts
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
  lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

lData = list(resp=ifelse(Contraception$use == 'N', 0, 1), mModMatrix=model.matrix(use ~ age + I(age^2) + urban + children, data=Contraception))
start = c(rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.2 = laplace(mylogpost, start, lData)
fit.2
data.frame(coef(fit.1), fit.2$mode)
se = sqrt(diag(fit.2$var))


### lets take a sample from this 
## parameters for the multivariate t density
tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 5000, lData)
colnames(s) = colnames(lData$mModMatrix)
apply(s, 2, mean)
apply(s, 2, sd)
pairs(s, pch=20)
fit.2$sir = s

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='BayesRegression/binomialRegression.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas'))
# some diagnostics for stan
print(fit.stan, digits=3)
plot(fit.stan)
traceplot(fit.stan, ncol=1, nrow=6, inc_warmup=F)
#stan_diag(fit.stan)
## some sample diagnostic plots
library(coda)
oCoda = As.mcmc.list(fit.stan)
xyplot(oCoda[[1]])
autocorr.plot(oCoda[[1]])

## how do the samples from sir and stan compare
s1 = fit.2$sir
s2 = extract(fit.stan)$betas
par(mfrow=c(2,3))

plot(density(s2[,1]), col=1, main='intercept')
lines(density(s1[,1]), col=2)

plot(density(s2[,2]), col=1, main='age')
lines(density(s1[,2]), col=2)

plot(density(s2[,3]), col=1, main='age^2')
lines(density(s1[,3]), col=2)

plot(density(s2[,4]), col=1, main='urbanY')
lines(density(s1[,4]), col=2)

plot(density(s2[,5]), col=1, main='childrenY')
lines(density(s1[,5]), col=2)


### calculate model fits
## first write the log predictive density function
lpd = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## likelihood function with posterior theta
  return(sum(dbinom(resp, 1, iFitted, log=T)))
}

## averages of posterior from stan sample
post = apply(s2, 2, mean)

# AIC
iAIC = (lpd(post, lData) - 5) * -2
AIC(fit.1)

i = sample(1:10000, size = 1000, replace = F)
## get a distribution of observed log predictive density
s = s2[i,]
lpdSample = sapply(1:1000, function(x) lpd(s[x,], lData))
summary(lpdSample)
max(lpdSample) - mean(lpdSample)
## 5/2 = 2.5
# The mean of the posterior distribution of the log predictive density and the difference
# between the mean and the maximum is close to the value of 5/2 that would
# be predicted from asymptotic theory, given that 5 parameters are being estimated. [Gelman 2013]

## calculate WAIC
## DIC 
## pDIC are the effective number of parameters
## 2 * [lpd(Expectation(theta)) - Expectation(lpd(Sample of thetas from posterior))]
# calculate E(lpd(theta))
eLPD = mean(sapply(1:10000, function(x) lpd(s2[x,], lData)))
# calculate lpd(E(theta)) and pDIC
pDIC = 2 *(lpd(post, lData) - eLPD)
iDIC = (lpd(post, lData) - pDIC) * -2
# calculate ELPD

## calculate one data point at a time
## log pointwise predictive density
lppd = function(theta, data){
  betas = t(theta) # matrix of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = matrix(data$mModMatrix, nrow = 1, byrow = T)
  # calculate fitted value
  iFitted = as.vector(mModMatrix %*% betas)
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## likelihood function with posterior theta
  return(mean(dbinom(resp, 1, iFitted, log=F)))
}

ilppd = sum(log(sapply(seq_along(lData$resp), function(x) {
  d = list(resp=lData$resp[x], mModMatrix = lData$mModMatrix[x,])
  lppd(s2, d)
})))

## effective numbers of parameters pWAIC1
pWAIC1 = 2 * (ilppd - eLPD)

iWAIC = -2 * (ilppd - pWAIC1)




#######################################################################
########### lets fit a more complex random effects model
fit.3 = glmer(use ~ age + I(age^2) + urban + children + (1 | district), data=Contraception, family = binomial)
summary(fit.3)
## ignore the convergence error

lData$groupIndex = as.numeric(Contraception$district)
############################# test with stan and compare results
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='BayesRegression/binomialRegressionRandomEffects.stan')

lStanData = list(Ntotal=length(lData$resp), Nclusters=length(unique(lData$groupIndex)), 
                 NgroupMap=lData$groupIndex, Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, 
                    pars=c('betas', 'sigmaRan', 'rGroupsJitter'), cores=3)
print(fit.stan, digits=3)

## extract the data
s2 = extract(fit.stan) 
colnames(s2$betas) = paste('betas', 1:ncol(s2$betas), sep='')
colnames(s2$rGroupsJitter) =  paste('ran', 1:ncol(s2$rGroupsJitter), sep='')
s2 = cbind(s2$sigmaRan, s2$betas, s2$rGroupsJitter)
colnames(s2)[1] = 'sigmaRan'
dim(s2)

### calculate model fits
## first write the log predictive density function
lpd = function(theta, data){
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  groupIndex = data$groupIndex  ## mapping variable to map each random effect to its respective response variable
  ## parameters to track/estimate
  sigmaRan = exp(theta['sigmaRan']) # random effect scale/sd
  betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
  iGroupsJitter = theta[grep('ran', names(theta))]# random effects jitters for the group deflections
  
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
  # write the likelihood function
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  return(lik)
}

## averages of posterior from stan sample
post = apply(s2, 2, mean)

# AIC
iAIC = (lpd(post, lData) - length(post)) * -2
AIC(fit.3)

## get a distribution of observed log predictive density
lpdSample = sapply(1:nrow(s2), function(x) lpd(s2[x,], lData))
summary(lpdSample)
max(lpdSample) - mean(lpdSample)
## 38.5/2 = 19.25 (38 are effective number of parameters calculated from pDIC)
# The mean of the posterior distribution of the log predictive density and the difference
# between the mean and the maximum is close to the value of 5/2 that would
# be predicted from asymptotic theory, given that 5 parameters are being estimated. [Gelman 2013]

## calculate WAIC and DIC
## DIC 
## pDIC are the effective number of parameters
## 2 * [lpd(Expectation(theta)) - Expectation(lpd(Sample of thetas from posterior))]
# calculate E(lpd(theta))
eLPD = mean(sapply(1:nrow(s2), function(x) lpd(s2[x,], lData)))
# calculate lpd(E(theta)) and pDIC
pDIC = 2 *(lpd(post, lData) - eLPD)
iDIC = (lpd(post, lData) - pDIC) * -2
# calculate ELPD

## calculate one data point at a time
## log pointwise predictive density
lppd = function(theta, data){
  betas = t(theta) # matrix of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = matrix(data$mModMatrix, nrow = 1, byrow = T)
  # calculate fitted value
  iFitted = as.vector(mModMatrix %*% betas)
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  ## likelihood function with posterior theta
  return(mean(dbinom(resp, 1, iFitted, log=F)))
}

ilppd = sum(log(sapply(seq_along(lData$resp), function(x) {
  d = list(resp=lData$resp[x], mModMatrix = lData$mModMatrix[x,])
  lppd(s2, d)
})))

## effective numbers of parameters pWAIC1
pWAIC1 = 2 * (ilppd - eLPD)

iWAIC = -2 * (ilppd - pWAIC1)




####################### lets do this in R, or try to

myloglike.random = function(theta, data){
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  groupIndex = data$groupIndex  ## mapping variable to map each random effect to its respective response variable
  ## parameters to track/estimate
  sigmaRan = exp(theta['sigmaRan']) # random effect scale/sd
  betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
  iGroupsJitter = theta[grep('ran', names(theta))]# random effects jitters for the group deflections
  
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
  lhp = dunif(sigmaRan, 0, 2, log=T)
  lran = sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T))
  lp = sum(dcauchy(betas, 0, 10, log=T))
  # write the likelihood function
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lhp + lran + lp + lik
  return(val)
}

# define a modification of the laplace function from learnbayes
library(numDeriv)
library(optimx)

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

## variance for the random effect 
v = sd(tapply(logit(lData$resp), lData$groupIndex, mean) - mean(logit(lData$resp)))

start = c(sigmaRan=log(v), betas=rep(0, times=ncol(lData$mModMatrix)), 
          ran=rep(0, times=length(unique(lData$groupIndex))))

myloglike.random(start, lData)

fit.4 = laplace(myloglike.random, start, lData)
fit.4 = mylaplace(myloglike.random, start, lData)

data.frame(fixef(fit.3), fit.4$mode[2:6])  ## pretty close

## they both seem to converge close enough, lets try a few other optimisers, this
## may take a while
op = optimx(start, myloglike.random, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)
summary(op) ## nelder-mead seems to converge better

## try with another starting value
v = sd(tapply(logit(lData$resp), lData$groupIndex, mean) - mean(logit(lData$resp)))
d = tapply(logit(lData$resp), lData$groupIndex, mean) - mean(logit(lData$resp))
start = c(sigmaRan=log(v), betas=fit.2$mode, 
          ran=d)

fit.5 = laplace(myloglike.random, start, lData)
fit.5 = mylaplace(myloglike.random, start, lData)

data.frame(fixef(fit.3), fit.4$mode[2:6], fit.5$mode[2:6])  ## pretty close


# ## second version of the function
# # but this time restricting the search space and 
# ## not tracking all the parameters i.e. random effects
# # write a new log posterior function with random effects
# myloglike.random2 = function(theta, data){
#   ## data
#   resp = data$resp # resp
#   mModMatrix = data$mModMatrix
#   groupIndex = data$groupIndex  ## mapping variable to map each random effect to its respective response variable
#   ## we are restricting the random effects deviations to the ones observed in the data
#   ## calculating jitters using the group means subtracted by the grand mean
#   iGroupsJitter = tapply(logit(resp), groupIndex, mean) - mean(logit(resp))
#   ## parameters to track/estimate
#   sigmaRan = exp(theta['sigmaRan']) # random effect scale/sd
#   betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
#   
#   ## random effect jitter for the population intercept
#   # each group contributes a jitter centered on 0
#   # population slope + random jitter
#   ivBetaRand = betas[1] + iGroupsJitter
#   # create a matrix of betas with the new interceptr/unique intercept for each random effect
#   ivIntercept = ivBetaRand[groupIndex] # expand this intercept
#   iFitted = mModMatrix[,2:ncol(mModMatrix)] %*% betas[2:ncol(mModMatrix)]
#   iFitted = as.vector(ivIntercept) + as.vector(iFitted)
#   # using logit link so use inverse logit
#   iFitted = logit.inv(iFitted)
#   ## first level, log priors
#   lhp = dunif(sigmaRan, 0, 2, log=T)
#   lran = sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T))
#   lp = sum(dcauchy(betas, 0, 10, log=T))
#   # write the likelihood function
#   lik = sum(dbinom(resp, 1, iFitted, log=T))
#   val = lhp + lran + lp + lik
#   return(val)
# }
# 
# start = c(sigmaRan=log(v), betas=rep(0, times=ncol(lData$mModMatrix)))
# start = c(sigmaRan=log(v), betas=fit.2$mode)
# 
# myloglike.random2(start, lData)
# 
# fit.5 = laplace(myloglike.random2, start, lData)
# fit.5 = mylaplace(myloglike.random2, start, lData)

#op = optimx(start, myloglike.random, control = list(maximize=T, maxit=500, all.methods=T, usenumDeriv=T), data=lData)

# start = coef(op)[2,]
# 
# 
# op = optim(start, myloglike.random, gr=NULL,
#            control = list(fnscale = -1, maxit=1000), method='SANN', data=lData)
# start = op$par
# 
# 
# start = c(sigmaRan=log(mean(m$sigmaRan)), betas=apply(m$betas, 2, mean))
# start = c(sigmaRan=log(sd(lData$resp)), betas=rep(0, times=ncol(lData$mModMatrix)))
# myloglike.random2(start, lData)
# 
# fit.5 = laplace(myloglike.random2, start, lData)
# 
# 
# ### define a custom laplace function 
# 
# 
# 
# mylaplace = function (logpost, mode, data) 
# {
#   options(warn = -1)
#   fit = optim(mode, logpost, gr = NULL,  
#               control = list(fnscale = -1, maxit=1000), method='CG', data=data)
#   # calculate hessian
#   fit$hessian = (hessian(logpost, fit$par, data=data))
#   colnames(fit$hessian) = names(mode)
#   rownames(fit$hessian) = names(mode)
#   options(warn = 0)
#   mode = fit$par
#   h = -solve(fit$hessian)
#   stuff = list(mode = mode, var = h, converge = fit$convergence == 
#                  0)
#   return(stuff)
# }
# 
# fit.6 = mylaplace(myloglike.random2, start, lData)
# 
# 
# 
# start['mu'] = mean(extract(fit.stan)$mu)
# start['tau'] = log(mean(extract(fit.stan)$tau))
# start
# 
# mylaplace(mylogpost, start, lData)
# 


