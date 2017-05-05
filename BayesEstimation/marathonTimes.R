# Name: marathonTimes.R
# Auth: uhkniazi
# Date: 04/05/2017
# Desc: various parameter estimation methods using examples from a couple of books. 
#       see the main blog description for the references and details.

library(LearnBayes)

attach(marathontimes)

str(marathontimes)
summary(marathontimes)

### we will try a few methods for parameter estimation.
### first we will use a conjugate analysis and calculate the parameters
### for the data analytically as shown in Gelman 2008 - see page 62 
### also see page 64 of bayesian computation with R
## We are not using any prior information, but to make this function more general 
## I am adding some prior information with 0s
## look at page 68 of Bayesian Data Analysis (Gelman) for formula
sim.post = function(dat.grp){
  # prior variance
  sigma.0 = 0
  # prior observations
  k.0 = 0
  # prior degrees of freedom
  v.0 = k.0 - 1
  # prior mean
  mu.0 = 0
  # calculate conjugate posterior
  # number of observations in the data
  n = length(dat.grp)
  # prior observations + data points = posterior number of observations
  k.n = k.0 + n
  # posterior degrees of freedom
  v.n = v.0 + n
  # mean and sd for the data
  y.bar = mean(dat.grp)
  s = sd(dat.grp)
  # posterior mean
  mu.n = (k.0/k.n * mu.0) + (n/k.n * y.bar)
  # posterior var
  sigma.n = (( v.0*sigma.0 ) + ( (n-1)*(s^2) ) + ( (k.0*n/k.n)*((y.bar-mu.0)^2) )) / v.n
  #post.scale = ((prior.dof * prior.scale) + (var(dat.grp) * (length(dat.grp) - 1))) / post.dof
  ## simulate posterior variance and conditional posterior mean
  sigma = (sigma.n * v.n)/rchisq(1000, v.n)
  mu = rnorm(1000, mu.n, sqrt(sigma)/sqrt(k.n))
  return(list(mu=mu, var=sigma))
}

## simulate the posterior variance and get estimates for mean and variance
## the results are very approximately equal to what we see in the book.
sim.1 = sim.post(time)
quantile(sim.1$mu, c(0.025, 0.975))
quantile(sqrt(sim.1$var), c(0.025, 0.975))
quantile(sim.1$var, c(0.025, 0.975))
#################################### Repeat the same analysis but this time
########################### on using a grid based approximation

## define a discrete grid of parameters
Mu = seq(250, 302, by=0.5)
Sigma2 = seq(500, 9000, by=10)

## define a log posterior function
lp = function(theta, data){
  s = theta[2]
  m = theta[1]
  d = data$time # data from marathon times
  log.lik = sum(dnorm(d, m, sqrt(s), log=T))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=Mu[1], 'sigma2'=Sigma2[1])
lData = list('time'=time)
lp(start, lData)
# draw the contour plot 
mycontour(lp, c(220, 330, 500, 9000), lData)
# how do the simulated values earlier fit in this contour plot
points(sim.1$mu, sim.1$var, pch=20)

## calculate the log posterior at each grid point 
lm = matrix(NA, nrow=length(Mu), ncol = length(Sigma2))
for (r in seq_along(Mu)){
  for (c in seq_along(Sigma2)){
    s = c('mu'=Mu[r], 'sigma2'=Sigma2[c])
    lm[r, c] = lp(s, lData)
  }
}

# convert from log posterior 
# subtract max before exponentiating
lm = lm-max(lm)
lm = exp(lm)
# set total prob of grid to 1
lm = lm/sum(colSums(lm))

## take samples from this distribution
pSigma2 = colSums(lm)
pMu = rowSums(lm)

## add a uniform jitter centered at zero and half of grid spacing to make distribution continuous
sig2Sample = runif(1000, -1*10/2, 10/2) +  sample(size = 1000, Sigma2, replace = T, prob = pSigma2)
# for each sampled value of sig2Sample, get the posterior mu, conditioned on sig2Sample
muSample = rnorm(1000, mean(time), sqrt(sig2Sample)/sqrt(length(time)))
## how do they fit on the contour plot
points(muSample, sig2Sample, pch=20, col=2)
## pretty close
quantile(muSample, c(0.025, 0.975))
quantile(sqrt(sig2Sample), c(0.025, 0.975))
quantile(sig2Sample, c(0.025, 0.975))

# take a sample directly from the grid for mean
muSample2 = runif(1000, -1*0.5/2, 0.5/2) + sample(size = 1000, Mu, replace = T, prob = pMu)
quantile(muSample2, c(0.025, 0.975))

################################################################
#### try an optimization based approach
## define a log posterior function
lp2 = function(theta, data){
  # we define the sigma on a log scale as optimizers work better
  # if scale parameters are well behaved
  s = exp(theta[2])
  m = theta[1]
  d = data$time # data from marathon times
  log.lik = sum(dnorm(d, m, s, log=T))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=mean(time), 'sigma'=log(var(time)))
lData = list('time'=time)
lp2(start, lData)

op = optim(start, lp2, control = list(fnscale = -1), data=lData)
op$par
exp(op$par[2])

## try the laplace function from LearnBayes
fit = laplace(lp2, start, lData)
fit
se = sqrt(diag(fit$var))
se
fit$mode+1.96*se
fit$mode-1.96*se

## lets take a sample from this
sigSample.op = rnorm(1000, fit$mode['sigma'], se['sigma'])
muSample.op = rnorm(1000, mean(time), exp(sigSample.op)/sqrt(length(time)))
muSample2.op = rnorm(1000, fit$mode['mu'], se['mu'])

########################## estimation using stan
library(rstan)
stanDso = rstan::stan_model(file='BayesEstimation/marathonTimes.stan')

lStanData = list(Ntotal=length(time), y=time)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('mu', 'sig'))
print(fit.stan)
m = extract(fit.stan)
muStan = m$mu
sigStan = m$sig

##### lets plot all the various approximations together
mycontour(lp, c(220, 330, 500, 9000), lData, xlab='Posterior Mean', ylab='Posterior Variance')
# how do the simulated values earlier fit in this contour plot
points(muStan, sigStan^2, pch=20, col='cyan')
points(sim.1$mu, sim.1$var, pch=20, col='black')
points(muSample, sig2Sample, pch=20, col='red')
points(muSample.op, exp(sigSample.op)^2, pch=20, col='blue')

par(mfrow=c(1,2))
plot(density(sim.1$mu), col=1, type='l', main='Posterior Mean Density')
lines(density(muSample), col='red')
lines(density(muSample.op), col='blue')
lines(density(muSample2.op), col='green')
lines(density(muSample2), col='darkgreen')
lines(density(muStan), col='cyan')

plot(density(sim.1$var), col=1, type='l', main='Posterior Variance Density')
lines(density(sig2Sample), col='red')
lines(density(exp(sigSample.op)^2), col='blue')
lines(density(sigStan^2), col='cyan')

