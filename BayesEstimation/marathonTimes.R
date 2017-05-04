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
  # we define the sigma on a log scale as optimizers work better
  # if scale parameters are well behaved
  s = exp(theta['sigma2'])
  m = theta['mu']
  d = data$time # data from marathon times
  log.lik = sum(dnorm(d, m, sqrt(s), log=T))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=Mu[1], 'sigma2'=log(Sigma2[1]))
lData = list('time'=time)
lp(start, lData)

op = optim(start, lf, control = list(fnscale = -1), data=X, hessian = T)
exp(op$par[2])

## calculate the log posterior at each grid point 
lm = matrix(NA, nrow=length(Mu), ncol = length(Sigma2))
for (r in seq_along(Mu)){
  for (c in seq_along(Sigma2)){
    s = c('mu'=Mu[r], 'sigma2'=log(Sigma2[c]))
    lm[r, c] = lp(s, lData)
  }
}

# convert from log posterior 
# subtract max before exponentiating
lm = lm-max(lm)
lm = exp(lm)
# set total prob of grid to 1
lm = lm/sum(colSums(lm))

## plot contour
contour(Mu, Sigma, lm, xlim=c(9, 10.5), ylim=c(2.2, 3.5))

## take samples from this distribution
pSigma = colSums(lm)
pMu = rowSums(lm)

sigSample = runif(1000, -0.05, 0.05) +  sample(size = 1000, Sigma, replace = T, prob = pSigma)
muSample = runif(1000, -0.05, 0.05) + sample(size = 1000, Mu, replace = T, prob = pMu)
contour(Mu, Sigma, (lm)/max(lm),  levels = c(0, 0.5, 0.3, 0.1, 0.01, 0.001, 0.0001))
points(muSample, sigSample, col=1, pch=20, cex=0.5)
