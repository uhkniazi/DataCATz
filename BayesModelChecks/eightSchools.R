# File: eightSchools.R
# Auth: uhkniazi
# Date: 28/04/2017
# Desc: Fitting a hierarchical model to the 8 schools data set


library(LearnBayes)
library(numDeriv)

## The 8 Schools data set described in Section 5.5 of Gelman et al (2003), which studied coaching effects from eight schools. 

y = c(28,  8, -3,  7, -1,  1, 18, 12)
sig = c(15, 10, 16, 11,  9, 11, 10, 18)

## define the log posterior function
mylogpost = function(theta, data){
  ## parameters and data that are known
  # mean for each group, i.e. y_bar
  y_bar = data$y
  # standard deviation for each group
  sig_group = data$sig
  
  ## parameters to track
  ## prior for the hierarchical hyperparameters
  mu = theta['mu']
  tau = exp(theta['tau'])
  eta = theta[-c(1,2)]
  
  ## likelihood function
  lf = function(dat, m, s){
    return(dnorm(dat, m, s, log=T))
  }
  th_grps = mu + tau * eta
  val = sum(lf(y_bar, th_grps, sig_group))
  val = val + sum(dnorm(eta, 0, 1, log=T)) + dunif(tau, 0, 1e+3, log = T)
  return(val)
}

# define starting values
start = c('mu'=mean(y), 'tau'=log(mean(sig)), 'eta'=rep(0, times=8))
lData = list('y'=y, 'sig'=sig)

op = optim(start, mylogpost, gr = NULL, control = list(fnscale = -1, maxit=10000), data=lData, method='SANN')
op$par


mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), data=data)
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

mylaplace(mylogpost, start, lData)



############################# test with stan
library(rstan)
stanDso = rstan::stan_model(file='BayesModelChecks/eightSchools.stan')

lStanData = list(Ntotal=8, y=y, sig=sig)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('mu', 'tau', 'eta', 'theta'))
print(fit.stan)

start['mu'] = mean(extract(fit.stan)$mu)
start['tau'] = log(mean(extract(fit.stan)$tau))
start

mylaplace(mylogpost, start, lData)

