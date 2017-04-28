# File: eightSchools.R
# Auth: uhkniazi
# Date: 28/04/2017
# Desc: Fitting a hierarchical model to the 8 schools data set


library(LearnBayes)


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
  th_grps = theta[-c(1,2)]
  
  ## likelihood function
  lf = function(dat, m, s){
    return(dnorm(dat, m, s, log=T))
  }
  val = sum(lf(y_bar, th_grps, sig_group))
  val = val + sum(dnorm(th_grps, mu, tau, log=T)) + dunif(tau, 1, 1e+3, log = T)
  return(val)
}

# define starting values
start = c('mu'=10, 'tau'=log(10), 'th_grps'=rep(0, times=8))
lData = list('y'=y, 'sig'=sig)

optim(start, mylogpost, gr = NULL, control = list(fnscale = -1, maxit=10000), data=lData, method='SANN')