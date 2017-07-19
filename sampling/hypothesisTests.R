# File: hypothesisTests.R
# Auth: uhkniazi
# Date: 19/07/2017
# Desc: testing binary and multiple hypotheses


### define three machines
## model m1 - machine makes less defective widgets
m1 = function(th) dbeta(th, 1, 5, log = T)

## model m2 - machine makes between average defective widgets
m2 = function(th) dbeta(th, 5, 5, log = T)

## model m3 - machine makes almost everything bad
m3 = function(th) dbeta(th, 5, 1, log = T)

## define an array that represents number of models in our parameter space
## each index has a prior weight/probability of being selected
## this can be thought of coming from a categorical distribution 
mix.prior = c(m1=10/11 *(1-1e-6),m2= 1/11 *(1-1e-6),m3= 1e-6)

## define a joint prior space for these 3 models
g.mix = function(theta) {
  log(mix.prior['m1']*exp(m1(theta)) + mix.prior['m2']*exp(m2(theta)) +
    mix.prior['m3']*exp(m3(theta)))}

## some plots to look at the shapes of the priors
curve(m1, 0, 1, main='machine 1')
curve(m2, 0, 1, main='machine 2')
curve(m3, 0, 1, main='machine 3')
par(p.old)
curve(g.mix, 0, 1, main='Joint Prior Mixture of 3 models')

# we observe data from a machine, how many bad do we see
data = rbinom(20, 1, 1/6)
b = sum(data)
data = c(fail=length(data)-b, success=b)

############################
## try using log posterior function
library(LearnBayes)
library(car)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

mylogpost_m1 = function(theta, data){
  ## theta contains parameters we wish to track
  th = logit.inv(theta['theta'])
  suc = data['success']
  fail = data['fail']
  
  # define likelihood function
  lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
  
  if (th >=0.17) return(-Inf)
  # calculate log posterior
  val = lf(suc, fail, th) + m1(th) #dunif(th, 0, 0.17, log=T)
  return(val)
}

mylogpost_m2 = function(theta, data){
  ## theta contains parameters we wish to track
  th = logit.inv(theta['theta'])
  suc = data['success']
  fail = data['fail']
  
  # define likelihood function
  lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
  
  # calculate log posterior
  val = lf(suc, fail, th) + m2(th)
  return(val)
}

mylogpost_m3 = function(theta, data){
  ## theta contains parameters we wish to track
  th = logit.inv(theta['theta'])
  suc = data['success']
  fail = data['fail']
  
  # define likelihood function
  lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
  
  # calculate log posterior
  val = lf(suc, fail, th) + m3(th)
  return(val)
}

## choose starting values
# failure rate for first machine
start = c(theta=logit(1/6))
mylogpost_m1(start, data)
fit_m1 = laplace(mylogpost_m1, start, data)

# failure rate for second machine
start = c(theta=logit(1/3))
mylogpost_m2(start, data)
fit_m2 = laplace(mylogpost_m2, start, data)

# failure rate for 3rd machine
start = c(theta=logit(99/100))
mylogpost_m3(start, data)
fit_m3 = laplace(mylogpost_m3, start, data)

# values for theta for each model
logit.inv(fit_m1$mode)
logit.inv(fit_m2$mode)
logit.inv(fit_m3$mode)
# approximate posterior predictive distribution
exp(fit_m1$int)
exp(fit_m2$int)
exp(fit_m3$int)

BFm1 = exp(fit_m1$int) / exp(fit_m2$int)
BFm3 = exp(fit_m3$int) / exp(fit_m2$int)

# posterior prob for the models
mix.post = c(m1=0, m2=0, m3=0)
mix.post[1] = exp(fit_m1$int) * mix.prior[1] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2]
                                                + exp(fit_m3$int) * mix.prior[3])

mix.post[2] = exp(fit_m2$int) * mix.prior[2] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2]
                                                + exp(fit_m3$int) * mix.prior[3])

mix.post[3] = exp(fit_m3$int) * mix.prior[3] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2]
                                                + exp(fit_m3$int) * mix.prior[3])

mix.post