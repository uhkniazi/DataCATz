# File: betaMixtureCoinBias.R
# Auth: Umar Niazi,
# Date: 18/4/2017
# Desc:

### define two models that explain the data
## model m1
## it is a beta distribution with a weight on the left side
g1 = function(theta) dbeta(theta, 3.5, 8.5)

## model m2
## it is a beta distribution with a weight on the right side
g2 = function(theta) dbeta(theta, 8.5, 3.5)

## define an array that represents number of models in our parameter space
## each index has a prior weight/probability of being selected
## this can be thought of coming from a categorical distribution 
mix.prior = c(50, 50)
mix.prior = mix.prior/sum(mix.prior)

## define a joint prior space for these 2 models
## (mix.prior[1] AND g1) OR (mix.prior[2] AND g2)
g.mix = function(theta) mix.prior[1]*g1(theta) + mix.prior[2]*g2(theta)

## some plots to look at the shapes of the priors
p.old = par(mfrow=c(2,1))
curve(g1, 0, 1, main='Model m1')
curve(g2, 0, 1, main='Model m2')
par(p.old)
curve(g.mix, 0, 1, main='Joint Prior Mixture of 2 models')

## we flip a coin, perform a binary experiment 
## our prior beliefs are that the experiment tends to follow two extreme models
## we will either see more 1s or more 0s, depending on which model (m1 or m2) is 
## more representative of the data
data = c(success=6, fail=9-6)

## the posterior is proportional to Likelihood * Prior
## P(Data | Theta, Model) X P(Theta, Model)

## define the Likelihood function, both models share the same likelihood functional form
lik = function(data, theta) dbinom(data['success'], sum(data), theta)

## since this is a conjugate analysis, as prior is beta distributed and likelihood is binomial
## the posterior can be derived analytically 
## using a similar analogy as before when describing priors for each model
## we can define a posterior for each of the 2 models
g1.post = function(theta, data) dbeta(theta, 3.5+data['success'], 8.5+data['fail'])
g2.post = function(theta, data) dbeta(theta, 8.5+data['success'], 3.5+data['fail'])

## however in order to calculate a mixture probability, i.e. what is the probabilty
## each model m1 or m2 being selected from the categorical distribution after we 
## see the data i.e. posterior for m1 and m2
## integrate across the theta parameter for the selected model
## see page P-266-268 Equations 10.1 to 10.3, Doing Bayesian Data Analysis [REF]
## P(Data, Model[1]) = P(Data, Model[1])
## P(Model[1] | Data) X P(Data) = P(Data | Model[1]) X P(Model[1])
## P(Model[1] | Data) = P(Data | Model[1]) X P(Model[1]) / P(Data)
## where P(Data) can be expanded using summation across all models
## P(Data) = Sum for each Model P(Data, Model[i])
## P(Data) = Sum for each Model P(Data | Model[i]) X P(Model[i])

## we need to calculate a few things
## P(Data | Model[i])
## this is the prior predictive probability for the data given the selected model
## P(Data, Theta) = P(Data | Theta) X P(Theta)
## P(Data) = P(Data | Theta) X P(Theta) / P(Theta | Data)
## Prior predictive probability for Data = Likelihood X Prior / Posterior
## for model 1
data.prior.g1 = function(data, theta){
  ret = log(lik(data, theta)) + log(g1(theta)) - log(g1.post(theta, data))
  return(exp(ret))
}
## for model 2
data.prior.g2 = function(data, theta){
  ret = log(lik(data, theta)) + log(g2(theta)) - log(g2.post(theta, data))
  return(exp(ret))
}

## P(Data | Model) for each model should be the same for any value of theta
## you can use that as a sanity check
th = seq(0.01, 0.99, by=0.01)
data.g1 = mean(data.prior.g1(data, th))
data.g2 = mean(data.prior.g2(data, th))
## we have P(Data | Model[1]) - prior predictive probability of data for parameters in model 1
## P(Model[1]) - prior mixture probability for model 1
## P(Data | Model[2]); P(Model[2])
## P(Data) = (P(Data | Model[1]) X P(Model[1])) + (P(Data | Model[2]) X P(Model[2])
## we can calculate the posterior for Model 1
## P(Model[1] | Data) = P(Data | Model[1]) X P(Model[1]) / P(Data)
mix.post = data.g1 * mix.prior[1] / (data.g1 * mix.prior[1] + data.g2 * mix.prior[2])
mix.post = c(mix.post, 1-mix.post)

## Bayes factor for the ratios of posterior predictive distribution
## of the 2 models
## P(Data | Model[1]) / P(Data | Model[2])
BF = data.g1 / data.g2
## OR posterior odds for the models / prior odds for the 2 models
mix.post[1]/mix.post[2]/(mix.prior[1]/mix.prior[2])


## (mix.post[1] AND g1.post) OR (mix.post[2] AND g2.post)
g.mix.post = function(theta, data) mix.post[1]*g1.post(theta, data) + mix.post[2]*g2.post(theta, data)

## make some curves to look at the theta
par(mfrow=c(2,1))
# use the grid of thetas to draw the curves
th = seq(0.01, 0.99, by=0.01)
plot(th, g1.post(th, data), type='l', main='Posterior Model m1')
plot(th, g2.post(th, data), type='l', main='Posterior Model m2')
par(p.old)
plot(th, g.mix.post(th, data), type='l', main='Joint Posterior mixture of 2 Models')
lines(th, g.mix(th), lty=2)
legend('topleft', legend = c('prior', 'posterior'), lty=c(2, 1))

## approximate the posterior theta on the grid
p.th = g.mix.post(th, data)
p.th = p.th/sum(p.th)
th.sam = sample(th, 10000, replace = T, prob=p.th)
th.sam = th.sam + runif(10000, 0.01/2 * -1, 0.01/2)
summary(th.sam); sd(th.sam)

## plot the sample histogram and density together
hist(th.sam, prob=T)
lines(th, g.mix.post(th, data))

############################
## try calculating the same parameters using a single function and optimizer
library(LearnBayes)
#library(DirichletReg)
library(car)
#logit = function(p) log(p/(1-p))
logit.inv = function(p) {exp(p)/(exp(p)+1) }

mylogpost = function(theta, data){
  ## theta contains parameters we wish to track
  th = logit.inv(theta['theta'])
  m1 = logit.inv(theta['m1'])
  m2 = 1-m1
  suc = data['success']
  fail = data['fail']

  # define likelihood function
  lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))

  # define mixture distribution
  g1 = function(t) dbeta(t, 3.5, 8.5)
  g2 = function(t) dbeta(t, 8.5, 3.5)
  g.mix = log(m1*g1(th) + m2*g2(th))
  ms = matrix(c(m1, m2), nrow = 1)
  # calculate log posterior
  val = lf(suc, fail, th) + g.mix + dunif(m1, log=T)#ddirichlet(ms, c(1, 1), log=T)
  return(val)
}

start = c(theta=logit(0.5), m1=logit(0.9))
data = c(success=6, fail=3)

mylogpost(start, data)

fit = laplace(mylogpost, start, data)
logit.inv(fit$mode)

mylogpost_m1 = function(theta, data){
  ## theta contains parameters we wish to track
  th = logit.inv(theta['theta'])
  suc = data['success']
  fail = data['fail']
  
  # define likelihood function
  lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
  
  # calculate log posterior
  val = lf(suc, fail, th) + dbeta(th, 3.5, 8.5, log = T)
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
  val = lf(suc, fail, th) + dbeta(th, 8.5, 3.5, log = T) 
  return(val)
}

start = c(theta=logit(0.5))
data = c(success=6, fail=3)

mylogpost_m1(start, data)
mylogpost_m2(start, data)


fit_m1 = laplace(mylogpost_m1, start, data)
fit_m2 = laplace(mylogpost_m2, start, data)

fit_m1
fit_m2

logit.inv(fit_m1$mode)
logit.inv(fit_m2$mode)

exp(fit_m1$int)
exp(fit_m2$int)

exp(fit_m1$int) / exp(fit_m2$int)

exp(fit_m1$int) * mix.prior[1] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2])
