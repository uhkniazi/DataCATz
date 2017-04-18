# File: betaMixtureCoinBias.R
# Auth: Umar Niazi,
# Date: 18/4/2017
# Desc:

### define two models that explain the data
## model m1
## it is a beta distribution with a weight on the left side
g1 = function(theta) dbeta(theta, 6, 14)

## model m2
## it is a beta distribution with a weight on the right side
g2 = function(theta) dbeta(theta, 14, 6)

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
data = c(success=7, fail=3)

## the posterior is proportional to Likelihood * Prior
## P(Data | Theta, Model) X P(Theta, Model)

## define the Likelihood function, both models share the same likelihood functional form
lik = function(data, theta) dbinom(data['success'], sum(data), theta)

## since this is a conjugate analysis, as prior is beta distributed and likelihood is binomial
## the posterior can be derived analytically 
## using a similar analogy as before when describing priors for each model
## we can define a posterior for each of the 2 models
g1.post = function(theta, data) dbeta(theta, 6+data['success'], 14+data['fail'])
g2.post = function(theta, data) dbeta(theta, 14+data['success'], 6+data['fail'])

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
  ret = lik(data, theta) * g1(theta) / g1.post(theta, data)
  return(ret)
}
## for model 2
data.prior.g2 = function(data, theta){
  ret = lik(data, theta) * g2(theta) / g2.post(theta, data)
  return(ret)
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




