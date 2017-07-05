# File: finiteMixtureRegression.R
# Auth: Umar Niazi,
# Date: 4/7/2017
# Desc: Fitting a finite mixture regression model


library(flexmix)
data("NPreg")
str(NPreg)

### how is the response variable distributed
yresp = density(NPreg$yn)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Response variable', ylab='scaled density')
########### this part of the script is looking at if the mixture model is appropriate
############################### fit a model using stan to estimate mixture parameters
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='BayesMixtureModel/fitNormalMixture.stan')

## take a subset of the data
lStanData = list(Ntotal=length(NPreg$yn), y=NPreg$yn, iMixtures=2)

## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(10, 30), sigma = c(11, 112), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=4000, chains=4, init=initf, cores=4)
print(fit.stan, digi=3)
traceplot(fit.stan)

## check if labelling degeneracy has occured
## see here: http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
params1 = as.data.frame(extract(fit.stan, permuted=FALSE)[,1,])
params2 = as.data.frame(extract(fit.stan, permuted=FALSE)[,2,])
params3 = as.data.frame(extract(fit.stan, permuted=FALSE)[,3,])
params4 = as.data.frame(extract(fit.stan, permuted=FALSE)[,4,])

## check if the means from different chains overlap
## Labeling Degeneracy by Enforcing an Ordering
par(mfrow=c(2,2))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
plot(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)
plot(params3$`mu[1]`, params3$`mu[2]`, pch=20, col=4)
plot(params4$`mu[1]`, params4$`mu[2]`, pch=20, col=5)

par(mfrow=c(1,1))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2, xlim=c(3, 28), ylim=c(28, 42))
points(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)
points(params3$`mu[1]`, params3$`mu[2]`, pch=20, col=4)
points(params4$`mu[1]`, params4$`mu[2]`, pch=20, col=5)

############# extract the mcmc sample values from stan
mStan = do.call(cbind, extract(fit.stan))
mStan = mStan[,-(ncol(mStan))]
colnames(mStan) = c('mu1', 'mu2', 'sigma1', 'sigma2', 'mix1', 'mix2')
dim(mStan)
## get a sample for this distribution
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(NPreg$yn), ncol=200)

for (i in 1:200){
  p = sample(1:nrow(mStan), size = 1)
  mix = mean(mStan[,'mix1'])
  ## this will take a sample from a normal mixture distribution
  sam = function() {
    ind = rbinom(1, 1, prob = mix)
    return(ind * rnorm(1, mStan[p, 'mu1'], mStan[p, 'sigma1']) + 
             (1-ind) * rnorm(1, mStan[p, 'mu2'], mStan[p, 'sigma2']))
  }
  mDraws[,i] = replicate(length(NPreg$yn), sam())
}

mDraws.normMix = mDraws

plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2)
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)

## perform some checks
## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
# The procedure for carrying out a posterior predictive model check requires specifying a test
# quantity, T (y) or T (y, θ), and an appropriate predictive distribution for the replications
# y rep [Gelman 2008]
## variance
T1_var = function(Y) return(var(Y))

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 


t1 = apply(mDraws, 2, T1_var)
ivTestQuantities = getPValue(t1, var(lStanData$y))

t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lStanData$y)
ivTestQuantities = c(ivTestQuantities, getPValue(t1, t2))

t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lStanData$y)
ivTestQuantities = c(ivTestQuantities, getPValue(t1, t2))

t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lStanData$y)
ivTestQuantities = c(ivTestQuantities, getPValue(t1, t2))

names(ivTestQuantities) = c('variance', 'minimum', 'maximum', 'mean')

ivTestQuantities

## the mixture model with 2 components appears to be appropriate, lets start fitting
###################################################################################

fit.flex = flexmix(yn ~ x, data=NPreg, k=2)
summary(fit.flex)
## fitted coefficients
parameters(fit.flex)
## predicted values, both return the same values
p = predict(fit.flex, newdata=NPreg)
p = do.call(cbind, p)
f = fitted(fit.flex)

## they will however return 2 components, 
head(f)

## if we want to aggregate these values, we need a weighted average, based on
## the weights assigned to the mixture components
## following the source code from the flexmix package
## https://rdrr.io/cran/flexmix/src/R/flexmix.R
# if (aggregate) {
#   prior_weights <- prior(object, newdata)
#   z <- lapply(x, function(z) matrix(rowSums(do.call("cbind", z) * prior_weights), nrow = nrow(z[[1]])))
# }
## lets check
p2 = predict(fit.flex, newdata=NPreg, aggregate=T)
p2 = do.call(cbind, p2)
head(p2)

## take a simple averga
head(rowMeans(p)) ## not the same
pr = prior(fit.flex)
p.agg = sweep(p, 2, pr, '*')
p.agg = rowSums(p.agg)
head(p.agg)
identical(as.numeric(p.agg), as.numeric(p2))


#plot(density(NPreg$yn))
############### fit a similar model using stan and MCMC approach
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='BayesMixtureModel/fitNormalMixtureRegression.stan')
m = model.matrix(yn ~ x, data=NPreg)
## prepare data
lStanData = list(Ntotal=nrow(NPreg), y=NPreg$yn, iMixtures=2, Ncol=ncol(m), X=m)

## give initial values
initf = function(chain_id = 1) {
  list(sigma = c(10, 10), iMixWeights=c(0.5, 0.5), mu=c(10, 40))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=4000, chains=1, init=initf, cores=4, pars=c('sigma', 'iMixWeights', 'betasMix1',
                                                                                              'betasMix2', 'mu'))
print(fit.stan, digi=3)
traceplot(fit.stan)
parameters(fit.flex)

### fit a second model with more restrictive priors 
stanDso2 = rstan::stan_model(file='BayesMixtureModel/fitNormalMixtureRegression2.stan')
m = model.matrix(yn ~ x, data=NPreg)
## prepare data
lStanData = list(Ntotal=nrow(NPreg), y=NPreg$yn, iMixtures=2, Ncol=ncol(m), X=m)

## give initial values
initf = function(chain_id = 1) {
  list(sigma = c(10, 10), iMixWeights=c(0.5, 0.5), mu=c(0, 35))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan2 = sampling(stanDso2, data=lStanData, iter=10000, chains=1, init=initf, cores=4, pars=c('sigma', 'iMixWeights', 'betasMix1',
                                                                                              'betasMix2', 'mu'))
print(fit.stan2)

## get fitted values
m = extract(fit.stan2)
names(m)
intercepts = apply(m$mu, 2, mean)
betas = c(mean(m$betasMix1), mean(m$betasMix2))
iMixWeights = c(0.44, 0.56)

iPred1 = lStanData$X %*% c(intercepts[1], betas[1])
iPred2 = lStanData$X %*% c(intercepts[2], betas[2])

## compare with fit.flex
head(f)
head(cbind(iPred1, iPred2))

## get aggregate
iAggregate = cbind(iPred1, iPred2)
iAggregate = sweep(iAggregate, 2, iMixWeights, '*')
iAggregate = rowSums(iAggregate)

head(data.frame(iAggregate, p.agg))

dfPredicted = data.frame(stan=iAggregate, flexmix=p.agg)

## calculate MSE
mean((dfPredicted$flexmix - NPreg$yn)^2)
mean((dfPredicted$stan - NPreg$yn)^2)
## both pretty close



