# File: finiteMixtureRegression.R
# Auth: Umar Niazi,
# Date: 4/7/2017
# Desc: Fitting a finite mixture regression model


library(flexmix)
data("NPreg")
str(NPreg)

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


plot(density(NPreg$yn))
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
  list(sigma = c(10, 10), iMixWeights=c(0.5, 0.5))
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
head(iPred1)
head(iPred2)

## get aggregate
iAggregate = cbind(iPred1, iPred2)
iAggregate = sweep(iAggregate, 2, iMixWeights, '*')
iAggregate = rowSums(iAggregate)

head(data.frame(iAggregate, p.agg))

dfPredicted = data.frame(iAggregate, p.agg)

## calculate MSE
mean((NPreg$yn - dfPredicted$p.agg)^2)
mean((NPreg$yn - dfPredicted$iAggregate)^2)

## both pretty close



