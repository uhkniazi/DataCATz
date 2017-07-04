# File: finiteMixtureRegression.R
# Auth: Umar Niazi,
# Date: 4/7/2017
# Desc: Fitting a finite mixture regression model

library(flexmix)

data("NPreg")
str(NPreg)

fit.flex = flexmix(yn ~ x, data=NPreg, k=2)
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
  list(sigma = c(12, 12), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, init=initf, cores=4, pars=c('sigma', 'iMixWeights', 'betas'))
print(fit.stan, digi=3)
traceplot(fit.stan)

## calculate fitted values
m = extract(fit.stan)

mBetas.Comp1 = m$betas[,,1]
mBetas.Comp2 = m$betas[,,2]

f1 = lStanData$X %*% colMeans(mBetas.Comp1)
f2 = lStanData$X %*% colMeans(mBetas.Comp2)


