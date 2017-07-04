# File: finiteMixtureRegression.R
# Auth: Umar Niazi,
# Date: 4/7/2017
# Desc: Fitting a finite mixture regression model


### generate some test data
set.seed(201503031)

## 2 intercepts
gSigmaInt = 2
gIntercepts = rnorm(n = 2, mean = 0, sd = gSigmaInt)

# mixture probability
iMix = c(0.6, 0.4)

## error for each distribution
gSigmaPop = c(2, 4)
# we have 5 groups X 6 per group = 30 observations
mErrors = matrix(c(rnorm(n = 100, mean = 0, sd = gSigmaPop[1]), 
                   rnorm(n = 100, mean = 0, sd = gSigmaPop[2])), ncol = 2, byrow = F)

## We generate some random predictor and
## generate some data
## predictor variable
gPredictor = runif(n = 100, min = 0, max = 20)
## Population slopes for each distribution.
gSlope = c(1.5, 3)
## Outcome
gResponse1 = gIntercepts[1] + gSlope[1]*gPredictor + mErrors[,1]
gResponse2 = gIntercepts[2] + gSlope[2]*gPredictor + mErrors[,2]
gResponseAggregate = cbind(gResponse1, gResponse2)
gResponseAggregate = sweep(gResponseAggregate, 2, iMix, '*')
gResponseAggregate = rowSums(gResponseAggregate)

dfData = data.frame(resp = gResponseAggregate, pred=gPredictor)

plot(density(dfData$resp))
lines(density(gResponse1))
lines(density(gResponse2))

library(flexmix)

fit.flex = flexmix(resp ~ pred, data=dfData, k=2)
## fitted coefficients
parameters(fit.flex)
## predicted values, both return the same values
p = predict(fit.flex, newdata=dfData)
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
p2 = predict(fit.flex, newdata=dfData, aggregate=T)
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
m = model.matrix(resp ~ pred, data=dfData)
## prepare data
lStanData = list(Ntotal=nrow(dfData), y=dfData$resp, iMixtures=2, Ncol=ncol(m), X=m)

## give initial values
initf = function(chain_id = 1) {
  list(sigma = c(11, 12), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, init=initf, cores=4, pars=c('sigma', 'iMixWeights', 'betasMix1',
                                                                                              'betasMix2', 'mu'))
print(fit.stan, digi=3)
traceplot(fit.stan)

## calculate fitted values
m = extract(fit.stan)

mBetas.Comp1 = m$betas[,,1]
mBetas.Comp2 = m$betas[,,2]

f1 = lStanData$X %*% colMeans(mBetas.Comp1)
f2 = lStanData$X %*% colMeans(mBetas.Comp2)


