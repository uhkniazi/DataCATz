# Name: estimateBBparam.R
# Auth: uhkniazi@gmail.com
# Date: 7/3/2023
# Desc: Load a data set and estimate the parameters for the beta-binomial distribution

# instead of using Stan directly we are using the rethinking library, however you can see 
# the stan code by calling rethinking::stancode()

library(rethinking)

# generate some dummy data
#y = rbetabinom(n = 200, size = 1000, prob = 0.1, theta = 10)



#d = data.frame(y, t=1000)

dat = list( A= , N=)

m1 <- ulam(
  alist(
    # likelihood
    A ~ dbetabinom( N , pbar , rDispersion ),
    # proportion estimated on a logit scale 
    logit(pbar) <- rProportion,
    rProportion ~ normal(-10, 3),
    rDispersion ~ cauchy(0, 1)
  ), data=dat , chains=4, cores = 4, constraints=list(rDispersion="lower=1,upper=2000", rProportion='upper=0')  )

# check traceplot for chains
plot(m1)
traceplot_ulam(m1)
# check Rhat - should not be much larger than 1
precis(m1)
lPosterior = extract.samples(m1)
precis(lPosterior)

par(mfrow=c(1,1))
plot(density(lPosterior$rProportion), xlab='')
plot(density(lPosterior$rDispersion), xlab='', main='m1')
plot(density(plogis(lPosterior$rProportion)), xlab='', main='m1')

## second parameterisation
m2 <- ulam(
  alist(
    # likelihood
    A ~ dbetabinom( N , pbar , rDispersion ),
    # proportion estimated on a logit scale 
    logit(pbar) <- rProportion,
    rDispersion <- (1-rRhop)/rRhop,
    logit(rRhop) <- rRho, 
    rProportion ~ dunif(-14, 0),
    rRho ~ dunif(-14, 0)
  ), data=dat , chains=4, cores = 4, iter = 1000)

traceplot_ulam(m2)
# check Rhat - should not be much larger than 1
precis(m2)
lPosterior2 = extract.samples(m2)
precis(lPosterior2)

par(mfrow=c(1,1))
plot(density(lPosterior2$rProportion))
plot(density(plogis(lPosterior2$rRho)), xlab='')
plot(density((1-plogis(lPosterior2$rRho))/plogis(lPosterior2$rRho)), xlab='', main='m2')
plot(density(plogis(lPosterior2$rProportion)), xlab='', main='m2')

