# Name: newcombSpeedOfLight.R
# Auth: uhkniazi
# Date: 11/05/2017
# Desc: posterior predictive checks for the model using Newcomb's speed of light data


library(LearnBayes)
## the example is taken from Simon Newcomb's experiment to measure the speed of light
## we follow Gelman [2008] Page 67 for the analysis.

dfData = read.csv('BayesModelChecks/newcombMeasurements.csv', header=T)
str(dfData)

ivTime = dfData$Transformed
summary(ivTime)
sd(ivTime)

## define a log posterior function
lp = function(theta, data){
  # we define the sigma on a log scale as optimizers work better
  # if scale parameters are well behaved
  s = exp(theta[2])
  m = theta[1]
  d = data$vector # observed data vector
  log.lik = sum(dnorm(d, m, s, log=T))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=mean(ivTime), 'sigma'=log(sd(ivTime)))
lData = list('vector'=ivTime)
lp(start, lData)

op = optim(start, lp, control = list(fnscale = -1), data=lData)
op$par
exp(op$par[2])

## try the laplace function from LearnBayes
fit = laplace(lp, start, lData)
fit
se = sqrt(diag(fit$var))
se
fit$mode+1.96*se
fit$mode-1.96*se

## lets take a sample from this
sigSample.op = rnorm(1000, fit$mode['sigma'], se['sigma'])
muSample.op = rnorm(1000, mean(lData$vector), exp(sigSample.op)/sqrt(length(lData$vector)))
# second way of taking the sample
muSample2.op = rnorm(1000, fit$mode['mu'], se['mu'])

#### the posterior interval reported on P 67 Gelman [2008] is
## y ± 1.997s/ 66 = [23.6, 28.8]
## compare with our intervals
fit$mode['mu']-1.96*se['mu']; fit$mode['mu']+1.96*se['mu']
quantile(muSample.op, c(0.025, 0.975))

### if we look at the histogram of the 66 measurements
hist(ivTime, xlab='Speed of light measurements', main='', breaks=50)
## we can see the outlier measurements and the normal model is inappropriate for this problem

########### Model checking
## External validation --> use model to make predictions about future data and collect new data to compare with predictions
## OR use the existing data we have at the moment i.e. recycle

## Which quantities do we wish to compare (defining predictive quantities)
## POSTERIOR PREDICTIVE CHECKING
## observed data should look plausible under the posterior predictive distribution
## Draw simulated values from the joint posterior of Yrep and compare to Yobs and look for systematic differences

## Gelman [2008] P 144 -
## sample 66 values, 20 times, each time drawing a fresh draw of sd and mean from the joint posterior
mDraws = matrix(NA, nrow = 66, ncol=20)

for (i in 1:20){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  mDraws[,i] = rnorm(66, m, s)
}

p.old = par(mfrow=c(3, 3))
garbage = apply(mDraws, 2, function(x) hist(x, main='', xlab='', ylab=''))
hist(ivTime, xlab='Speed of light measurements', main='')

# One way to measure the discrepancy is to
# compare the smallest value in each hypothetical replicated dataset to Newcomb’s
# smallest observation, −44.


## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
T1_var = function(Yrep) return(var(Yrep))

########## simulate 200 test quantities
mDraws = matrix(NA, nrow = 66, ncol=200)

for (i in 1:200){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  mDraws[,i] = rnorm(66, m, s)
}

### get the test quantity from the test function
t1 = apply(mDraws, 2, T1_var)
par(p.old)
hist(t1, xlab='Test Quantity - Variance', main='', breaks=50)
abline(v = var(lData$vector), lwd=2)
getPValue(t1, var(lData$vector))