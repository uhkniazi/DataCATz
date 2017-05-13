# Name: newcombSpeedOfLight.R
# Auth: uhkniazi
# Date: 11/05/2017
# Desc: posterior predictive checks for the model using Newcomb's speed of light data


library(LearnBayes)
## the example is taken from Simon Newcomb's experiment to measure the speed of light
## we follow Gelman [2008] Page 67 for the analysis.

set.seed(123) # for replication

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
# The procedure for carrying out a posterior predictive model check requires specifying a test
# quantity, T (y) or T (y, θ), and an appropriate predictive distribution for the replications
# y rep [Gelman 2008]
## variance
T1_var = function(Y) return(var(Y))
## is the model adequate except for the extreme tails
T1_symmetry = function(Y, th){
  Yq = quantile(Y, c(0.90, 0.10))
  return(abs(Yq[1]-th) - abs(Yq[2]-th))
} 

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

## mChecks
mChecks = matrix(NA, nrow=5, ncol=3)
rownames(mChecks) = c('Variance', 'Symmetry', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('Normal', 'NormalCont', 'T')
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = 66, ncol=200)
mThetas = matrix(NA, nrow=200, ncol=2)
colnames(mThetas) = c('mu', 'sd')

for (i in 1:200){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  mDraws[,i] = rnorm(66, m, s)
  mThetas[i,] = c(m, s)
}

mDraws.norm = mDraws
### get the test quantity from the test function
t1 = apply(mDraws, 2, T1_var)
par(p.old)
hist(t1, xlab='Test Quantity - Variance', main='', breaks=50)
abline(v = var(lData$vector), lwd=2)
mChecks['Variance', 1] = getPValue(t1, var(lData$vector))
# 0.48, the result from Figure 6.4 Gelman [2008]
# The sample variance does not make a good test statistic because it is a sufficient statistic of
# the model and thus, in the absence of an informative prior distribution, the posterior
# distribution will automatically be centered near the observed value. We are not at all
# surprised to find an estimated p-value close to 1/2 . [Gelman 2008]

## test for symmetry
t1 = sapply(seq_along(1:200), function(x) T1_symmetry(mDraws[,x], mThetas[x,'mu']))
t2 = sapply(seq_along(1:200), function(x) T1_symmetry(lData$vector, mThetas[x,'mu']))
plot(t2, t1, xlim=c(-12, 12), ylim=c(-12, 12), pch=20, xlab='Realized Value T(Yobs, Theta)',
     ylab='Test Value T(Yrep, Theta)', main='Symmetry Check')
abline(0,1)
mChecks['Symmetry', 1] = getPValue(t1, t2) # we should see somewhere around 0.2 on repeated simulations
# The estimated p-value is 0.26, implying that any observed asymmetry in the middle of the distribution can easily be
# explained by sampling variation. [Gelman 2008]

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lData$vector)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lData$vector)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lData$vector)
mChecks['Mean', 1] = getPValue(t1, t2)
# Major failures of the model, typically corresponding to extreme tail-area probabilities (less
# than 0.01 or more than 0.99), can be addressed by expanding the model appropriately. [Gelman 2008]
# The relevant goal is not to answer the question, ‘Do the data come from the assumed
# model?’ (to which the answer is almost always no), but to quantify the discrepancies between
# data and model, and assess whether they could have arisen by chance, under the model’s
# own assumptions. [Gelman 2008]

## define a second log posterior function for mixture with contaminated normal distribution
lp2 = function(theta, data){
  # we define the sigma on a log scale as optimizers work better
  # if scale parameters are well behaved
  s = exp(theta[2])
  m = theta[1]
  mix = 0.9
  cont = theta[3]
  d = data$vector # observed data vector
  log.lik = sum(log(dnorm(d, m, s) * mix + dnorm(d, m, s*cont) * (1-mix)))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=mean(ivTime), 'sigma'=log(sd(ivTime)), 'cont'=1)
lp2(start, lData)

op = optim(start, lp2, control = list(fnscale = -1), data=lData)
op$par
exp(op$par[2])

## try the laplace function from LearnBayes
fit2 = laplace(lp2, start, lData)
fit2
se2 = sqrt(diag(fit2$var))
se2
fit2$mode+1.96*se2
fit2$mode-1.96*se2

sigSample.op = rnorm(1000, fit2$mode['sigma'], se2['sigma'])
muSample.op = rnorm(1000, mean(lData$vector), exp(sigSample.op)/sqrt(length(lData$vector)))

########## simulate 200 test quantities
mDraws = matrix(NA, nrow = 66, ncol=200)
mThetas = matrix(NA, nrow=200, ncol=2)
colnames(mThetas) = c('mu', 'sd')

for (i in 1:200){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  ## this will take a sample from a contaminated normal distribution
  sam = function() {
    ind = rbinom(1, 1, 0.9)
    return(ind * rnorm(1, m, s) + (1-ind) * rnorm(1, m, s*fit2$mode['cont']))
  }
  mDraws[,i] = replicate(66, sam())
  mThetas[i,] = c(m, s)
}

mDraws.normCont = mDraws
## get the p-values for the test statistics
t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 2] = getPValue(t1, var(lData$vector))

## test for symmetry
t1 = sapply(seq_along(1:200), function(x) T1_symmetry(mDraws[,x], mThetas[x,'mu']))
t2 = sapply(seq_along(1:200), function(x) T1_symmetry(lData$vector, mThetas[x,'mu']))
plot(t2, t1, xlim=c(-12, 12), ylim=c(-12, 12), pch=20, xlab='Realized Value T(Yobs, Theta)',
     ylab='Test Value T(Yrep, Theta)', main='Symmetry Check')
abline(0,1)
mChecks['Symmetry', 2] = getPValue(t1, t2) 

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lData$vector)
mChecks['Min', 2] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lData$vector)
mChecks['Max', 2] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lData$vector)
mChecks['Mean', 2] = getPValue(t1, t2)

######################### try a third distribution, t with a low degrees of freedom
lp3 = function(theta, data){
  # function to use to use scale parameter
  ## see here https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
  dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
  ## likelihood function
  lf = function(dat, pred){
    return(log(dt_ls(dat, nu, pred, sigma)))
  }
  nu = exp(theta['nu']) ## normality parameter for t distribution
  sigma = exp(theta['sigma']) # scale parameter for t distribution
  m = theta[1]
  d = data$vector # observed data vector
  if (exp(nu) < 1) return(-Inf)
  log.lik = sum(lf(d, m))
  log.prior = 1
  log.post = log.lik + log.prior
  return(log.post)
}

# sanity check for function
# choose a starting value
start = c('mu'=mean(ivTime), 'sigma'=log(sd(ivTime)), 'nu'=log(2))
lp3(start, lData)

op = optim(start, lp3, control = list(fnscale = -1), data=lData)
op$par
exp(op$par[2:3])

## try the laplace function from LearnBayes
fit3 = laplace(lp3, start, lData)
fit3
se3 = sqrt(diag(fit3$var))

sigSample.op = rnorm(1000, fit3$mode['sigma'], se3['sigma'])
nuSample = exp(rnorm(1000, fit3$mode['nu'], se3['nu']))
# threshold the sample values to above or equal to 1
nuSample[nuSample < 1] = 1

## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
muSample.op = rnorm(1000, mean(lData$vector), exp(sigSample.op)/sqrt(length(lData$vector)))

########## simulate 200 test quantities
mDraws = matrix(NA, nrow = 66, ncol=200)
mThetas = matrix(NA, nrow=200, ncol=3)
colnames(mThetas) = c('mu', 'sd', 'nu')

for (i in 1:200){
  p = sample(1:1000, size = 1)
  s = exp(sigSample.op[p])
  m = muSample.op[p]
  n = nuSample[p]
  mDraws[,i] = rt_ls(66, n, m, s)
  mThetas[i,] = c(m, s, n)
}

mDraws.t = mDraws
## get the p-values for the test statistics
t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 3] = getPValue(t1, var(lData$vector))

## test for symmetry
t1 = sapply(seq_along(1:200), function(x) T1_symmetry(mDraws[,x], mThetas[x,'mu']))
t2 = sapply(seq_along(1:200), function(x) T1_symmetry(lData$vector, mThetas[x,'mu']))
plot(t2, t1, xlim=c(-12, 12), ylim=c(-12, 12), pch=20, xlab='Realized Value T(Yobs, Theta)',
     ylab='Test Value T(Yrep, Theta)', main='Symmetry Check')
abline(0,1)
mChecks['Symmetry', 3] = getPValue(t1, t2) 

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws, 2, T1_min)
t2 = T1_min(lData$vector)
mChecks['Min', 3] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws, 2, T1_max)
t2 = T1_max(lData$vector)
mChecks['Max', 3] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws, 2, T1_mean)
t2 = T1_mean(lData$vector)
mChecks['Mean', 3] = getPValue(t1, t2)

mChecks

par(mfrow=c(2,2))
hist(lData$vector, prob=T, breaks=30, xlim=c(-50, 60))
garbage = apply(mDraws.norm, 2, function(x) lines(density(x), lwd=0.5, lty=2, col='darkgrey'))

hist(lData$vector, prob=T, breaks=30, xlim=c(-50, 60))
garbage = apply(mDraws.normCont, 2, function(x) lines(density(x), lwd=0.5, lty=2, col='darkgrey'))

hist(lData$vector, prob=T, breaks=30, xlim=c(-50, 60))
garbage = apply(mDraws.t, 2, function(x) lines(density(x), lwd=0.5, lty=2, col='darkgrey'))
