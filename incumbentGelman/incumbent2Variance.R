# Name: incumbent2Variance.R
# Auth: uhkniazi
# Date: 26/06/2018
# Desc: fitting various regression models, following the example from Gelman 2013

## load the merged data
load('incumbentGelman/dataExternal/dfMerged.rdata')
dfData = dfMerged[dfMerged$year == 1988, ]

#####################################################################################################
############# Data sets produced, now fit model
#####################################################################################################
fit.1 = lm(response ~ treatment + previousProportion + incumbentParty, data=dfData)
summary(fit.1)

#################### use stan to generate MCMC sample
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='incumbentGelman/linearRegression2Variance.stan')

m = model.matrix(response ~ treatment + previousProportion + incumbentParty, data=dfData)

lStanData = list(Ntotal=nrow(dfData), NgroupMap=as.numeric(dfData$treatment), Ncol=ncol(m), X=m,
                 y=dfData$response)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'mu', 'sigmaPop'),
                    cores=2)
print(fit.stan, c('betas', 'sigmaPop'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('betas', 'sigmaPop'), ncol=1, inc_warmup=F)

## some sample diagnostic plots
library(coda)
library(lattice)
oCoda = As.mcmc.list(fit.stan, c('betas', 'sigmaPop'))
xyplot(oCoda[[1]])
autocorr.plot(oCoda[[1]])

m = extract(fit.stan, 'mu')
names(m)
dim(m$mu)
fitted = apply(m$mu, 2, mean)

m = extract(fit.stan, 'betas')
betas = colMeans(m$betas)
names(betas) = colnames(lStanData$X)
# compare with lm 
data.frame(coef(fit.1), betas)

s = cbind(extract(fit.stan)$betas, extract(fit.stan)$sigmaPop)
colnames(s) = c(colnames(lStanData$X), 'sigmaOpen', 'sigmaIncumbent')
pairs(s, pch=20)

## table 14.1
mTable = matrix(NA, nrow = 6, ncol = 5)
rownames(mTable) = c('incumbency', 'vote proportion 1986', 'incumbent party', 'intercept', 'residual sd Open', 'residual sd Incumbent')
mTable[1, ] = round(quantile(s[,2], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[2, ] = round(quantile(s[,3], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[3, ] = round(quantile(s[,4], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[4, ] = round(quantile(s[,1], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[5, ] = round(quantile((s[,5]), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[6, ] = round(quantile((s[,6]), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)

#########################################################
########## some model checks on 1988
#########################################################
mStan = s
# get the fitted value
iBetas = colMeans(mStan[,1:4])
fitted = lStanData$X %*% iBetas

plot(dfData$response, fitted, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5)
# get residuals that is response minus fitted values
iResid = (dfData$response - fitted)
plot(dfData$response, iResid, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5)

plot(fitted, iResid, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='SIR')
lines(lowess(fitted, iResid), col=2, lwd=2)

plot(predict(fit.1), resid(fit.1), pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='lm')
lines(lowess(predict(fit.1), resid(fit.1)), col=2, lwd=2)

## calculate standardized residuals
library(MASS)
plot(fitted, stdres(fit.1), pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='standardized residuals - lm')
lines(lowess(fitted, stdres(fit.1)), col=2, lwd=2)
## see equation 14.7 in Gelman 2013
s = colMeans(mStan[,5:6])
s = s[as.numeric(dfData$treatment)]
plot(fitted, iResid/s, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='standardized residuals - SIR')
lines(lowess(fitted, iResid/s), col=2, lwd=2)

plot(stdres(fit.1), iResid/s, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5)

### generate some posterior predictive data
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(betas, sigma, mModMatrix, groupMap){
  f = mModMatrix %*% betas
  yrep = rnorm(length(f), f, sigma[groupMap])
  return(yrep)
}

runRegression = function(yrep, lStanData){
  lStanData$y = yrep
  f.s = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'mu', 'sigmaPop'),
                      cores=2)
  return(f.s)
}

# returns standardized residuals
getResiduals = function(fit.object, mModMatrix, groupMap, yrep){
  b = colMeans(extract(fit.object)$betas)
  f = mModMatrix %*% b
  r = (yrep - f)/extract(fit.object)$sigmaPop[groupMap]
  return(r)
}


## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = nrow(dfData), ncol=1000)
mDraws.fitted = matrix(NA, nrow = nrow(dfData), ncol=1000)
mDraws.res = matrix(NA, nrow = nrow(dfData), ncol=1000)

for (i in 1:1000){
  p = sample(1:nrow(mStan), 1)
  sigma = mStan[p,5:6]
  betas = mStan[p, 1:4]
  mDraws.sim[,i] = simulateOne(betas, sigma, lStanData$X, as.numeric(dfData$treatment))
  f.s = runRegression(mDraws.sim[,i], lStanData)
  mDraws.fitted[,i] = apply(extract(f.s)$mu, 2, mean)
  mDraws.res[,i] = getResiduals(f.s, lStanData$X, groupMap = as.numeric(dfData$treatment),
                                mDraws.sim[,i])
}

## add a section here to get residuals not standardized

### visual checks
ivResp = dfData$response
yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)#, ylim=c(0, 1))
temp = apply(mDraws.sim, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
lines(yresp)

### plot the residuals
plot(dfData$response, iResid, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='SIR')
lines(lowess(dfData$response, iResid))
temp = sapply(1:1000, function(x) lines(lowess(dfData$response, mDraws.res[,x]), lwd=0.5, col=2))

plot(density(iResid))
g = apply(mDraws.res, 2, function(x) lines(density(x), lwd=0.5, col=2))
lines(density(iResid))

###############################################################
######### test quantities for model fits
###############################################################
## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}

## test for variance
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

## proportion of absolute value residuals that exceed limit iResLimit
T1_proportion = function(Y) {
  return(sum(abs(Y) > iResLimit)/length(Y))
}

T1_proportion.open = function(Y) {
  return(sum(abs(Y)[dfData$treatment == 'open'] > iResLimit)/length(Y[dfData$treatment == 'open']))
}

T1_proportion.inc = function(Y) {
  return(sum(abs(Y)[dfData$treatment == 'incumbent'] > iResLimit)/length(Y[dfData$treatment == 'incumbent']))
}


### search for outliers by checking the simulated residuals
#s = exp(mean(fit.2$sir[,'sigmaPop']))
## 3 standard deviations - this is the sd for the normal distribution used in the likelihood
iResLimit = 2; #round(s*3, 1)

## generate the test statistics and the observed value of that statistic and p-values
# observed numbers of outliers in residuals
iProp.observed = sum(abs(iResid) > iResLimit)/length(iResid)
iProp.observed.open = sum(abs(iResid)[dfData$treatment == 'open'] > iResLimit)/length(iResid[dfData$treatment == 'open'])
iProp.observed.inc = sum(abs(iResid)[dfData$treatment == 'incumbent'] > iResLimit)/length(iResid[dfData$treatment == 'incumbent'])
iVar.observed = var(dfData$response)
iMin.observed = min(dfData$response)
iMax.observed = max(dfData$response)
iMean.observed = mean(dfData$response)

# simulated quantites using posterior predictive data
iProp.sim = apply(mDraws.res, 2, T1_proportion)
iProp.sim.open = apply(mDraws.res, 2, T1_proportion.open)
iProp.sim.inc = apply(mDraws.res, 2, T1_proportion.inc)
iVar.sim = apply(mDraws.sim, 2, T1_var)
iMin.sim = apply(mDraws.sim, 2, T1_min)
iMax.sim = apply(mDraws.sim, 2, T1_max)
iMean.sim = apply(mDraws.sim, 2, T1_mean)

mChecks = matrix(NA, nrow=7, ncol=1)
rownames(mChecks) = c('Variance', 'Residuals', 'Residual.open', 'Residual.inc', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('Model')

mChecks['Variance',] = getPValue(iVar.sim, iVar.observed)
mChecks['Residuals', ] = getPValue(iProp.sim, iProp.observed)
mChecks['Max',] = getPValue(iMax.sim, iMax.observed)
mChecks['Min',] = getPValue(iMin.sim, iMin.observed)
mChecks['Mean',] = getPValue(iMean.sim, iMean.observed)
mChecks['Residual.open', ] = getPValue(iProp.sim.open, iProp.observed.open)
mChecks['Residual.inc', ] = getPValue(iProp.sim.inc, iProp.observed.inc)
mChecks


######################################################################################
### repeat this on all the years 
#####################################################################################

ivYears = unique(dfMerged$year)
lYears = vector('list', length = length(ivYears))

for (yr in seq_along(ivYears)){
  dfData = dfMerged[dfMerged$year == ivYears[yr], ]
  
  ## fit models for each year
  lData = list(resp=dfData$response, mModMatrix=model.matrix(response ~ treatment + previousProportion + incumbentParty, data=dfData))
  start = c('sigmaPop'=log(2), 'betas'=rep(0, times=ncol(lData$mModMatrix)))

  fit.2 = mylaplace(mylogpost, start, lData)
  ### lets take a sample from this using SIR
  ## parameters for the multivariate t density
  tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
  ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
  s = sir(mylogpost, tpar, 1000, lData)
  colnames(s)[-1] = colnames(lData$mModMatrix)
  lYears[[yr]] = s
}

names(lYears) = ivYears

## utility function for plotting
getms = function(f){
  m = mean(f)
  se = sd(f)
  m.up = m+1.96*se
  m.down = m-1.96*se
  ret= c(m, m.up, m.down)
  names(ret) = c('m', 'm.up', 'm.down')
  return(ret)
}

mTreatment = sapply(lYears, function(x) x[,3])

df = apply(mTreatment, 2, getms)
x = 1:ncol(mTreatment)

## reproduce figure 14.2 from gelman 2013 
plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='Year', 
     ylab='estimated incumbency advantage', main='', xaxt='n')
axis(1, at = x, labels = colnames(mTreatment), las=2)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}



