# Name: incumbent.R
# Auth: uhkniazi
# Date: 26/06/2018
# Desc: fitting various regression models, following the example from Gelman 2013

######################### fit a simple regression model first after loading the 1988 dataset
dfData = dfData = read.table('incumbentGelman/dataExternal/1988.asc', header=F)
colnames(dfData) = c('state', 'district', 'incumbency', 'democratVotes', 'republicanVotes')
str(dfData)

## previous election
dfData.prev = read.table('incumbentGelman/dataExternal/1986.asc', header=F)
colnames(dfData.prev) = c('state', 'district', 'incumbency', 'democratVotes', 'republicanVotes')


getWinner = function(d, r){
  w = -1
  if ((d == -9 | r == -9) & (d > 0 | r > 0)) {
    w = ifelse(d > r, 'd', 'r')
  } else if ((d > 0 & r > 0)) {
    w = ifelse(d > r, 'd', 'r')
  }
  return(w)
}

getTreatment = function(inc) {
  ifelse(inc == 0, 'open', 'incumbent')
}

# getResponse = function(inc, d, r, winner){
#   y = -1
#   if (inc == 0) {
#     # open seat
#     if (winner == 'd') y = d/(d+r) else y = r/(d+r)
#   } else if (inc == 1){
#     # not open seat
#     y = d / (d+r)
#   } else if (inc == -1) {
#     y = r / (d + r)
#   }
#   return(y)
# }

getResponse = function(inc, d, r){
  y = -1
  if (inc == 'r') {
    y = r / (d + r)
  } else if (inc == 'd'){
    y = d / (d+r)
  }
  return(y)
}

# check if the samples i.e. units of analysis are concordant in the 2 datasets
identical(dfData$state, dfData.prev$state)
identical(dfData$district, dfData.prev$district)
## get the winner for previous year
w = sapply(1:nrow(dfData.prev), function(x) getWinner(dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x]))
table(w)
## drop -1s as these are not concordant results
i = which(w == -1)
dfData = dfData[-i,]
dfData.prev = dfData.prev[-i,]
identical(dfData$state, dfData.prev$state)
identical(dfData$district, dfData.prev$district)

## perform same check for the current year
w = sapply(1:nrow(dfData), function(x) getWinner(dfData$democratVotes[x], dfData$republicanVotes[x]))
table(w)
## drop -1s as these are not concordant results
i = which(w == -1)
dfData = dfData[-i,]
dfData.prev = dfData.prev[-i,]
identical(dfData$state, dfData.prev$state)
identical(dfData$district, dfData.prev$district)

## get the winner of previous year as that is the current incumbent party
w = sapply(1:nrow(dfData.prev), function(x) getWinner(dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x]))
table(w)
dfData$incumbentParty = factor(w)

## get the treatment variable i.e. if open seat or incumbent candidate running
t = getTreatment(dfData$incumbency)
table(t)
dfData$treatment = factor(t, levels = c('open', 'incumbent'))

## get previous election results for current incumbent party, i.e. winner of previous elections
w = sapply(1:nrow(dfData.prev), function(x) getWinner(dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x]))
table(w)
p = sapply(1:nrow(dfData.prev), function(x){
  getResponse(w[x], dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x])
}) 
summary(p)
dfData$previousProportion = p

## get the current year proportion for the incumbent party
p = sapply(1:nrow(dfData), function(x){
  getResponse(dfData$incumbentParty[x], dfData$democratVotes[x], dfData$republicanVotes[x])
}) 
summary(p)
dfData$response = p

## incumbent party variable 
i = ifelse(as.character(dfData$incumbentParty) == 'd', 1, -1)
table(i, dfData$incumbentParty)
dfData$incumbentParty = i
# get the democratic vote for current and previous years
prev = dfData.prev$democratVotes / (dfData.prev$democratVotes + dfData.prev$republicanVotes)
cur = dfData$democratVotes / (dfData$democratVotes + dfData$republicanVotes)
summary(prev)
summary(cur)


## produce figure 14.1
plot(prev, cur, pch=c(2,20)[as.numeric(dfData$treatment)], xlab='1986', ylab='1988',
     xlim=c(0,1), ylim=c(0,1), main='Democratic Vote in 1988 vs 1986')


## load the merged data
load('incumbentGelman/dataExternal/dfMerged.rds')
dfData = dfMerged[dfMerged$year == 1988, ]

#####################################################################################################
############# Data sets produced, now fit model
#####################################################################################################
fit.1 = lm(response ~ treatment + previousProportion + incumbentParty, data=dfData)
summary(fit.1)

## we can try the models in 2 ways, using a log posterior function in R and later in Stan
library(LearnBayes)

library(numDeriv)
#library(optimx) optional optimiser to use
# op = optimx(start, mylogpost, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), method='Nelder-Mead', data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}
## write the log posterior function

mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta[grep('betas', names(theta))] # vector of betas i.e. regression coefficients for population
  sigmaPop = exp(theta['sigmaPop']) # population level sd
  ## data
  resp = data$resp # response variable
  mModMatrix = data$mModMatrix # design matrix 
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # write the priors and likelihood 
  lp = dcauchy(betas[1], 0, 5, log=T) + sum(dcauchy(betas[-1], 0, 2, log=T))
  lik = sum(dnorm(resp,mean = iFitted, sd=sigmaPop, log=T))
  val = lik + lp
  return(val)
}

lData = list(resp=dfData$response, mModMatrix=model.matrix(response ~ treatment + previousProportion + incumbentParty, data=dfData))
start = c('sigmaPop'=log(2), 'betas'=rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.2 = mylaplace(mylogpost, start, lData)
fit.2
# compare with lm 
data.frame(coef(fit.1), fit.2$mode[-1])
se = sqrt(diag(fit.2$var))


### lets take a sample from this using SIR
## parameters for the multivariate t density
tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 5000, lData)
colnames(s)[-1] = colnames(lData$mModMatrix)
apply(s, 2, mean)
apply(s, 2, sd)
exp(mean(s[,'sigmaPop']))
pairs(s, pch=20)
fit.2$sir = s

## table 14.1
mTable = matrix(NA, nrow = 5, ncol = 5)
rownames(mTable) = c('incumbency', 'vote proportion 1986', 'incumbent party', 'intercept', 'residual sd')
mTable[1, ] = round(quantile(s[,3], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[2, ] = round(quantile(s[,4], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[3, ] = round(quantile(s[,5], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[4, ] = round(quantile(s[,2], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)
mTable[5, ] = round(quantile(exp(s[,1]), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)

#########################################################
########## some model checks on 1988
#########################################################

# get the fitted value
iBetas = colMeans(fit.2$sir)[-1]
fitted = lData$mModMatrix %*% iBetas

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
s = exp(mean(fit.2$sir[,'sigmaPop']))
plot(fitted, iResid/s, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='standardized residuals - SIR')
lines(lowess(fitted, iResid/s), col=2, lwd=2)

### generate some posterior predictive data
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(betas, sigma, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rnorm(length(f), f, sigma)
  return(yrep)
}

runRegression = function(yrep, df){
  df$response = yrep
  return(lm(response ~ treatment + previousProportion + incumbentParty, data=df))
}

getResiduals = function(fit.object){
  return(resid(fit.object))
}


## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = nrow(dfData), ncol=1000)
mDraws.fitted = matrix(NA, nrow = nrow(dfData), ncol=1000)
mDraws.res = matrix(NA, nrow = nrow(dfData), ncol=1000)

for (i in 1:1000){
  p = sample(1:nrow(fit.2$sir), 1)
  sigma = exp(fit.2$sir[p,'sigmaPop'])
  betas = fit.2$sir[p,-1]
  mDraws.sim[,i] = simulateOne(betas, sigma, lData$mModMatrix)
  mDraws.fitted[,i] = fitted(runRegression(mDraws.sim[,i], dfData))
  mDraws.res[,i] = getResiduals(runRegression(mDraws.sim[,i], dfData))
}

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
s = exp(mean(fit.2$sir[,'sigmaPop']))
## 3 standard deviations - this is the sd for the normal distribution used in the likelihood
iResLimit = round(s*3, 1)

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



