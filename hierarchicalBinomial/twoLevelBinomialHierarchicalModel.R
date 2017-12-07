# Name: twoLevelBinomialHierarchicalModel.R
# Auth: uhkniazi
# Date: 09/11/2017
# Desc: models from [1] Kruschke, J. K. (2014). Doing Bayesian data analysis: Chapter 9
#       and [2] Gelman (2013). Bayesian Data Analysis: Chapter 5


library(LearnBayes)

## start by analysing the data from the book examples
## the data is from the book
dfData = read.csv('BattingAverage.csv')
## sort the data by the two levels of interest
##dfData = dfData[order(dfData$Player, dfData$PriPos),]

## this is a 2 level model, see figure 9.13 page 252 [1]
# level 0 = highest level
# level 1 = population level with multiple groups
# level 2 = lowest level of hierarchy
#               Level 0
#                 |
#       -------------------               LEVEL 1
#       |         |       |
#       L1.1      L1.2    L1.3
# ---------    --------      --------     LEVEL 2
# |       |    |      |      |      |
# L.1.1  L.1.2 L.2.1  L.2.2  L3.1   L3.2
# |
# data
## one Level 2 parameter for each data point
str(dfData)

## we try and do this in different ways, first model using STAN
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='twoLevelBinomialHierarchicalModel.stan')

## set up starting data
lStanData = list(Ntotal=nrow(dfData), NgroupsLvl1=nlevels(dfData$PriPos), 
                 NgroupsLvl2Map=as.numeric(dfData$PriPos),
                 y=dfData$Hits,
                 N=dfData$AtBats)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2,
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('omega1', 'kappa1', 'kappa0', 'omega0'), digits=3)

traceplot(fit.stan, 'omega0')
traceplot(fit.stan, 'kappa0')
traceplot(fit.stan, 'omega1')
traceplot(fit.stan, 'kappa1')

mPositions = extract(fit.stan)$omega1
dim(mPositions)
colnames(mPositions) = levels(dfData$PriPos)

## figure 9.14 from [1]
par(mfrow=c(2,2))
hist(mPositions[,'Pitcher'], main='Pitcher', xlab='', ylab='')
hist(mPositions[,'Pitcher'] - mPositions[,'Catcher'], main='Pitcher - Catcher', xlab='', ylab='')
plot(mPositions[,'Pitcher'], mPositions[,'Catcher'], pch=20, main='', xlab='Pitcher', ylab='Catcher')
hist(mPositions[,'Catcher'], main='Catcher', xlab='', ylab='')

## second panel of the figure
hist(mPositions[,'Catcher'])
hist(mPositions[,'Catcher'] - mPositions[,'1st Base'])
plot(mPositions[,'Catcher'], mPositions[,'1st Base'], pch=20)
hist(mPositions[,'1st Base'])

mPlayers = extract(fit.stan)$theta
colnames(mPlayers) = as.character(dfData$Player)

## figure 9.15 [1]
hist(mPlayers[,'Kyle Blanks'])
hist(mPlayers[,'Kyle Blanks'] - mPlayers[,'Bruce Chen'])
plot(mPlayers[,'Kyle Blanks'], mPlayers[,'Bruce Chen'], pch=20)
hist(mPlayers[,'Bruce Chen'])

## example of shrinkage
## difference between two players 
dfData[c(494, 754),]
## both have same number of trials but different number of successes
## difference includes 0 - due to shrinkage towards the general trend for pitchers
hist(mPlayers[,'Mike Leake'] - mPlayers[,'Wandy Rodriguez'], xlab='', ylab='', main='Same Trials, Different Success',
     sub='Both Same Position')

## another set of players
dfData[c(573, 428),]
## both have larger number of trials but different successes
## data influences the shrinkage and it is less than in previous case
hist(mPlayers[,'Andrew McCutchen'] - mPlayers[,'Brett Jackson'], xlab='', ylab='', main='Large Trials, Different Success',
     sub='Both Same Position')

##############################################################################
### merge the lvl 2 data i.e. players, as we are only interested in 
### positions, trying a different parameterization of the model
d1 = tapply(dfData$Hits, dfData$PriPos, sum)
d2 = tapply(dfData$AtBats, dfData$PriPos, sum)
dfData.2 = data.frame(PriPos = factor(levels(dfData$PriPos)), Hits=d1, AtBats=d2)

stanDso.2 = rstan::stan_model(file='oneLevelBinomialHierarchicalModel.stan')

## set up starting data
lStanData = list(Ntotal=nrow(dfData.2), Ngroups1=nlevels(dfData.2$PriPos), 
                 NgroupsMap=as.numeric(dfData.2$PriPos),
                 y=dfData.2$Hits,
                 N=dfData.2$AtBats)

fit.stan.2 = sampling(stanDso.2, data=lStanData, iter=2000, chains=2,
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan.2, digits=3)

mPositions.2 = extract(fit.stan.2)$theta
dim(mPositions.2)
colnames(mPositions.2) = levels(dfData.2$PriPos)

## figure 9.14 from [1]
par(mfrow=c(2,2))
hist(mPositions.2[,'Pitcher'])
hist(mPositions.2[,'Pitcher'] - mPositions.2[,'Catcher'])
plot(mPositions.2[,'Pitcher'], mPositions.2[,'Catcher'], pch=20)
hist(mPositions.2[,'Catcher'])

## compare the 2 models
par(mfrow=c(1,2))
plot(colMeans(mPositions), colMeans(mPositions.2), pch=20, xlim=c(0.23, 0.27), ylim=c(0.23, 0.27))

########### try a conjugate approach 
## get alpha and beta parameters for the beta distribution from the data
## see gelman [2] p 583
getalphabeta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(c(alpha=al, beta=be))
}

## take a sample from beta posterior 
getFittedTheta = function(param){
  iAlpha = param[1];
  iBeta = param[2];
  iSuc = param[3];
  iFail = param[4];
  return(rbeta(2000, iSuc+iAlpha, iFail+iBeta))
}


## get the hyperparameters using all the data
ivThetas.data = dfData.2$Hits/dfData.2$AtBats
## get hyperparameters using population data
l = getalphabeta(mean(ivThetas.data), var(ivThetas.data))

## use the data to calculate success and failures
suc = dfData.2$Hits
fail = dfData.2$AtBats - dfData.2$Hits
# put data in matrix form to cycle the function over each row
m = cbind(l['alpha'], l['beta'], suc, fail)

mPositions.conj = apply(m, 1, getFittedTheta)

plot(colMeans(mPositions), colMeans(mPositions.conj), pch=20, xlim=c(0.23, 0.27), ylim=c(0.23, 0.27))
points(colMeans(mPositions), colMeans(mPositions.2), pch=20, col=2)

data.frame(m.2lvl = round(colMeans(mPositions), 2),
           m.1lvl = round(colMeans(mPositions.2), 2),
           m.conj = round(colMeans(mPositions.conj), 2))


##################################################################
### lets try an optimisation based approach
library(car) ## for logit function
library(LearnBayes)
library(numDeriv)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

### define log posterior
mylogpost = function(theta, data){
  ## parameters to track/estimate
  iTheta = logit.inv(theta[grep('theta', names(theta))])
  iAlDivBe = logit.inv(theta[grep('logitAlDivBe', names(theta))])
  iAlPlusBe = exp(theta[grep('logAlPlusBe', names(theta))])
  
  ## data, binomial
  iSuc = data$success # resp
  iTotal = data$trials
  groupIndex = data$groupIndex # mapping variable
  iTheta = iTheta[groupIndex]
  
  ### we solve two equations together to get the values for alpha and beta
  # log(alpha/beta) = iLogAlDivBe
  # log(alpha+beta) = iLogAlPlusBe
  # solve them together
  A = iAlDivBe
  B = iAlPlusBe
  iAlpha = B/(1+1/A)
  iBeta = B - iAlpha
  
  # checks on parameters
  if (iAlpha < 0 | iBeta < 0) return(-Inf)
  # write the priors and likelihood
  lp = sum(dbeta(iTheta, iAlpha, iBeta, log=T)) + dbeta(A, 0.5, 0.5, log=T) + dgamma(B, 0.5, 1e-4, log=T)
  lik = sum(dbinom(iSuc, iTotal, iTheta, log=T))
  val = lik + lp
  return(val)
}

lData = list(success=dfData.2$Hits, trials=dfData.2$AtBats, groupIndex=as.numeric(dfData.2$PriPos))
## set starting values
s1 = lData$success / lData$trials
names(s1) = rep('theta', times=length(s1))
l = getalphabeta(mean(s1), var(s1))

## logit(alpha/(alpha+beta)) = log(alpha/beta)
s2 = c(logit(l['alpha']/sum(l)), log(sum(l)))
names(s2) = c('logitAlDivBe', 'logAlPlusBe')
start = c(logit(s1), s2)

## test function
mylogpost(start, lData)

fit.1 = laplace(mylogpost, start, lData)

# define a modification of the laplace function from learnbayes
library(numDeriv)
library(optimx)

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=20000), method='Nelder-Mead', data=data)
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

fit.2 = mylaplace(mylogpost, start, lData)

logit.inv(fit.2$mode[1:10])
exp(fit.2$mode[11])

## lets take a sample from this
## parameters for the multivariate t density
tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 5000, lData)
mPositions.opt = logit.inv(s[,1:9])
colnames(mPositions.opt) = levels(dfData.2$PriPos)

data.frame(m.2lvl = round(colMeans(mPositions), 2),
           m.1lvl = round(colMeans(mPositions.2), 2),
           m.conj = round(colMeans(mPositions.conj), 2),
           m.opt = round(colMeans(mPositions.opt), 2))

data.frame(m.2lvl = round(apply(mPositions, 2, sd), 5),
           m.1lvl = round(apply(mPositions.2, 2, sd), 5),
           m.conj = round(apply(mPositions.conj, 2, sd), 5),
           m.opt = round(apply(mPositions.opt, 2, sd), 5))



### lets compare the results from the 4 different analyses
p1 = density(mPositions[,'Pitcher'])
p1$y = p1$y/max(p1$y)
par(mfrow=c(1,1))
plot(p1, main='Pitcher', lwd=2) 
p1 = density(mPositions.2[,'Pitcher'])
p1$y = p1$y/max(p1$y)
lines(p1, col=2, lwd=2)
p1 = density(mPositions.conj[,'Pitcher'])
p1$y = p1$y/max(p1$y)
lines(p1, col=3, lwd=2)
p1 = density(mPositions.opt[,'Pitcher'])
p1$y = p1$y/max(p1$y)
lines(p1, col=4, lwd=2)

legend('topleft', legend = c('2 Level', '1 Level', 'Conjugate', 'SIR'), fill=c(1, 2, 3, 4))


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

df = apply(mPositions, 2, getms)
x = jitter(dfData.2$Hits/dfData.2$AtBats, factor = 10)

## reproduce figure from book [2] 
par(mfrow=c(1,2))
plot(x, df['m',], pch=20, xlab=expression(paste('Observed Rate  ',frac(Hits[i],AtBats[i]))), xlim=c(0.243, 0.268), cex.axis=0.5,
     ylim=c(0.23, 0.28),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='2 level Hierarchical Model')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)

df = apply(mPositions.2, 2, getms)

plot(x, df['m',], pch=20, xlab=expression(paste('Observed Rate  ',frac(Hits[i],AtBats[i]))), xlim=c(0.243, 0.268), cex.axis=0.5,
     ylim=c(0.23, 0.28),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='1 level Hierarchical Model')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)

df = apply(mPositions.conj, 2, getms)

plot(x, df['m',], pch=20, xlab=expression(paste('Observed Rate  ',frac(Hits[i],AtBats[i]))), xlim=c(0.243, 0.268), cex.axis=0.5,
     ylim=c(0.23, 0.28),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='Conjugate EBayes Model')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)

df = apply(mPositions.opt, 2, getms)

plot(x, df['m',], pch=20, xlab=expression(paste('Observed Rate  ',frac(Hits[i],AtBats[i]))), xlim=c(0.243, 0.268), cex.axis=0.5,
     ylim=c(0.23, 0.28),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='SIR Model')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)



##########################################################################################
##########################################################################################
## application to the mouse data set from [2]
rm(dfData)
rm(dfData.2)

dfData = read.csv('mouseTumor.csv')
str(dfData)
dfData = rbind(dfData, c(4, 14))
dfData$group = factor(1:nrow(dfData))
str(dfData)

## setup the stan data
lStanData = list(Ntotal=nrow(dfData), Ngroups1=nlevels(dfData$group), 
                 NgroupsMap=as.numeric(dfData$group),
                 y=dfData$success,
                 N=dfData$trials)

fit.stan.2 = sampling(stanDso.2, data=lStanData, iter=2000, chains=2,
                      cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan.2, digits=3)

mMutations = extract(fit.stan.2)$theta
colnames(mMutations) = levels(dfData$group)

### try a conjugate, ebayes approach
nrow(dfData)
ivThetas.data = dfData$success[1:77]/dfData$trials[1:77]
## get hyperparameters using population data
l = getalphabeta(mean(ivThetas.data), var(ivThetas.data))

## use the data to calculate success and failures
suc = dfData$success[78]
fail = dfData$trials[78] - dfData$success[78]

ivMutations.conj = getFittedTheta(c(l['alpha'], l['beta'], suc, fail))

############################ perform on all data
### try a conjugate, ebayes approach
nrow(dfData)
ivThetas.data = dfData$success/dfData$trials
## get hyperparameters using population data
l = getalphabeta(mean(ivThetas.data), var(ivThetas.data))

## use the data to calculate success and failures
suc = dfData$success
fail = dfData$trials - dfData$success
trials = dfData$trials
# put data in matrix form to cycle the function over each row
m = cbind(l['alpha'], l['beta'], suc, fail)

mMutations.conj = apply(m, 1, getFittedTheta)

data.frame(m.1lvl = round(colMeans(mMutations), 2),
           m.conj = round(colMeans(mMutations.conj), 2))

data.frame(m.1lvl = round(apply(mMutations, 2, sd), 3),
           m.conj = round(apply(mMutations.conj, 2, sd), 3))

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

df = apply(mMutations, 2, getms)
x = jitter(suc/trials, factor = 10)

## reproduce figure from book [2] 
plot(x, df['m',], ylim=c(0, 0.42), xlim=c(0, 0.42), pch=20, xlab=expression(paste('Observed Rate  ',frac(Y[i],N[i]))),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='1 level Hierarchical Model')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)

df = apply(mMutations.conj, 2, getms)

plot(x, df['m',], ylim=c(0, 0.42), xlim=c(0, 0.42), pch=20, xlab=expression(paste('Observed Rate  ',frac(Y[i],N[i]))),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='EBayes Conjugate Model')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)

### plot both together
df = apply(mMutations, 2, getms)
x = jitter(suc/trials, factor = 10)

plot(x, df['m',], ylim=c(0, 0.42), xlim=c(0, 0.42), pch=20, xlab=expression(paste('Observed Rate  ',frac(Y[i],N[i]))),
     ylab=expression(paste('Model Estimated  ', theta[i])), main='Fits from Models')
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(0, 1)

df = apply(mMutations.conj, 2, getms)

points(x, df['m',], col=2, pch=20)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5, col=2)
}

legend('topleft', legend = c('Hierarchical', 'EBayes'), fill=c('black', 'red'))




