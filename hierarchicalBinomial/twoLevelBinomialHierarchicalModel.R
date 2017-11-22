# Name: twoLevelBinomialHierarchicalModel.R
# Auth: uhkniazi
# Date: 09/11/2017
# Desc: model from [1] Kruschke, J. K. (2014). Doing Bayesian data analysis: Chapter 9


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

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
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
hist(mPositions[,'Pitcher'])
hist(mPositions[,'Pitcher'] - mPositions[,'Catcher'])
plot(mPositions[,'Pitcher'], mPositions[,'Catcher'], pch=20)
hist(mPositions[,'Catcher'])

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
hist(mPlayers[,'Mike Leake'] - mPlayers[,'Wandy Rodriguez'])

## another set of players
dfData[c(573, 428),]
## both have larger number of trials but different successes
## data influences the shrinkage and it is less than in previous case
hist(mPlayers[,'Andrew McCutchen'] - mPlayers[,'Brett Jackson'])

########### load second smaller data set with less parameters at lvl 2
## the data comes from a Bisulphite sequencing experiment 
## we look at the number of methylated cytosine bases at a certain position in the genome
## and the total number of coverage for that region.

dfData = read.csv('sampleData.csv', header=T)
str(dfData)

## our interest is to estimate the level 1 (group level) and sample level (level 2) parameters
dfData$interaction = factor(dfData$group1:dfData$group2)
xtabs(~ interaction + group3, data=dfData)

## add a one to each total to avoid zeros
dfData$total = dfData$total + 2
dfData$methylated = dfData$methylated + 1

## we try and do this 2 ways, model using STAN and using a log posterior function
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='twoLevelBinomialHierarchicalModel.stan')

## set up starting data
lStanData = list(Ntotal=nrow(dfData), NgroupsLvl1=nlevels(dfData$interaction), 
                 NgroupsLvl2Map=as.numeric(dfData$interaction),
                 y=dfData$methylated,
                 N=dfData$total)


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, digits=3)


########################################## 
## define a log posterior function and a different parameterization
#########################################

stanDso2 = rstan::stan_model(file='oneLevelBinomialHierarchicalModel.stan')

lStanData = list(Ntotal=nrow(dfData), NgroupsLvl1=nlevels(dfData$interaction), 
                 NgroupsLvl2Map=as.numeric(dfData$interaction),
                 y=dfData$methylated,
                 N=dfData$total)


fit.stan2 = sampling(stanDso2, data=lStanData, iter=1000, chains=2,
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan2, digits=3)

library(car) ## for logit function
library(LearnBayes)
library(numDeriv)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

## get alpha and beta for beta distribution from the data
## see gelman p 583
getalphabeta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(c(alpha=al, beta=be))
}

## try a conjugate approach
# getTheta = function(iAlpha, iBeta, iSuc, iFail){
#   return(mean(rbeta(1000, iSuc+iAlpha, iFail+iBeta)))
# }
getTheta = function(param){
  iAlpha = param[1];
  iBeta = param[2];
  iSuc = param[3];
  iFail = param[4];
  return(mean(rbeta(1000, iSuc+iAlpha, iFail+iBeta)))
}



## get prior for each group
s1 = dfData$methylated / dfData$total
m = tapply(s1, dfData$interaction, mean)
v = tapply(s1, dfData$interaction, var)
m = cbind(m, v)
l = sapply(1:nrow(m), function(x) getalphabeta(m[x,1], m[x,2]))
alpha = l['alpha',]
beta = l['beta',]
alpha = alpha[as.numeric(dfData$interaction)]
beta = beta[as.numeric(dfData$interaction)]
suc = dfData$methylated
fail = dfData$total - dfData$methylated
m = cbind(alpha, beta, suc, fail)
apply(m, 1, getTheta)



mylogpost = function(theta, data){
  ## parameters to track/estimate
  iTheta = logit.inv(theta[grep('theta', names(theta))])
  iAlpha1 = (theta[grep('alpha1', names(theta))])
  iBeta1 = (theta[grep('beta1', names(theta))])

  ## data, binomial
  iSuc = data$success # resp
  iTotal = data$trials
  groupIndex = data$groupIndex # mapping variable

  # checks on parameters
  if (any(iAlpha1 < 0) | any(iBeta1 < 0)) return(-Inf)
  # write the priors and likelihood
  lp = sum(dbeta(iTheta, iAlpha1[groupIndex], iBeta1[groupIndex], log=T)) + 1
  lik = sum(dbinom(iSuc, iTotal, iTheta, log=T))
  val = lik + lp
  return(val)
}

s1 = dfData$methylated / dfData$total
names(s1) = rep('theta', times=length(s1))
# means and variances for level 1
m = mean(s1); v = var(s1)
l = getalphabeta(m, v)
s1 = logit(s1)
s2 = rep(l$alpha, times=nlevels(dfData$interaction))
names(s2) = rep('alpha1', times=length(s2))
s3 = rep(l$beta, times=nlevels(dfData$interaction))
names(s3) = rep('beta1', times=length(s2))
start = c(s1, s2, s3)
lData = list('success' = dfData$methylated,
'trials' = dfData$total,
'groupIndex' = as.numeric(dfData$interaction))
start
mylogpost(start, lData)



# mylogpost = function(theta, data){
#   ## parameters to track/estimate
#   iTheta = logit.inv(theta[grep('theta', names(theta))])
#   iOmega1 = logit.inv(theta[grep('omega1', names(theta))])
#   iKappa1 = exp(theta[grep('kappa1', names(theta))])
#   
#   ## data, binomial
#   iSuc = data$success # resp
#   iTotal = data$trials
#   groupIndex = data$groupIndex # mapping variable
#   
#   # checks on parameters
#   if (any(iKappa1 < 2)) return(-Inf)
#   # write the priors and likelihood 
#   lp1 = sum(dbeta(iOmega1, 0.5, 0.5, log=T)) + sum(dgamma(iKappa1, 0.5, 1e-4, log=T))
#   ## calculate the parameters for the next level - i.e. first grouping level
#   a = iOmega1*(iKappa1-2)+1
#   b = (1-iOmega1)*(iKappa1-2)+1
#   lp2 = sum(dbeta(iTheta, a[groupIndex], b[groupIndex], log=T))
#   lik = sum(dbinom(iSuc, iTotal, iTheta, log=T))
#   val = lik + lp1 + lp2
#   return(val)
# }
# 
# s1 = dfData$methylated / dfData$total
# names(s1) = rep('theta', times=length(s1))
# s1 = logit(s1)
# 
# s2 = tapply(s1, dfData$interaction, mean)
# names(s2) = rep('omega1', times=length(s2))
# 
# s3 = rep(log(1000), times=nlevels(dfData$interaction))
# names(s3) = rep('kappa1', times=length(s3))
# 
# start = c(s1, s2, s3)
# lData = list('success' = dfData$methylated,
#              'trials' = dfData$total,
#              'groupIndex' = as.numeric(dfData$interaction))

mylogpost(start, lData)
library(optimx)
op = optimx(start, mylogpost, control = list(maximize=T, usenumDeriv=T, all.methods=T), data=lData)



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



fit = mylaplace(mylogpost, start, lData)


