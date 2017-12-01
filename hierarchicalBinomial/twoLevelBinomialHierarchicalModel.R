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

fit.stan.2 = sampling(stanDso.2, data=lStanData, iter=1000, chains=2,
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
  return(rbeta(1000, iSuc+iAlpha, iFail+iBeta))
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
  iAlpha1 = (theta[grep('alpha1', names(theta))])
  iBeta1 = (theta[grep('beta1', names(theta))])
  
  ## data, binomial
  iSuc = data$success # resp
  iTotal = data$trials
  groupIndex = data$groupIndex # mapping variable
  iTheta = iTheta[groupIndex]
  
  # checks on parameters
  if (any(iAlpha1 < 0) | any(iBeta1 < 0)) return(-Inf)
  # write the priors and likelihood
  lp = sum(dbeta(iTheta, iAlpha1, iBeta1, log=T)) + 1
  lik = sum(dbinom(iSuc, iTotal, iTheta, log=T))
  val = lik + lp
  return(val)
}

lData = list(success=dfData.2$Hits, trials=dfData.2$AtBats, groupIndex=as.numeric(dfData.2$PriPos))
## set starting values
s1 = logit(lData$success / lData$trials)
names(s1) = rep('theta', times=length(s1))
s2 = c(alpha1=0.5, beta1=0.5)
start = c(s1, s2)

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

fit.2 = mylaplace(mylogpost, start, lData)

## lets take a sample from this
## parameters for the multivariate t density
tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 1000, lData)













#####################################################################
library(car) ## for logit function
library(LearnBayes)
library(numDeriv)
logit.inv = function(p) {exp(p)/(exp(p)+1) }


mylogpost = function(theta, data){
  ## parameters to track/estimate
  iTheta = logit.inv(theta[grep('theta', names(theta))])
  iOmega1 = logit.inv(theta[grep('omega1', names(theta))])
  iKappa1 = exp(theta[grep('kappa1', names(theta))])

  ## data, binomial
  iSuc = data$success # resp
  iTotal = data$trials
  groupIndex = data$groupIndex # mapping variable

  # checks on parameters
  if (any(iKappa1 < 2)) return(-Inf)
  # write the priors and likelihood
  lp1 = sum(dbeta(iOmega1, 0.5, 0.5, log=T)) + sum(dgamma(iKappa1, 0.5, 1e-4, log=T))
  ## calculate the parameters for the next level - i.e. first grouping level
  a = iOmega1*(iKappa1-2)+1
  b = (1-iOmega1)*(iKappa1-2)+1
  lp2 = sum(dbeta(iTheta, a[groupIndex], b[groupIndex], log=T))
  lik = sum(dbinom(iSuc, iTotal, iTheta, log=T))
  val = lik + lp1 + lp2
  return(val)
}

s1 = dfData$Hits / dfData$AtBats
names(s1) = rep('theta', times=length(s1))
s1 = logit(s1)

s2 = tapply(s1, dfData$PriPos, mean)
names(s2) = rep('omega1', times=length(s2))

s3 = rep(log(1000), times=nlevels(dfData$PriPos))
names(s3) = rep('kappa1', times=length(s3))

start = c(s1, s2, s3)
lData = list('success' = dfData$Hits,
             'trials' = dfData$AtBats,
             'groupIndex' = as.numeric(dfData$PriPos))

mylogpost(start, lData)
laplace(mylogpost, start, lData)



########### try a conjugate approach 
## get alpha and beta parameters for the beta distribution from the data
## see gelman p 583
getalphabeta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(c(alpha=al, beta=be))
}

getTheta = function(param){
  iAlpha = param[1];
  iBeta = param[2];
  iSuc = param[3];
  iFail = param[4];
  return(rbeta(1000, iSuc+iAlpha, iFail+iBeta))
}



## get the hyperparameters from all the data
ivThetas.data = dfData$Hits/dfData$AtBats + 0.5
## get hyperparameters using population data
l = getalphabeta(mean(ivThetas.data), var(ivThetas.data))

## get the parameters for the first level 
s = tapply(dfData$Hits, dfData$PriPos, sum)
t = tapply(dfData$AtBats, dfData$PriPos, sum)
m = cbind(l['alpha'], l['beta'], s, t-s)

mLevel1 = apply(m, 1, getTheta)

## use these to estimate the parameters for each level 2
## lets put this into a function to make it simpler
getBetaBinomialPosterior = function(suc, trials, grouping, prior.alpha, prior.beta){
  s = tapply(suc, grouping, sum)
  t = tapply(trials, grouping, sum)
  m = cbind(prior.alpha, prior.beta, s, t-s)
  return(apply(m, 1, getTheta))
}

mLevel2 = getBetaBinomialPosterior(dfData$Hits, dfData$AtBats, dfData$Player, l['alpha'], l['beta'])


# calculate priors for each sub group
mLevel2 = sapply(seq_along(levels(dfData$PriPos)), function(x) {
  lvl = levels(dfData$PriPos)[x]
  df = dfData[dfData$PriPos == lvl,]
  df = droplevels.data.frame(df)
  ivThetas.data = df$Hits/df$AtBats + 0.5
  ## get hyperparameters using population data
  l = getalphabeta(mean(ivThetas.data), var(ivThetas.data))
  return(getBetaBinomialPosterior(df$Hits, df$AtBats, df$Player, l['alpha'], l['beta']))
})

mLevel2 = do.call(cbind, mLevel2)

l = getalphabeta(mean(colMeans(mLevel1)), var(colMeans(mLevel1)))
mLevel2 = getBetaBinomialPosterior(dfData$Hits, dfData$AtBats, dfData$Player, l['alpha'], l['beta'])




## try a conjugate approach
# getTheta = function(iAlpha, iBeta, iSuc, iFail){
#   return(mean(rbeta(1000, iSuc+iAlpha, iFail+iBeta)))
# }



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


