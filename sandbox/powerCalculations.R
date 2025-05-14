# Name: powerCalculations.R
# Auth: u.niazi@soton.ac.uk
# Date: 12/05/25
# Desc: Example of power calculations from a Bayesian Perspective


###### utility function from Doing Bayesian Data Analysis book
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}


dfData = read.csv('hierarchicalBinomial/mouseTumor.csv')
str(dfData)
dfData$group = factor(1:nrow(dfData))
str(dfData)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='hierarchicalBinomial/oneLevelBinomialHierarchicalModel.stan')


## setup the stan data
lStanData = list(Ntotal=nrow(dfData), Ngroups1=nlevels(dfData$group), 
                 NgroupsMap=as.numeric(dfData$group),
                 y=dfData$success,
                 N=dfData$trials)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, thin=20,
                      cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, digits=3)

## what is average theta from historical data
mThetas = do.call(cbind, extract(fit.stan, 'theta'))
dim(mThetas)
hist(mThetas)
ivThetas.historical = colMeans(mThetas)
hist(ivThetas.historical)
# mHDI.orig = apply(mThetas, 2, HDIofMCMC)
# iHDI.orig.diff = mHDI.orig[2, ] - mHDI.orig[1, ]
# iHDI.orig.width = mean(iHDI.orig.diff)

## extract hyperparameters
traceplot(fit.stan, c('alpha', 'beta'))

mPriors = extract(fit.stan, c('alpha', 'beta'))
mPriors = do.call(cbind, mPriors)
dim(mPriors)
head(mPriors)
plot(mPriors)

## simulate new experiment of given size or samples
## sample alpha and beta, then sample new thetas
simulateOne = function(alpha, beta, subjects, size){
  theta = rbeta(subjects, alpha, beta)
  return(rbinom(subjects, size, theta))
}

mNewData = sapply(sample(1:nrow(mPriors), size = 100),
                  function(x) simulateOne(mPriors[x,'alpha'], mPriors[x, 'beta'], subjects=20, size = 50)) 
dim(mNewData)

## fit model on new data sets (posterior predictive checks)
bGoal = NULL

for (i in 1:100){
  fExperiments = gl(20, 1)
  lStanData = list(Ntotal=nrow(mNewData), Ngroups1=nlevels(fExperiments), 
                   NgroupsMap=as.numeric(fExperiments),
                   y=mNewData[,i],
                   N=rep(50, times=nrow(mNewData)))
  fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=4,  
                        cores=4)#, control=list(adapt_delta=0.99, max_treedepth = 15))
  
  mThetas = do.call(cbind, extract(fit.stan.2, 'theta'))
  mHDI = apply(mThetas, 2, HDIofMCMC)
  # iHDI.diff = mHDI[2,] - mHDI[1,]
  # iHDI.width = mean(iHDI.diff)
  # bGoal = c(bGoal, any(mean(colMeans(mThetas)) > mean(ivThetas.historical)))
  bGoal = c(bGoal, all(mHDI[2,] < 0.3))
}

sum(bGoal)/length(bGoal)
# [1] 0.35 with 20 subjects and size 25
# [1] 0.93 with subjects=20 and size=50

par(mfrow=c(1,2))
hist(mPriors[,'alpha'], main='Hyperparameter Alpha', xlab='')
hist(mPriors[,'beta'], main='Hyperparameter Beta', xlab='')
plot(mPriors, pch=20, main='Correlation of hyperparameters')
hist(iNewData, main='Simulated Data Sets', xlab='')



