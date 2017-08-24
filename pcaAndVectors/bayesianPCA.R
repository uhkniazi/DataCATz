# Name: bayesianPCA.R
# Auth: uhkniazi
# Date: 18/08/2017
# Desc: perform a usual and bayesian version of the pca


## sample data
# mData = cbind(c(2.5,0.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1),
#               c(2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9))
# colnames(mData) = c('X', 'Y')

data("faithful")
mData = as.matrix(faithful)
 
## structure of the data
str(mData)

plot(mData, pch=20)
## now standardize the data with 0 mean and variance 1
plot(scale(mData), pch=20)

## now perform PCA and data whitening
ivMeans = colMeans(mData)
# centered data
mData.s = sweep(mData, 2, ivMeans, '-')

# covariance matrix for the data
mCov = cov(mData)
lEigens = eigen(mCov)
colnames(lEigens$vectors) = c('Vector1', 'Vector2')

## inputs are the variables i.e. the data, one column at a time
t(mData.s)
## operations/transformation matrix is the matrix of eigen vectors
t(lEigens$vectors)

## get the transformed points after running them through the matrix
mData.rotated = t(lEigens$vectors) %*% t(mData.s)
rownames(mData.rotated) = c('Comp1', 'Comp2')
mData.rotated

plot(t(mData.rotated), main='Whitened')

## get the original data back
## rowDataMatrix = (inverse(rowEigenVectors) * rotated Data) + original Means
solve(t(lEigens$vectors)) ## this equals the transpose of the rowEigenVectors matrix
mData.original.s = (solve(t(lEigens$vectors)) %*% mData.rotated)
mData.original.s

## add the mean to un-center the data
mData.original = sweep(mData.original.s, 1, ivMeans, '+')
t(mData)

# perform PCA using the prcomp
# the vectors are in columns
pr.out = prcomp(mData, scale=F, center = T)
pr.out$rotation
pr.out$x
plot(pr.out$x)
pr.out = prcomp(scale(mData), scale=F, center = F)
plot(pr.out$x)

######### lets try a bayesian approach
### first get the estimates using stan and mcmc
library(rstan)
stanDso = rstan::stan_model(file='pcaAndVectors/bayesianPCA.stan')

lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=ncol(mData),
                 y=t(mData))
fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=2)
print(fit.stan)

lResults = extract(fit.stan)
mComp = lResults$mComponents
mComp.1 = mComp[,,1]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,2]
ivComp.2 = colMeans(mComp.2)

plot(ivComp.1, ivComp.2)

mStan = do.call(cbind, extract(fit.stan))
colnames(mStan) = c('beta0', 'beta1', 'sigma')
apply(mStan, 2, mean)
apply(mStan, 2, sd)





library(LearnBayes)

# log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate and data
  # data
  mResp = data$resp
  # fix parameters in initial model
  mEigens = eigen(cov(mResp))$vectors
  ivMu = colMeans(mResp)
  # parameter to track
  mData.rot = rbind(X=theta[grep('X', names(theta))], Y=theta[grep('Y', names(theta))]) 
  iSigma = exp(theta[grep('iSigma', names(theta))])
  
  # calculate fitted value
  mFitted = mEigens %*% mData.rot
  mFitted = sweep(mFitted, 1, ivMu, '+')
  mFitted = t(mFitted)
  # write the priors and likelihood 
  lp1 = sum(apply(mData.rot, 2, function(x) dmnorm(x, mean = c(0, 0), diag(1, 2, 2), log=T)))
  lp2 = lp1 + dcauchy(iSigma, 0, 2.5, log=T)
  lik = sum(sapply(1:nrow(mResp), function(x) sum(dnorm(mResp[x,], mFitted[x,], iSigma, log=T))))
  val = lik + lp2
  return(val)
}

lData = list(resp=mData)
start = c(X=rnorm(nrow(mData), 0, 1), Y=rnorm(nrow(mData), 0, 1), iSigma=log(1))

mylogpost(start, lData)
fit.lap = laplace(mylogpost, start, lData)
tpar = list(m=fit.lap$mode, var=fit.lap$var*2, df=3)
s = sir(mylogpost, tpar, 1000, lData)

library(numDeriv)
mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), data=data)
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

fit.lap = mylaplace(mylogpost, start, lData)











## lets write a custom glm using a bayesian approach
lData = list(resp = mData)
start = c(mEigens=eigen(cov(mData))$vectors,
          mData.rot = scale(t(mData)),
          ivMu = rep(0, 2),
          mSigma = log(cov(mData)))
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  # get eigens
  i = grep('eigens', names(theta))
  mEigens = theta$mEigens
  mData.rot = theta$mData.rot
  ivMu = theta$ivMu
  mSigma = exp(theta$mSigma)
  ## data
  mResp = data$resp # resp
  
  # calculate fitted value
  mFitted = mEigens %*% mData.rot
  mFitted = sweep(mFitted, 1, ivMu, '+')
  mFitted = t(mFitted)
  # write the priors and likelihood 
  lp = sum(dmnorm(t(mData.rot), mean = 0, diag(1, 2, 2), log=T))
  lik = sum(sapply(1:nrow(mResp), function(x) dmnorm(mResp[x,], mFitted[x,], mSigma, log=T)))
  val = lik + lp
  return(val)
}

myl
