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


p.old = par(mfrow=c(3,2))
plot(mData, pch=20, main='raw')
## now standardize the data with 0 mean and variance 1
plot(scale(mData), pch=20, main='Scaled and Centered')

## now perform PCA and data whitening
ivMeans = colMeans(mData)
# centered data
mData.s = sweep(mData, 2, ivMeans, '-')

# covariance matrix for the data without scaling
mCov = cov(mData)
lEigens = eigen(mCov)
colnames(lEigens$vectors) = c('Vector1', 'Vector2')

## inputs are the variables i.e. the data, one column at a time
t(mData.s)[,1:6]
## operations/transformation matrix is the matrix of eigen vectors
t(lEigens$vectors)

## get the transformed points after running them through the matrix, i.e. input output system
mData.rotated = t(lEigens$vectors) %*% t(mData.s)
rownames(mData.rotated) = c('Comp1', 'Comp2')
mData.rotated[,1:6]

plot(t(mData.rotated), main='Centered and decorrelated', pch=20)
plot(prcomp(mData, center = T, scale. = T)$x, main='Centered, Scaled and decorrelated', pch=20)

## get the original data back
## rowDataMatrix = (inverse(rowEigenVectors) * rotated Data) + original Means
solve(t(lEigens$vectors)) ## this equals the transpose of the rowEigenVectors matrix
mData.original.s = (solve(t(lEigens$vectors)) %*% mData.rotated)
mData.original.s[,1:6]

## add the mean to un-center the data
mData.original = sweep(mData.original.s, 1, ivMeans, '+')
mData.original[,1:6]
t(mData)[,1:6]

## whiten the data
whiten = function(mData){
  ## center the data
  ivMeans = colMeans(mData)
  # centered data
  mData.s = sweep(mData, 2, ivMeans, '-')
  ## calculate covariance matrix
  mCov = cov(mData)
  ## see bishop 2006 chapter 12 page 568 for formula
  # y = 1/sqrt(L) * t(U) * centered data
  ## get the eigen vectors and values
  lEigens = eigen(mCov)
  L = diag(lEigens$values)
  U = lEigens$vectors
  # invert after taking square root
  Z = solve(sqrt(L))
  Z = Z %*% t(U)
  yn = Z %*% t(mData.s)
  rownames(yn) = colnames(mData)
  return(t(yn))
}

mData.white = whiten(mData)

plot(mData.white, main='Whitening', pch=20, xlab='Var 1', ylab='Var 2')
par(p.old)


## repeat but with scaling and centring as well
mData.s = scale(mData)
mCov = cov(mData.s)
lEigens = eigen(mCov)
colnames(lEigens$vectors) = c('Vector1', 'Vector2')

## get the transformed points after running them through the matrix, i.e. input output system
mData.rotated = t(lEigens$vectors) %*% t(mData.s)
rownames(mData.rotated) = c('Comp1', 'Comp2')
mData.rotated[,1:6]



# perform PCA using the prcomp
# the vectors are in columns
pr.out = prcomp(mData, scale=F, center = T)
pr.out$rotation
pr.out$x
plot(pr.out$x)
pr.out = prcomp(mData, scale=T, center = T)
plot(pr.out$x)

######### lets try a bayesian approach
### first get the estimates using stan and mcmc
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## load the dataset
data("iris")
mData = as.matrix(iris[,-5])

## compile the first stan script
stanDso = rstan::stan_model(file='pcaAndVectors/bayesianPCA.stan')

## set stan data
lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=2,
                 y=scale(mData))

# ## calculate initial values
# initf = function(chain_id = 1) {
#   list(mEigens=pr.out$rotation)
# } 

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=3, cores=3)#, init=initf) 
                    #control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
print(fit.stan)

## extract the results
lResults = extract(fit.stan)
## extract the latent variables
mComp = lResults$mComponents
dim(mComp)
## get the first 2 latent variables for plotting
mComp.1 = mComp[,,1]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,2]
ivComp.2 = colMeans(mComp.2)

ivCols = rainbow(n = nlevels(iris$Species))
ivCols = ivCols[as.numeric(iris$Species)]
par(mfrow=c(2,2))
plot(ivComp.2, ivComp.1, pch=20, col=ivCols, xlab='PC2', ylab='PC1', main='MCMC 1')

pr.out = prcomp(mData, scale=T, center = T)
plot(pr.out$x, pch=20, col=ivCols, main='Analytical')


## try a second approach with an additional parameter for variance of eigen vectors
stanDso2 = rstan::stan_model(file='pcaAndVectors/bayesianPCA2.stan')

## set stan data
lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=2,
                 y=mData)

fit.stan2 = sampling(stanDso2, data=lStanData, iter=5000, chains=3, cores=3)#, #init=initf, 
                    #control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
print(fit.stan2)

## extract the results
lResults = extract(fit.stan2)
## extract the variances
mSigma2 = lResults$sigma2
summary(mSigma2)
## extract the latent variables
mComp = lResults$mComponents
dim(mComp)
## get the first 2 latent variables for plotting
mComp.1 = mComp[,,1]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,2]
ivComp.2 = colMeans(mComp.2)
plot(ivComp.1, ivComp.2, pch=20, col=ivCols, xlab='PC1', ylab='PC2', main='MCMC 2')

## try a third approach with a different way to sample from multi_normal distribution
stanDso3 = rstan::stan_model(file='pcaAndVectors/bayesianPCA3.stan')

## set stan data
lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=2,
                 y=scale(mData))

fit.stan3 = sampling(stanDso3, data=lStanData, iter=5000, chains=3, cores=3)#, #init=initf, 
                     #control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
print(fit.stan3)

## extract the results
lResults = extract(fit.stan3)
## extract the variances
mSigma2 = lResults$sigma2
summary(mSigma2)
## extract the latent variables
mComp = lResults$mComponents
dim(mComp)
## get the first 2 latent variables for plotting
mComp.1 = mComp[,,3]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,1]
ivComp.2 = colMeans(mComp.2)
plot(ivComp.1, ivComp.2, pch=20, col=ivCols, xlab='PC1', ylab='PC2', main='MCMC 3')


################ repeat with a high dimensional dataset
## load a high dimension data set with zeros
load(file.choose())
mData = t(lData_first$data)
dim(mData)

## try a second approach with an additional parameter for variance of eigen vectors
stanDso2 = rstan::stan_model(file='pcaAndVectors/bayesianPCA2.stan')

## set stan data
lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=min(dim(mData)),
                 y=mData)

fit.stan2 = sampling(stanDso2, data=lStanData, iter=10000, chains=2, cores=2)#, #init=initf, 
#control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
print(fit.stan2)

## extract the results
lResults = extract(fit.stan2)
## extract the variances
mSigma2 = lResults$sigma2
summary(mSigma2)
iSigma2 = colMeans(mSigma2)

## plot the scale parameters 
iOrder = order(iSigma2, decreasing = T)
# par(mfrow=c(2,2))
# plot(iSigma2[iOrder], type='b', xlab='Index', ylab='Standard Deviation', main='MCMC SD Eigen Vec')
# pr.out = prcomp(mData, scale=F, center = T)
# plot(pr.out$sdev, type='b', xlab='Index', ylab='Standard Deviation', main='Analytical SD')

## extract the latent variables
mComp = lResults$mComponents
dim(mComp)
## get the first 2 latent variables for plotting
mComp.1 = mComp[,,iOrder[1]]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,iOrder[2]]
ivComp.2 = colMeans(mComp.2)

ivCols = rainbow(n = nlevels(lData_first$batch))
ivCols = ivCols[as.numeric(lData_first$batch)]
par(mfrow=c(1,2))
plot(1*ivComp.2, 1*ivComp.1, pch=20, col=ivCols, xlab='PC1', ylab='PC2', main='MCMC 2')
text(1*ivComp.2, 1*ivComp.1, labels = rownames(mData), adj =c(0.5, -0.5) , cex = 0.6)
pr.out = prcomp(mData, scale=F, center = T)
plot(pr.out$x, col=ivCols, pch=20, main='Analytical')
text(pr.out$x, labels = rownames(mData), adj = 1.5, cex = 0.6)









### remove this 4th approach 
# ## try a 4th approach but not using multi_normal 
# stanDso4 = rstan::stan_model(file='pcaAndVectors/bayesianPCA4.stan')
# 
# ## set stan data
# lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=2,
#                  y=scale(mData))
# 
# fit.stan4 = sampling(stanDso4, data=lStanData, iter=5000, chains=3, cores=3)#, #init=initf, 
#                      #control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
# print(fit.stan4)
# 
# ## extract the results
# lResults = extract(fit.stan4)
# ## extract the variances
# mSigma2 = lResults$sigma2
# summary(mSigma2)
# ## extract the latent variables
# mComp = lResults$mComponents
# dim(mComp)
# ## get the first 2 latent variables for plotting
# mComp.1 = mComp[,,1]
# ivComp.1 = colMeans(mComp.1)
# mComp.2 = mComp[,,2]
# ivComp.2 = colMeans(mComp.2)
# plot(ivComp.1, ivComp.2, pch=20)



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
