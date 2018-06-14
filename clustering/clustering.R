# Name: clustering.R
# Auth: uhkniazi
# Date: 7/6/2018
# Desc: playing with some clustering methods and multiple datasets

library(iCluster)

data("breast.chr17")
?breast.chr17

mData.rna = breast.chr17$mRNA.data
mData.dna = breast.chr17$DNA.data

################################ some basic exploratory analysis of the data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

mData = t(mData.rna)
dim(mData)
colnames(mData)
fBatch = gsub('(\\w+)\\..+', '\\1', colnames(mData))
fBatch[!(fBatch %in% c('NORWAY', 'STANFORD'))] = 'Lines'
fBatch = factor(fBatch)

oDiag.1 = CDiagnosticPlots(mData, 'RNA')

boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.5, cex.main=0.7)

## second data set
mData = t(mData.dna)
dim(mData)
colnames(mData)
fBatch = gsub('(\\w+)\\..+', '\\1', colnames(mData))
fBatch[!(fBatch %in% c('NORWAY', 'STANFORD'))] = 'Lines'
fBatch = factor(fBatch)

oDiag.2 = CDiagnosticPlots(mData, 'DNA')

boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.2, fBatch, cex.main=1)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.5, cex.main=0.7)

### cluster both data sets together using iCluster
fit=iCluster(breast.chr17, k=4, lambda=c(0.2,0.2))
plotiCluster(fit=fit, label=rownames(breast.chr17[[2]]))

######## try clustering the data using PCA model

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

mData.white = whiten(t(mData.dna))

oDiag.3 = CDiagnosticPlots(mData.white, 'white dna')

boxplot.median.summary(oDiag.3, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.3, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.3, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.3, fBatch, cex.main=1)
plot.dendogram(oDiag.3, fBatch, labels_cex = 0.5, cex.main=0.7)


######### lets try a bayesian approach
### first get the estimates using stan and mcmc
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mData = mData.dna

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
save(fit.stan, file='Temp/fit.stan.rds')
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

ivCols = rainbow(n = nlevels(fBatch))
ivCols = ivCols[as.numeric(fBatch)]
par(mfrow=c(1,2))
plot(ivComp.2, ivComp.1, pch=20, col=ivCols, xlab='PC2', ylab='PC1', main='MCMC 1')
text(ivComp.2, ivComp.1, labels = rownames(mData), adj = 1, cex = 0.6)
pr.out = prcomp(mData, scale=T, center = T)
plot(pr.out$x, col=ivCols, pch=20, main='Analytical')
text(pr.out$x, labels = rownames(mData), adj = 1, cex = 0.6)


## try a second approach with an additional parameter for variance of eigen vectors
stanDso2 = rstan::stan_model(file='pcaAndVectors/bayesianPCA2.stan')

## set stan data
lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=min(dim(mData)),
                 y=scale(mData))

fit.stan2 = sampling(stanDso2, data=lStanData, iter=5000, chains=3, cores=3)#, #init=initf, 
#control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
save(fit.stan2, file='Temp/fit.stan2.rds')
print(fit.stan2)

## extract the results
lResults = extract(fit.stan2)
## extract the variances
mSigma2 = lResults$sigma2
summary(mSigma2)
iSigma2 = colMeans(mSigma2)

## plot the scale parameters 
iOrder = order(iSigma2, decreasing = T)
par(mfrow=c(1,2))
plot(iSigma2[iOrder], type='b', xlab='Index', ylab='Standard Deviation', main='MCMC SD Eigen Vec')
pr.out = prcomp(mData, scale=T, center = T)
plot(pr.out$sdev, type='b', xlab='Index', ylab='Standard Deviation', main='Analytical SD')

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
pr.out = prcomp(mData, scale=T, center = T)
plot(pr.out$x, col=ivCols, pch=20, main='Analytical')
text(pr.out$x, labels = rownames(mData), adj = 1, cex = 0.6)

library(scatterplot3d)
## get the first 2 latent variables for plotting
mComp.1 = mComp[,,iOrder[1]]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,iOrder[2]]
ivComp.2 = colMeans(mComp.2)
mComp.3 = mComp[,,iOrder[3]]
ivComp.3 = colMeans(mComp.3)

s3d = scatterplot3d(ivComp.1, ivComp.2, ivComp.3, color = ivCols, pch=20)
text(s3d$xyz.convert(ivComp.1, ivComp.2, ivComp.3), labels = rownames(mData), 
     cex= 0.7, col = "steelblue")

m = cbind(ivComp.1, ivComp.2, ivComp.3)
d = dist(m)
d2 = as.matrix(d)
dim(d2)
image(1:ncol(d2), 1:ncol(d2), d2, axes=F, xlab='', ylab='', col=heat.colors(3))
axis(1, 1:ncol(d2), colnames(d2), cex.axis=0.5, las=3)
axis(2, 1:ncol(d2), colnames(d2), cex.axis=0.4, las=1)




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