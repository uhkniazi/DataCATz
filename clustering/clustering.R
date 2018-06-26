# Name: clustering.R
# Auth: uhkniazi
# Date: 7/6/2018
# Desc: playing with some clustering methods and multiple datasets

library(iCluster)

p.old = par()

# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

lData.train = f_LoadObject('~/Data/R/CGraphClust/workflow/results/lData.train.rds')
names(lData.train)

oExp = lData.train$data

## take a sample
set.seed(123)
i = sample(1:537, 60)
oExp = oExp[,i]
dim(oExp)

mCounts = t(exprs(oExp))
dim(mCounts)
str(mCounts)
fGroups = oExp$fSamples
table(fGroups)
rownames(mCounts) = as.character(fGroups)

#### randomly split the data set into 2 parts
## sub select top genes sorted on p-values
cvTopGenes = rownames(lData.train$results)[1:50]
mCounts = mCounts[,cvTopGenes]
dim(mCounts)
set.seed(123)
i = sample(1:ncol(mCounts), 25)
mCounts.1 = mCounts[,i]
mCounts.2 = mCounts[,-i]
str(mCounts.1)
lData = list(mCounts.1, mCounts.2)
### cluster both data sets together using iCluster
fit=iCluster(lData, k=3, lambda=c(0.2,0.2))
plotiCluster(fit=fit, label=rownames(mCounts))

######## try clustering the data using PCA model
ivCols = rainbow(n = nlevels(fGroups))
ivCols = ivCols[as.numeric(fGroups)]
pr.comp = prcomp(mCounts)
plot(pr.comp$x, col=ivCols, pch=20, main='Analytical')
text(pr.comp$x, labels = rownames(mCounts), adj = 1, cex = 0.6)
legend('bottomright', legend = levels(fGroups), fill = rainbow(n = nlevels(fGroups)))
plot(pr.comp$sdev)
######### lets try a bayesian approach
### first get the estimates using stan and mcmc
mData = mCounts
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## try a second approach with an additional parameter for variance of eigen vectors
# stanDso2 = rstan::stan_model(file='pcaAndVectors/bayesianPCA2.stan')
# 
# dim(mData)
# ## set stan data
# lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=7,
#                  y=(mData))

# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

stanDso2 = rstan::stan_model(file='clustering/bayesianPCAwithBatch.stan')

l = gammaShRaFromModeSD(sd(mData), 2*sd(mData))

fGroups = oExp$fHIV

lStanData = list(Ntotal=nrow(mData), Nvars=ncol(mData), Neigens=3,
                 Nclusters1=nlevels(fGroups),
                 NgroupMap1=as.numeric(fGroups),
                 y=(mData),
                 gammaShape=l$shape, gammaRate=l$rate)

pr.out = prcomp(mData, scale=F, center = T)
plot(pr.out$x, col=ivCols, pch=20, main='Analytical')
text(pr.out$x, labels = rownames(mData), adj = 1, cex = 0.6)
plot(pr.out$sdev)

fit.stan2 = sampling(stanDso2, data=lStanData, iter=5000, chains=3, cores=3)#, #init=initf, 
#control=list(adapt_delta=0.99, max_treedepth = 10)) # some additional control values/not necessary usually
save(fit.stan2, file='Temp/fit.stan2.tb.batch.rds')
print(fit.stan2)

## extract the results
lResults = extract(fit.stan2)
names(lResults)
## extract the variances
mSigma2 = lResults$sigma2
summary(mSigma2)
iSigma2 = colMeans(mSigma2)

## plot the scale parameters 
iOrder = order(iSigma2, decreasing = T)
par(mfrow=c(1,2))
plot(iSigma2[iOrder], type='b', xlab='Index', ylab='Standard Deviation', main='MCMC SD Eigen Vec')
pr.out = prcomp(mData, scale=T, center = T)
plot(pr.out$sdev[1:3], type='b', xlab='Index', ylab='Standard Deviation', main='Analytical SD')

## extract the latent variables
mComp = lResults$mComponents
dim(mComp)
## get the first 2 latent variables for plotting
mComp.1 = mComp[,,iOrder[1]]
ivComp.1 = colMeans(mComp.1)
mComp.2 = mComp[,,iOrder[2]]
ivComp.2 = colMeans(mComp.2)

# ivCols = rainbow(n = nlevels(lData_first$batch))
# ivCols = ivCols[as.numeric(lData_first$batch)]
par(mfrow=c(1,2))
plot(1*ivComp.2, 1*ivComp.1, pch=20, col=ivCols, xlab='PC1', ylab='PC2', main='MCMC 2')
text(1*ivComp.2, 1*ivComp.1, labels = rownames(mData), adj =c(0.5, -0.5) , cex = 0.6)
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

par(mfrow=c(1,1))
s3d = scatterplot3d(ivComp.1, ivComp.2, ivComp.3, color = ivCols, pch=20)
text(s3d$xyz.convert(ivComp.1, ivComp.2, ivComp.3), labels = rownames(mData), 
     cex= 0.7, col = "steelblue")

m = cbind(ivComp.1, ivComp.2, ivComp.3)
rownames(m) = rownames(mData)
d = dist(m)
d2 = as.matrix(d)
dim(d2)
d2 = d2[nrow(d2):1,]
image(1:ncol(d2), 1:ncol(d2), d2, axes=F, xlab='', ylab='', col=heat.colors(3))
axis(1, 1:ncol(d2), colnames(d2), cex.axis=0.5, las=3)
axis(2, 1:ncol(d2), rownames(d2), cex.axis=0.4, las=1)

dim(mComp)
mComp.all = lapply(seq_along(1:dim(mComp)[3]), function(x){
  colMeans(mComp[,,iOrder[x]])
})

mComp.all = do.call(cbind, mComp.all)
dim(mComp.all)
m = mComp.all
rownames(m) = rownames(mData)
d = dist(m)
plot(hclust(d))
d2 = as.matrix(d)
dim(d2)
d2 = d2[nrow(d2):1,]
image(1:ncol(d2), 1:ncol(d2), d2, axes=F, xlab='', ylab='', col=heat.colors(3))
axis(1, 1:ncol(d2), colnames(d2), cex.axis=0.2, las=3)
axis(2, 1:ncol(d2), rownames(d2), cex.axis=0.2, las=1)

k = kmeans(m, 3)
k$cluster

fGroups = oExp$fSamples
table(k$cluster, fGroups)
table(k$cluster, oExp$fHIV)
table(k$cluster, oExp$fHIV:fGroups)
c = sort(k$cluster)
c
d2 = d2[names(c), names(c)]
dim(d2)
image(1:ncol(d2), 1:ncol(d2), d2, axes=F, xlab='', ylab='', col=heat.colors(6))
axis(1, 1:ncol(d2), colnames(d2), cex.axis=0.5, las=3)
axis(2, 1:ncol(d2), rownames(d2), cex.axis=0.4, las=1)

oDiag.3 = CDiagnosticPlots(t(m), 'projected')

boxplot.median.summary(oDiag.3, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.3, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.3, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.3, fBatch, cex.main=1)
plot.dendogram(oDiag.3, fBatch, labels_cex = 0.5, cex.main=0.7)

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