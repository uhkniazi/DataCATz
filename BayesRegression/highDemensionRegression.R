# Name: highDemensionRegression.R
# Auth: uhkniazi
# Date: 26/01/2018
# Desc: Perform high dimension regression in one model using hierarchical approach


# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}
# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}



lData = f_LoadObject(file.choose())
m = lData$data
## remove missing data vectors
dim(m)
table(rowSums(m) <= 0)
f = rowSums(m) <= 0
m = m[!f,]

dfData = data.frame(t(m))
dfData = stack(dfData)
dfData$fBatch = lData$batch
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData = dfData[order(dfData$Coef), ]
#### fit mixed effect model
library(lme4)
fit.lme1 = lmer(values ~ 1 + (1 | Coef), data=dfData, REML=F)
summary(fit.lme1)

plot(fitted(fit.lme1), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess(fitted(fit.lme1), resid(fit.lme1)), col=2)

## fit model with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponse2RandomEffectsNoFixed.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$values), 2*sd(dfData$values))

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 #Nclusters2=nlevels(dfData$Patient.ID),
                 NgroupMap1=as.numeric(dfData$Coef),
                 #NgroupMap2=as.numeric(dfData$Patient.ID),
                 Ncol=1, 
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, 
                    pars=c('betas', 'sigmaRan1', #'sigmaRan2',
                           'nu', 'sigmaPop','mu', 
                           'rGroupsJitter1'), #'rGroupsJitter2'),
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaPop', 'nu'), digits=3)

## model checks
m = extract(fit.stan, 'mu')
names(m)
dim(m$mu)
fitted = apply(m$mu, 2, mean)

plot(dfData$values, fitted, pch=20, cex=0.5)
plot(dfData$values, dfData$values - fitted, pch=20, cex=0.5)
iResid = (dfData$values - fitted)

par(mfrow=c(1,2))
plot(fitted, iResid, pch=20, cex=0.5, main='t model')
lines(lowess(fitted, iResid), col=2, lwd=2)

plot(predict(fit.lme1), resid(fit.lme1), pch=20, cex=0.5, main='normal')
lines(lowess(predict(fit.lme1), resid(fit.lme1)), col=2, lwd=2)

plot(fitted, predict(fit.lme1), pch=20, cex=0.5)

### plot the posterior predictive values
m = extract(fit.stan, c('mu', 'nu', 'sigmaPop'))
i = sample(1:5000, 5000)
muSample = m$mu[i,]
nuSample = m$nu[i]
sigSample = m$sigmaPop[i]

## t sampling functions
dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu

## use point-wise predictive approach to sample a new value from this data
ivResp = dfData$values
mDraws = matrix(NA, nrow = length(ivResp), ncol=2000)

# rppd = function(index){
#   f = muSample[,index]
#   return(rt_ls(length(f), nuSample, f, sigSample))
# }

for (i in 1:ncol(mDraws)){
  mDraws[,i] = rt_ls(length(ivResp), nuSample[i], muSample[i,], sigSample[i])
}

# 
# temp = sapply(1:length(ivResp), function(x) rppd(x))
# mDraws = t(temp)

yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)#, ylim=c(0, 1))
temp = apply(mDraws, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})