# Name: missingCovariateMixtureModel.R
# Auth: uhkniazi
# Date: 20/03/2018
# Desc: difference in results when using a mixture model or simple linear model when a covariate is missing

####### data loading
dfData = read.csv('BayesRegression/missingCovariateMixtureModelData.csv', header=T)
str(dfData)
dfData$fGroups = factor(dfData$fGroups)

library(lattice)
densityplot(~ x | fGroups, data=dfData)
densityplot(~ x, data=dfData, groups=fGroups, auto.key=T)

## fit a linear regression model first
fit.1 = lm(x ~ fGroups, data=dfData)
summary(fit.1)


## fit a mixture model 
library(flexmix)
fit.flex = flexmix(x ~ fGroups, data=dfData, k=2)
summary(fit.flex)
## fitted coefficients
parameters(fit.flex)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='BayesRegression/normResponseFiniteMixture1RandomEffect.stan')

# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$x), 2*sd(dfData$x))
  
### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$fGroups),
                 NgroupMap1=as.numeric(dfData$fGroups), 
                 y=dfData$x, iMixtures=2, 
                 gammaShape=l$shape, gammaRate=l$rate, 
                 iIntercepts=c(10, 15))
  
fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, 
                    pars=c('sigmaRan1', 'sigma', 'rGroupsJitter1', 'mu', 'iMixWeights', 'muFitted'),
                    cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))

print(fit.stan, c('sigmaRan1', 'sigma', 'mu', 'iMixWeights'), digits=3)

traceplot(fit.stan, c('sigmaRan1', 'sigma', 'mu', 'iMixWeights'))

## check if labelling degeneracy has occured
## see here: http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
params1 = as.data.frame(extract(fit.stan, permuted=FALSE)[,1,])
params2 = as.data.frame(extract(fit.stan, permuted=FALSE)[,2,])

## check if the means from different chains overlap
## Labeling Degeneracy by Enforcing an Ordering
par(mfrow=c(2,2))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
plot(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)

plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
points(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)

## extract the coefficients from the model
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(d=mean(d), z=z, p=p))
}

colnames(mCoef) = levels(dfData$fGroups)

getDifference(mCoef[,'1'], mCoef[,'0'])


################# find the batch using kmeans
km = kmeans(dfData$x, centers = 2)

dfData$cluster = factor(km$cluster)

densityplot(~ x | cluster, groups=fGroups, data=dfData)

fit.2 = lm(x ~ fGroups + cluster, data=dfData)
summary(fit.2)

# #### fit mixed effect model
library(lme4)
fit.lme1 = lmer(x ~ 1 + fGroups + (1 | cluster), data=dfData, REML=F)
summary(fit.lme1)




