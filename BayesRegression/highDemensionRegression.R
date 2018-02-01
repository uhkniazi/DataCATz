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

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$adj.P.Val < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logFC, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

lData = f_LoadObject(file.choose())
m = lData$data
dim(m)

dfData = data.frame(t(m))
dfData = stack(dfData)
dfData$fBatch = lData$cov1
dfData$fAdjust = lData$cov2
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj = factor(dfData$fAdjust:dfData$ind)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj), ]


###### fit a model using limma
library(limma)

#fBatch = lData$cov1; fAdjust=lData$cov2;
design = model.matrix(~ lData$cov1 + lData$cov2)
colnames(design) = levels(lData$cov1)
head(design)

fit = lmFit(lData$data, design)
fit = eBayes(fit)

dfLimmma.2 = topTable(fit, coef=2, adjust='BH', number = Inf)
dfLimmma.0 = topTable(fit, coef=3, adjust='BH', number = Inf)


# #### fit mixed effect model
library(lme4)
fit.lme1 = lmer(values ~ 1 + (1 | Coef) + (1 | Coef.adj), data=dfData, REML=F)
summary(fit.lme1)
ran = ranef(fit.lme1, condVar=T)

plot(fitted(fit.lme1), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess(fitted(fit.lme1), resid(fit.lme1)), col=2)

## fit model with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='normResponse2RandomEffectsNoFixed.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$values), 2*sd(dfData$values))
# ## set initial values
ran = ranef(fit.lme1)
r1 = ran$Coef
r2 = ran$Coef.adj
initf = function(chain_id = 1) {
  list(sigmaRan1 = 2, sigmaPop=1, rGroupsJitter1=r1, rGroupsJitter2=r2)
}

### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj),
                 Ncol=1, 
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(dfData$values), intercept_sd= sd(dfData$values))

fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2, 
                    pars=c('sigmaRan1', 'sigmaRan2', 'betas',
                           'sigmaPop', 'mu', 
                           'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=2, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop'), digits=3)


## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# ## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan)$betas)
# ## add the intercept to each random effect variable, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='12', deflection='2') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### compare the results from the 3 models
# dfResults = dfResults[order(dfResults$pvalue, decreasing = F), ]
# dfLimmma.2 = dfLimmma.2[order(dfLimmma.2$P.Value, decreasing = F), ]
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
dfLimmma.2$SYMBOL = as.character(rownames(dfLimmma.2))
dfResults$SYMBOL = as.character(rownames(dfResults))

## produce the plots 
f_plotVolcano(dfLimmma.2, 'limma 2M vs 12M')
f_plotVolcano(dfResults, 'Stan 2M vs 12M')

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(m, dfResults, 0.1, 'Stan 2M vs 12M')

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfLimmma.2), names(m))
m = m[i]
m = m[!is.na(m)]
i = match(names(m), rownames(dfLimmma.2))
dfLimmma.2 = dfLimmma.2[i,]
identical(names(m), rownames(dfLimmma.2))
plotMeanFC(m, dfLimmma.2, 0.1, 'limma 2M vs 12M')


i = match(rownames(dfResults), rownames(dfLimmma.2))
dfLimmma.2 = dfLimmma.2[i,]
i = match(rownames(dfLimmma.2), rownames(dfResults))
dfResults = dfResults[i,]
identical(rownames(dfResults), rownames(dfLimmma.2))

plot(dfResults$pvalue, dfLimmma.2$P.Value, pch=20, cex=0.6, col='grey')
plot(dfResults$adj.P.Val, dfLimmma.2$adj.P.Val, pch=20, cex=0.6, col='grey')
df = cbind(stan=dfResults$pvalue, limma=dfLimmma.2$P.Value)

write.csv(dfResults, file='temp/stan.csv', row.names = F)

################# t model
stanDso = rstan::stan_model(file='tResponse1RandomEffectsNoFixed.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$values), 2*sd(dfData$values))

## set initial values
# initf = function(chain_id = 1) {
#   gm = tapply(dfData$values, dfData$Coef, mean) - mean(dfData$values)
#   list(betas = mean(dfData$values), sigmaRan1 = sd(gm), sigmaPop=sd(dfData$values), nu=4, rGroupsJitter1=gm)
# }
## set initial values
ran = ranef(fit.lme1)
ran = ran$Coef
initf = function(chain_id = 1) {
  list(betas = 5.4, sigmaRan1 = 2, sigmaPop=0.5, rGroupsJitter1=ran)
}



### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 #Nclusters2=nlevels(dfData$Patient.ID),
                 NgroupMap1=as.numeric(dfData$Coef),
                 #NgroupMap2=as.numeric(dfData$Patient.ID),
                 Ncol=1,
                 y=dfData$values,
                 gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
                    pars=c('betas', 'sigmaRan1', #'sigmaRan2',
                           'nu', 'sigmaPop', #'mu',
                           'rGroupsJitter1'), #'rGroupsJitter2'),
                    cores=2, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaPop', 'nu'), digits=3)




# ## model checks
# m = extract(fit.stan, 'mu')
# names(m)
# dim(m$mu)
# fitted = apply(m$mu, 2, mean)
# 
# plot(dfData$values, fitted, pch=20, cex=0.5)
# plot(dfData$values, dfData$values - fitted, pch=20, cex=0.5)
# iResid = (dfData$values - fitted)
# 
# par(mfrow=c(1,2))
# plot(fitted, iResid, pch=20, cex=0.5, main='t model')
# lines(lowess(fitted, iResid), col=2, lwd=2)
# 
# plot(predict(fit.lme1), resid(fit.lme1), pch=20, cex=0.5, main='normal')
# lines(lowess(predict(fit.lme1), resid(fit.lme1)), col=2, lwd=2)
# 
# plot(fitted, predict(fit.lme1), pch=20, cex=0.5)
# 
# ### plot the posterior predictive values
# m = extract(fit.stan, c('mu', 'nu', 'sigmaPop'))
# i = sample(1:5000, 5000)
# muSample = m$mu[i,]
# nuSample = m$nu[i]
# sigSample = m$sigmaPop[i]
# 
# ## t sampling functions
# dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
# rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
# 
# ## use point-wise predictive approach to sample a new value from this data
# ivResp = dfData$values
# mDraws = matrix(NA, nrow = length(ivResp), ncol=2000)
# 
# # rppd = function(index){
# #   f = muSample[,index]
# #   return(rt_ls(length(f), nuSample, f, sigSample))
# # }
# 
# for (i in 1:ncol(mDraws)){
#   mDraws[,i] = rt_ls(length(ivResp), nuSample[i], muSample[i,], sigSample[i])
# }
# 
# # 
# # temp = sapply(1:length(ivResp), function(x) rppd(x))
# # mDraws = t(temp)
# 
# yresp = density(ivResp)
# plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)#, ylim=c(0, 1))
# temp = apply(mDraws, 2, function(x) {x = density(x)
# #x$y = x$y/max(x$y)
# lines(x, col='red', lwd=0.6)
# })