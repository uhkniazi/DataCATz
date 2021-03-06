# Name: radonLevels.R
# Auth: umar niazi 
# Date: 13/12/2018
# Desc: multilevel regression example with correlated intercept and slope from Gelman Multilevel modelling book 2006



# load and format the data
library(LearnBayes)
p.old = par()
# see instructions for data download here
# http://wiki.math.yorku.ca/index.php/Book:_Gelman_%26_Hill_%282007%29#Data_sets
# data formatting
# http://www.stat.columbia.edu/~gelman/arm/examples/radon/radon_setup.R
# Set up the radon data
setwd('multiLevelRegressionCorrelatedSlope/')
## this section is from http://www.stat.columbia.edu/~gelman/arm/examples/radon/radon_setup.R
# read in and clean the data

srrs2 <- read.table ("srrs2.dat", header=T, sep=",")
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
x <- floor

# get county index variable

county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

# get the county-level predictor

srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)
## end section for data loading and formatting from external source
county.name = gsub('^\\s+|\\s+$', '', county.name)
uniq <- unique(county.name)
dfData = data.frame(u.full = u[county], county.name, county,
                    floor, log.radon, radon, uranium.full=uranium[county], x, y, county.f = factor(county))

# Complete-pooling and no-pooling estimates of county radon levels
# Consider the goal of estimating the distribution of
# radon levels of the houses within each of the 85 counties in Minnesota.
# page 252 onwards Gelman 2006 book
## complete pooling
iCompletePooling = mean(dfData$y)

## no-pooling analysis
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='linearRegressionNoPooling.stan')

m = model.matrix(y ~ county.f - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'mu', 'sigmaPop'),
                    cores=2)
print(fit.stan, c('betas', 'sigmaPop'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)

# fitted coefficients
fit.stan.1.noPooling = fit.stan

mCoef = extract(fit.stan)$betas
colnames(mCoef) = uniq

## utility function for plotting
getms = function(f){
  m = mean(f)
  se = sd(f)
  m.up = m+1.96*se
  m.down = m-1.96*se
  ret= c(m, m.up, m.down)
  names(ret) = c('m', 'm.up', 'm.down')
  return(ret)
}

iOrder = tapply(dfData$y, dfData$county.f, length)
mCoef = mCoef[,order(iOrder)]
mCoef.noPooling = mCoef
## format for line plots
df = apply(mCoef, 2, getms)
x = 1:ncol(mCoef)
ylim=c(min(df), max(df))

par(p.old)
plot(x, df['m',], ylim=ylim, pch=20, xlab='', main='Average Radon Levels No Pooling',
     ylab='log Radon', xaxt='n')
axis(1, at = x, labels = colnames(mCoef), las=2, cex.axis=0.7)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(h = iCompletePooling, col='grey')

### partial pooling, multilevel model
stanDso = rstan::stan_model(file='linearRegressionPartialPooling_1.stan')

m = model.matrix(y ~ county.f - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan'),
                    cores=2)
print(fit.stan, c('betas', 'populationMean', 'sigmaPop', 'sigmaRan'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan'), ncol=1, inc_warmup=F)

# fitted coefficients
fit.stan.1.partialPooling = fit.stan

mCoef = extract(fit.stan)$betas
colnames(mCoef) = uniq

iOrder = tapply(dfData$y, dfData$county.f, length)
mCoef = mCoef[,order(iOrder)]
mCoef.partialPooling = mCoef
## format for line plots
df = apply(mCoef, 2, getms)
x = 1:ncol(mCoef)

par(p.old)
plot(x, df['m',], ylim=ylim, pch=20, xlab='', main='Average Radon Levels Partial Pooling',
     ylab='log Radon', xaxt='n')
axis(1, at = x, labels = colnames(mCoef), las=2, cex.axis=0.7)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(h = iCompletePooling, col='grey')

## plot together? todo

## check with lmer
library(lme4)
fit.lmer.1 = lmer(y ~ 1 + ( 1 | county.f), data=dfData)
summary(fit.lmer.1)
print(fit.stan.1.partialPooling, c('populationMean', 'sigmaPop', 'sigmaRan'), digits=3)
print(fit.stan.1.noPooling, c('sigmaPop'), digits=3)

################ adding continuous predictors and unit level predictors
# the general principle remains that multilevel models compromise between pooled and unpooled
# estimates, with the relative weights determined by the sample size in the group and
# the variation within and between groups. [Gelman 2006]

## using the variable x i.e. location in the house (basement=0, first floor=1) where measurement was made
## complete pooling model
fit.2.complete = lm(y ~ x, data=dfData)
summary(fit.2.complete)

## no-pooling of intercepts or country variable but complete pooling of basement variable
m = model.matrix(y ~ county.f + x - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 y=dfData$y)

stanDso = rstan::stan_model(file='linearRegressionNoPooling.stan')
fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'sigmaPop'),
                    cores=2)
print(fit.stan, c('betas', 'sigmaPop'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)

# fitted coefficients
fit.stan.2.noPooling = fit.stan

mCoef = extract(fit.stan)$betas
colnames(mCoef) = c(uniq, 'floor')
mCoef.noPooling = mCoef

par(mfrow=c(2,2))
ylim = c(range(dfData$y))

## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}


# choose some counties - figure 12.2 page 255
df = subset(dfData, dfData$county.name == 'AITKIN')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='AITKIN', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef)
m = m[c('AITKIN', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)

df = subset(dfData, dfData$county.name == 'DOUGLAS')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='DOUGLAS', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef)
m = m[c('DOUGLAS', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)

df = subset(dfData, dfData$county.name == 'CLAY')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='CLAY', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef)
m = m[c('CLAY', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)

df = subset(dfData, dfData$county.name == 'ST LOUIS')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='ST LOUIS', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef)
m = m[c('ST LOUIS', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)

# The no-pooling model includes county indicators, which can change the estimated coefficient for
# x, if the proportion of houses with basements varies among counties. a special case of the rule
# that adding new predictors in a regression can change
# the estimated coefficient of x, if these new predictors are correlated with x.
# [Gelman 2006]

# In the multilevel model, a “soft constraint” is applied to the α j ’s: they are as-
# signed a probability distribution, [Gelman 2006]

# this is a more general script, so the sigmaRan[2]
# will not be estimated well as there is only one beta that
# needs to be estimated for the slope as it is the population level
stanDso = rstan::stan_model(file='linearRegressionPartialPooling_2.stan')

m = model.matrix(y ~ county.f + x - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$county.f)), 2),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan'),
                    cores=2)
print(fit.stan, c('betas', 'populationMean', 'sigmaPop', 'sigmaRan'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan'), ncol=1, inc_warmup=F)

# fitted coefficients
fit.stan.2.partialPooling = fit.stan

mCoef = extract(fit.stan)$betas
colnames(mCoef) = c(uniq, 'floor')
mCoef.partialPooling = mCoef

par(mfrow=c(2,2))
ylim = c(range(dfData$y))

# choose some counties - figure 12.4 page 257
df = subset(dfData, dfData$county.name == 'AITKIN')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='AITKIN', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef.noPooling)
m = m[c('AITKIN', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)
m = colMeans(mCoef.partialPooling)
m = m[c('AITKIN', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m, col=2)



df = subset(dfData, dfData$county.name == 'DOUGLAS')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='DOUGLAS', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef.noPooling)
m = m[c('DOUGLAS', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)
m = colMeans(mCoef.partialPooling)
m = m[c('DOUGLAS', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m, col=2)


df = subset(dfData, dfData$county.name == 'CLAY')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='CLAY', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef.noPooling)
m = m[c('CLAY', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)
m = colMeans(mCoef.partialPooling)
m = m[c('CLAY', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m, col=2)

df = subset(dfData, dfData$county.name == 'ST LOUIS')
plot(jitter.binary(df$x), df$y, ylim=ylim, xlab='Floor', ylab='log Radon', main='ST LOUIS', pch=20)
abline(fit.2.complete, lty=2)
m = colMeans(mCoef.noPooling)
m = m[c('ST LOUIS', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m)
m = colMeans(mCoef.partialPooling)
m = m[c('ST LOUIS', 'floor')]
lines(c(0, 1), cbind(c(1, 1), c(0, 1)) %*% m, col=2)


# calculate variance ratio
print(fit.stan, c('sigmaPop', 'sigmaRan[1]'))

var.county = 0.33^2
var.pop = 0.76^2

var.ratio = round(var.county / var.pop, 2)
# the standard deviation of average radon levels between
# counties is the same as the standard
# deviation of the average of 5 measurements within a county [Gelman 2006]
# for a county with a sample size less than 5, there is more information
# in the group-level model than in the county’s data; for a county with more than 5
# observations, the within-county measurements are more informative (in the sense
# of providing a lower-variance estimate of the county’s average radon level). As a
# result, the multilevel regression line in a county is closer to the complete-pooling
# estimate when sample size is less than 5, and closer to the no-pooling estimate when
# sample size exceeds 5. [Gelman 2006]
iInformationContent = 100/ (var.ratio * 100) # around 1/5
0.76/(sqrt(iInformationContent))

# intra class correlation 
# The relative values of individual- and group-level variances are also sometimes expressed using the intraclass correlation,
# which ranges from 0 if the grouping conveys no information to 1 if
# all members of a group are identical. [Gelman 2006]
var.county/(var.pop + var.county)

## using lmer
fit.lmer.2 = lmer(y ~ 1 + x + ( 1 | county.f), data=dfData)
summary(fit.lmer.2)


#### add the floor predictor as a varying slope without correlation of slopes and intercepts
fit.lmer.3 = lmer(y ~ 1 + x + ( 1 | county.f) + (0 + x | county.f), data=dfData)
summary(fit.lmer.3)

## first use earlier more general script then use a different
## formulation for same model
m = model.matrix(y ~ county.f + x:county.f - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$county.f)), 
                                              rep(2, times=nlevels(dfData$county.f))),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan'),
                    cores=2)
print(fit.stan, c('populationMean', 'sigmaPop', 'sigmaRan'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan'), ncol=1, inc_warmup=F)

fit.stan.3.unCorCoef.1 = fit.stan
mCoef = extract(fit.stan)$betas
intercepts.1 = colMeans(mCoef[,1:85])
slopes.1 = colMeans(mCoef[, 86:ncol(mCoef)])

## try a second forumulation
stanDso = rstan::stan_model(file='linearRegressionPartialPooling_3.stan')

lStanData = list(Ntotal=nrow(dfData), 
                 Ngroup1 = nlevels(dfData$county.f),
                 Ngroup2 = nlevels(dfData$county.f),
                 X=dfData$x,
                 Ngroup1Map=as.numeric(dfData$county.f),
                 Ngroup2Map=as.numeric(dfData$county.f),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('coefGroup1', 'coefGroup2', 
                                                                         'MuPopGrp1', 'MuPopGrp2',
                                                                          'sigmaPop',
                                                                         'sigmaRan1', 'sigmaRan2'),
                    cores=2)
print(fit.stan, c('MuPopGrp1', 'MuPopGrp2', 'sigmaPop', 'sigmaRan1', 'sigmaRan2'), digits=3)

# fitted coefficients
fit.stan.3.unCorCoef.2 = fit.stan

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan1'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan2'), ncol=1, inc_warmup=F)

intercepts.2 = colMeans(extract(fit.stan)$coefGroup1)
slopes.2 = colMeans(extract(fit.stan)$coefGroup2)

# compare with lmer
c = coef(fit.lmer.3)
intercepts.lm = c$county.f[,1]
slopes.lm = c$county.f[,2]

par(mfrow=c(1,2))
plot(intercepts.lm, intercepts.1, pch=20, main='intercepts')
points(intercepts.lm, intercepts.2, pch=20, col=2)

plot(slopes.lm, slopes.1, pch=20, main='slopes')
points(slopes.lm, slopes.2, pch=20, col=2)


##### varying coefficients i.e. intercept and slope with correlation
fit.lmer.4 = lmer(y ~ 1 + x + ( 1 + x| county.f), data=dfData)
summary(fit.lmer.4)


stanDso = rstan::stan_model(file='linearRegressionPartialPooling_4.stan')

lStanData = list(Ntotal=nrow(dfData), 
                 Ngroup1 = nlevels(dfData$county.f),
                 X=dfData$x,
                 Ngroup1Map=as.numeric(dfData$county.f),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('coefGroup1', 
                                                                         'MuPopGrp1', 'rho',
                                                                         'sigmaPop',
                                                                         'sigmaRan1', 'sigmaRan2'),
                    cores=4, control=list(adapt_delta=0.99, max_treedepth = 10))
print(fit.stan, c('MuPopGrp1', 'rho', 'sigmaPop', 'sigmaRan1', 'sigmaRan2'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan1'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan2'), ncol=1, inc_warmup=F)

## try a different forumulation using cholesky factor decomposition of the correlation matrix
stanDso = rstan::stan_model(file='linearRegressionPartialPooling_4.2.stan')

lStanData = list(Ntotal=nrow(dfData), 
                 Ngroup1 = nlevels(dfData$county.f),
                 X=dfData$x,
                 Ngroup1Map=as.numeric(dfData$county.f),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=3000, chains=4, pars=c('coefGroup1_adjusted',
                                                                         'MuPopGrp1',
                                                                         'sigmaPop',
                                                                         'sigmaRan1'),
                    cores=4)
print(fit.stan, c('MuPopGrp1', 'sigmaPop', 'sigmaRan1'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('sigmaRan1'), ncol=1, inc_warmup=F)

# fitted coefficients
fit.stan.4.CorCoef = fit.stan

########## adding a group level predictor - uranium
## this predictor only has one measurement per county, and should effect the variation
## between counties i.e. the variance term sigmaRan
# use the uranium.full variable which expands/maps the uranium variable to map to individual observations
fit.lmer.5 = lmer(y ~ 1 + x + u.full + x:u.full + ( 1 + x| county.f), data=dfData)
summary(fit.lmer.5)

stanDso = rstan::stan_model(file='linearRegressionPartialPooling_5.stan')

lStanData = list(Ntotal=nrow(dfData), 
                 Ngroup1 = nlevels(dfData$county.f),
                 X=dfData$x,
                 Xg=dfData$u.full,
                 Ngroup1Map=as.numeric(dfData$county.f),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('coefGroup1_adjusted',
                                                                         'populationCoeff',
                                                                         'sigmaPop',
                                                                         'sigmaRan1'),
                    cores=4)

print(fit.stan, c('populationCoeff',
                  'sigmaPop',
                  'sigmaRan1'), digits=3)


mAdjustments = extract(fit.stan)$coefGroup1_adjusted
mAdjustments.intercepts = mAdjustments[,1,]
mAdjustments.slopes = mAdjustments[,2,]

# compare with lmer
c = ranef(fit.lmer.5)
intercepts.lm = c$county.f[,1]
slopes.lm = c$county.f[,2]

par(mfrow=c(1,2))
plot(intercepts.lm, colMeans(mAdjustments.intercepts), pch=20, xlab='lmer', ylab='stan', main='comparisons of intercepts')
plot(slopes.lm, colMeans(mAdjustments.slopes), pch=20, xlab='lmer', ylab='stan', main='comparisons of slopes')

# extract population coefficients
mPop = extract(fit.stan)$populationCoeff
colMeans(mPop)
colMeans(mAdjustments.intercepts)[85]
colMeans(mAdjustments.slopes)[85]
# county 85 results
# level 2 model - intercept
cbind(1, 1, u[85]) %*% rbind(1.46, -0.03, 0.8)
# level 2 model - slope
cbind(1, 1, u[85]) %*% rbind(-0.66, -0.001, 1*u[85]*-0.39)

# # intercepts for each group
mInt = mPop[,c(1,2)]
mX = cbind(rep(1, times=nlevels(dfData$county.f)), u)
mInt = mX %*% t(mInt)
iPop.int = colMeans(mPop)[1:2]
grid = seq(range(u)[1], range(u)[2], length.out = 100)
y = cbind(rep(1, times=100), grid) %*% cbind(iPop.int)
## format for line plots
mCoef = t(mInt)
plot(grid, y, type='l', ylim=range((mCoef)), xlab='county log uranium',
     ylab='county intercepts')

for (i in 1:85){
  mCoef[,i] = cbind(1, 1, u[i]) %*% rbind(mPop[,1], mAdjustments.intercepts[,i], mPop[,2])
}
# order according to increasing uranium
i = order(u)
mCoef = mCoef[,i]
df = apply(mCoef, 2, getms)
x = u[i]
points(x, colMeans(mCoef), pch=20)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}

var.county = 0.13^2
var.pop = 0.75^2

var.ratio = round(var.county / var.pop, 2)
iInformationContent = 100/ (var.ratio * 100)
