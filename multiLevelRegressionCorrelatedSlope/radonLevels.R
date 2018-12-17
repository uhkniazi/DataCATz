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
     ylab='Counties', xaxt='n')
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
     ylab='Counties', xaxt='n')
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

