# Name: biomarkers.R
# Auth: uhkniazi
# Date: 29/03/2018
# Desc: variable selection for biomarker discovery in TB dataset


# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

setwd('biomarkers/')

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')


lData.train = f_LoadObject('GSE19491_normalised_train.RList')
str(lData.train)

fGroups = rep('ATB', times=length(lData.train$grouping))
fGroups[lData.train$grouping != 'ATB'] = 'Other'
fGroups = factor(fGroups, levels = c('Other', 'ATB'))
table(fGroups)
table(lData.train$grouping)

## select a smaller subset of the data 
###### fit a model using limma
library(limma)

design = model.matrix(~ fGroups)
head(design)

fit = lmFit(lData.train$data, design)
fit = eBayes(fit)

dfLimmma = topTable(fit, coef=2, adjust='BH', number = Inf)
dfLimmma = dfLimmma[order(dfLimmma$adj.P.Val, decreasing = F), ]

cvTopVariables = rownames(dfLimmma)[1:2000]

dfData = data.frame(t(lData.train$data[cvTopVariables, ]))

## perform nested random forest
## adjust boot.num as desired
oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100, big.warn = F)
save(oVar.r, file='oVar.r.rds')

plot.var.selection(oVar.r)

############## use a binomial regression approach this time to rank variables
dfData = data.frame(scale(t(lData.train$data[cvTopVariables, ])))
dim(dfData)
dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == 'ATB', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

save(fit.stan, file='fit.stan.binom.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]

## function to calculate statistics for a coefficient
getDifference = function(ivData){
  # get the difference vector
  d = ivData
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(p)
}

ivPval = apply(mCoef, 2, getDifference)
hist(ivPval)
plot(colMeans(mCoef), ivPval, pch=19)
m = colMeans(mCoef)
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
m = abs(m)
m = sort(m, decreasing = T)
cvTopGenes.binomial = names(m)[1:30] #names(m[m > 0.25])
cvTopGenes.binomial[11] = "HLA-DPA1"

### test performance of both results, from Random Forest and Binomial Regression
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 genes to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes, ]))

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# use the top 30 genes to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes.binomial, ]))

oVar.sub2 = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub2)

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

### try a combination of top genes from both models
cvTopGenes.comb = NULL;
for (i in 1:10){
  cvTopGenes.comb = append(cvTopGenes.comb, CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)) 
  cat(i)
}
cvTopGenes.comb = unique(cvTopGenes.comb)

for (i in 1:10){
  cvTopGenes.comb = append(cvTopGenes.comb, CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i))
  cat(i)
}

cvTopGenes.comb = unique(cvTopGenes.comb)
length(cvTopGenes.comb)

# use these combined variables to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes.comb, ]))

oVar.subComb = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.subComb)

for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.subComb, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}


## load the test data and try these combinations
lData.test = f_LoadObject('GSE37250_normalised_subset_test.RList')
fGroups.test = rep('ATB', times=length(lData.test$grouping))
fGroups.test[lData.test$grouping != 'ATB'] = 'Other'
fGroups.test = factor(fGroups.test, levels = c('Other', 'ATB'))
table(fGroups.test)
table(lData.test$grouping)

dfData = data.frame(t(lData.test$data))
dim(dfData)

## subset the data first into test and training tests
test = sample(1:length(fGroups.test), size = 0.30*length(fGroups.test), replace = F)
table(fGroups.test[test])
table(fGroups.test[-test])

#pdf('cvFigures.pdf')
## 10 fold nested cross validation with various variable combinations
par(mfrow=c(1,2))
for (i in c(1,7)){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = data.frame(dfData[-test ,cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  dfData.test = data.frame(dfData[test ,cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups.test[test],
                             train.groups = fGroups.test[-test], level.predict = 'ATB', boot.num = 100)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

## 10 fold nested cross validation with various variable combinations
for (i in c(1,7)){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i)
  dfData.train = data.frame(dfData[-test ,cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  dfData.test = data.frame(dfData[test ,cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups.test[test],
                             train.groups = fGroups.test[-test], level.predict = 'ATB', boot.num = 100)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

## 10 fold nested cross validation with various variable combinations
for (i in c(1,7)){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.subComb, i)
  dfData.train = data.frame(dfData[-test ,cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  
  dfData.test = data.frame(dfData[test ,cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups.test[test],
                             train.groups = fGroups.test[-test], level.predict = 'ATB', boot.num = 100)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

#dev.off(dev.cur())

###################################
### test a 1 variable model with reject and accept regions
##################################
library(LearnBayes)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  return(iFitted)
}

## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix

  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood
  lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

dfData = data.frame(dfData[ , CVariableSelection.ReduceModel.getMinModel(oVar.sub2, 1)])
colnames(dfData) = CVariableSelection.ReduceModel.getMinModel(oVar.sub2, 1)
dim(dfData)
head(dfData)
dfData = data.frame(dfData, fGroups=fGroups.test)

lData = list(resp=ifelse(dfData$fGroups == 'ATB', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
start = c(rep(0, times=ncol(lData$mModMatrix)))
mylogpost(start, lData)

fit.2 = laplace(mylogpost, start, lData)
fit.2

fit.1 = glm(fGroups ~ ., data=dfData, family='binomial')
data.frame(coef(fit.1), fit.2$mode)
se = sqrt(diag(fit.2$var))
summary(fit.1)
## lets take a sample from this
# parameters for the multivariate t density
tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
## get a sample directly and using sir (sampling importance resampling with a t proposal density)
s = sir(mylogpost, tpar, 5000, lData)
colnames(s) = colnames(lData$mModMatrix)
apply(s, 2, mean)
apply(s, 2, sd)
pairs(s, pch=20)
fit.2$sir = s

stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=3000, chains=3, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', CVariableSelection.ReduceModel.getMinModel(oVar.sub2, 1))
pairs(mCoef, pch=20)

library(lattice)
## get the predicted values
dfData.new = dfData
str(dfData.new)
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
ivPredict = mypred(colMeans(mCoef), list(mModMatrix=X))
xyplot(ivPredict ~ fGroups.test, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)')
xyplot(ivPredict ~ lData.test$grouping, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)')
## choose an appropriate cutoff for accept and reject regions
ivTruth = fGroups.test == 'ATB'

p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')

ivTruth = !ivTruth
p = prediction(1-ivPredict, ivTruth)
perf.death = performance(p, 'tpr', 'fpr')
dfPerf.death = data.frame(c=perf.death@alpha.values, t=perf.death@y.values[[1]], f=perf.death@x.values[[1]], 
                          r=perf.death@y.values[[1]]/perf.death@x.values[[1]])
colnames(dfPerf.death) = c('c', 't', 'f', 'r')

plot(perf.alive)
plot(perf.death, add=T, col='red')
legend('bottomright', legend = c('Alive', 'Dead'), col = 1:2, lty=1)

fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.73062000] = 'ATB'
fPredict[ivPredict <= (1-0.94426070)] = 'Other'
table(fPredict, fGroups.test)

## draw these accept reject points
xyplot(ivPredict ~ fGroups.test, xlab='Actual Group', ylab='Predicted Probability of Being Alive (1)', groups=fPredict,
       auto.key = list(columns=3))



## fit a binomial model
fit.bin = glm(fGroups ~ ., data=dfData, family='binomial')
summary(fit.bin)
ivPredict.bin = predict(fit.bin, type = 'response')

m = data.frame(round(ivPredict, 2), round(ivPredict.bin, 2), fGroups)

## find the optimal point
ivTruth = fGroups == '1'

p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')


ivTruth = !ivTruth
p = prediction(1-ivPredict, ivTruth)
perf.death = performance(p, 'tpr', 'fpr')
dfPerf.death = data.frame(c=perf.death@alpha.values, t=perf.death@y.values[[1]], f=perf.death@x.values[[1]], 
                          r=perf.death@y.values[[1]]/perf.death@x.values[[1]])
colnames(dfPerf.death) = c('c', 't', 'f', 'r')

plot(perf.alive)
plot(perf.death, add=T, col='red')

dfPlot = data.frame(fGroups, ivPredict)

xyplot(ifelse(fGroups == '1', 1, 0) ~ ivPredict, type=c('p'), ylab='Survived')

plot(ivPredict, ifelse(fGroups == '1', 1, 0), type=c('p'), ylab='Class Prediction', pch=20)
#lines(lowess(dfPlot$ivPredict[dfPlot$tpr > 0.9], ifelse(dfPlot$fGroups[dfPlot$tpr > 0.9] == '1', 1, 0)))
lines(lowess(dfPlot$ivPredict, ifelse(dfPlot$fGroups == '1', 1, 0)))

points(ivPredict, ifelse(fGroups == '0', 1, 0), type=c('p'), col='red', pch=20)
lines(lowess(ivPredict, ifelse(fGroups == '0', 1, 0)), col='red')
plot(1-ivPredict, ifelse(fGroups == '0', 1, 0))






