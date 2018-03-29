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


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

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
### test performance of both models
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 genes to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes, ]))

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 30)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

## load the test data
lData.test = f_LoadObject('GSE37250_normalised_subset_test.RList')
fGroups.test = rep('ATB', times=length(lData.test$grouping))
fGroups.test[lData.test$grouping != 'ATB'] = 'Other'
fGroups.test = factor(fGroups.test, levels = c('Other', 'ATB'))
table(fGroups.test)
table(lData.test$grouping)
## 10 fold nested cross validation with various variable combinations
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 2:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  dfData.train = data.frame(t(lData.train$data[cvTopGenes.sub, ]))
  
  dfData.test = data.frame(t(lData.test$data[cvTopGenes.sub, ]))
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = fGroups.test,
                             train.groups = fGroups, level.predict = 'ATB', boot.num = 50)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}




