# Name: twoLevelBinomialHierarchicalModel.R
# Auth: uhkniazi
# Date: 09/11/2017
# Desc: model from [1] Kruschke, J. K. (2014). Doing Bayesian data analysis: Chapter 9


library(LearnBayes)

## the data comes from a Bisulphite sequencing experiment 
## we look at the number of methylated cytosine bases at a certain position in the genome
## and the total number of coverage for that region.

dfData = read.csv('sampleData.csv', header=T)
str(dfData)

## our interest is to look at the effect of group1 and group2, while controlling for gender
dfData$interaction = factor(dfData$group1:dfData$group2)
xtabs(~ interaction + group3, data=dfData)

## drop row 8 with missing data
dfData = dfData[-8,]

## we try and do this 2 ways, model using STAN and using a log posterior function
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='twoLevelBinomialHierarchicalModel.stan')

## set up starting data
lStanData = list(Ntotal=nrow(dfData), NgroupsLvl1=nlevels(dfData$group1), 
                 NgroupsLvl2=nlevels(dfData$group1:dfData$group2),
                 NgroupsLvl2Map=c(1, 1, 2, 2, 3, 3),
                 NdataMap = as.numeric(dfData$group1:dfData$group2),
                 y=dfData$methylated, 
                 N=dfData$total)

fit.stan = sampling(stanDso, data=lStanData, iter=4000, chains=4, pars=c('betas', 'sigmaRan1', 'sigmaRan2',
                                                                         'sigmaPop','nu', 'mu', 'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=4)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)
