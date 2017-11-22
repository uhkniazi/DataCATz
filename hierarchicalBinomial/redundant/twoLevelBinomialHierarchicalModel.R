# Name: twoLevelBinomialHierarchicalModel.R
# Auth: uhkniazi
# Date: 09/11/2017
# Desc: model from [1] Kruschke, J. K. (2014). Doing Bayesian data analysis: Chapter 9


library(LearnBayes)

## start by analysing the data from the book examples
## the data is from the book
dfData = read.csv('BattingAverage.csv')
## sort the data by the two levels of interest
dfData = dfData[order(dfData$Player, dfData$PriPos),]

## this is a 2 level model, see figure 9.13 page 252 [1]
# level 0 = highest level
# level 1 = population level with multiple groups
# level 2 = lowest level of hierarchy
#               Level 0
#                 |
#       -------------------               LEVEL 1
#       |         |       |
#       L1.1      L1.2    L1.3
# ---------    --------      --------     LEVEL 2
# |       |    |      |      |      |
# L.1.1  L.1.2 L.2.1  L.2.2  L3.1   L3.2
# |
# data







## the data comes from a Bisulphite sequencing experiment 
## we look at the number of methylated cytosine bases at a certain position in the genome
## and the total number of coverage for that region.

dfData = read.csv('sampleData.csv', header=T)
str(dfData)

## our interest is to look at the effect of group1 and group2, while controlling for gender
dfData$interaction = factor(dfData$group1:dfData$group2:dfData$group3)
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

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)
