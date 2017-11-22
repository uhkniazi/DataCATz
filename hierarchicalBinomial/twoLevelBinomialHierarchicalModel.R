# Name: twoLevelBinomialHierarchicalModel.R
# Auth: uhkniazi
# Date: 09/11/2017
# Desc: model from [1] Kruschke, J. K. (2014). Doing Bayesian data analysis: Chapter 9


library(LearnBayes)

## start by analysing the data from the book examples
## the data is from the book
dfData = read.csv('BattingAverage.csv')
## sort the data by the two levels of interest
##dfData = dfData[order(dfData$Player, dfData$PriPos),]

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
## one Level 2 parameter for each data point
str(dfData)

## we try and do this in different ways, first model using STAN
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='twoLevelBinomialHierarchicalModel.stan')

## set up starting data
lStanData = list(Ntotal=nrow(dfData), NgroupsLvl1=nlevels(dfData$PriPos), 
                 NgroupsLvl2Map=as.numeric(dfData$PriPos),
                 y=dfData$Hits,
                 N=dfData$AtBats)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 15))
print(fit.stan, c('omega1', 'kappa1', 'kappa0', 'omega0'), digits=3)

traceplot(fit.stan, 'omega0')
traceplot(fit.stan, 'kappa0')
traceplot(fit.stan, 'omega1')
traceplot(fit.stan, 'kappa1')

mPositions = extract(fit.stan)$omega1
dim(mPositions)
colnames(mPositions) = levels(dfData$PriPos)

## figure 9.14 from [1]
par(mfrow=c(2,2))
hist(mPositions[,'Pitcher'])
hist(mPositions[,'Pitcher'] - mPositions[,'Catcher'])
plot(mPositions[,'Pitcher'], mPositions[,'Catcher'], pch=20)
hist(mPositions[,'Catcher'])

## second panel of the figure
hist(mPositions[,'Catcher'])
hist(mPositions[,'Catcher'] - mPositions[,'1st Base'])
plot(mPositions[,'Catcher'], mPositions[,'1st Base'], pch=20)
hist(mPositions[,'1st Base'])

mPlayers = extract(fit.stan)$theta
colnames(mPlayers) = as.character(dfData$Player)

## figure 9.15 [1]
hist(mPlayers[,'Kyle Blanks'])
hist(mPlayers[,'Kyle Blanks'] - mPlayers[,'Bruce Chen'])
plot(mPlayers[,'Kyle Blanks'], mPlayers[,'Bruce Chen'], pch=20)
hist(mPlayers[,'Bruce Chen'])

## example of shrinkage
## difference between two players 
dfData[c(494, 754),]
## both have same number of trials but different number of successes
## difference includes 0 - due to shrinkage towards the general trend for pitchers
hist(mPlayers[,'Mike Leake'] - mPlayers[,'Wandy Rodriguez'])

## another set of players
dfData[c(573, 428),]
## both have larger number of trials but different successes
## data influences the shrinkage and it is less than in previous case
hist(mPlayers[,'Andrew McCutchen'] - mPlayers[,'Brett Jackson'])





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
