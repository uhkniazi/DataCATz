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

## we try and do this 2 ways, model using STAN and using a log posterior function
