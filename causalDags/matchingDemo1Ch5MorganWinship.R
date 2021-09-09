# Name: matchingDemo1Ch5MorganWinship.R
# Auth: uhkniazi
# Date: 06/09/2021
# Desc: example for matching demonstration 1 P145 Morgan & Winship 
#       book on causal models

# simulate the data for matching demonstration 1 page 145


# propensity score function
score = function(S){
  if (S == 1) {
    return(0.18)
  }
  else if (S == 2) {
    return(0.5)
  } else if (S == 3){
    return(0.625)
  }
  else {stop('unknown value')}
}

# generate a data vector
set.seed(123)
S = sample(c(1, 2, 3), size = 100, replace = T, prob = c(0.44, 0.24, 0.32))

# counter-factual data
y_0 = function(S){
  if (S == 1) {
    return(2)
  }
  else if (S == 2) {
    return(6)
  } else if (S == 3){
    return(10)
  }
  else {stop('unknown value')}
}

y_1 = function(S){
  if (S == 1) {
    return(4)
  }
  else if (S == 2) {
    return(8)
  } else if (S == 3){
    return(14)
  }
  else {stop('unknown value')}
}

y1 = sapply(S, y_1) 
y0 = sapply(S, y_0)

## select treatment variable by using propensity score
## more S3 have been assigned treatment i.e P(D=1 | S=3) = 0.625 
ps = sapply(S, score)
d = rbinom(100, 1, ps)
table(S, d)
# generate observed data
y = y1*(d) + (1-d)*y0

mean(y[d == 1])
mean(y[d == 0])
# naive estimator E[y | d = 1] - E[y | d = 0]
mean(y[d == 1]) - mean(y[d == 0])

## the same result after weighting the data by the
## conditional factor S
# under control state after weighting for conditional factor S 
# P(S=1 | D=0) * Y etc
0.6*2 + 0.2*6 + 0.2*10 # E(Y0|D=0) (after marginalising S)
# under treatment state 
# P(S=1| D=1) * Y etc
0.2*4 + 0.3*8 + 0.5*14 # E(Y1|D=1)
# naive estimator 
10.2 - 4.4

## now weight the differences between D=1 and D=0
## ATT
mean(y[d==1 & S==1])*0.2 - mean(y[d==0 & S==1]) * 0.2 # P(S=1|D=1)=P(S,D)/P(D) = 0.08/0.4
mean(y[d==1 & S==2])*0.3 - mean(y[d==0 & S==2]) * 0.3 # P(S=2|D=1)
mean(y[d==1 & S==3])*0.5 - mean(y[d==0 & S==3]) * 0.5 # P(S=3|D=1)
0.4 + 0.6 + 2

## ATC
mean(y[d==1 & S==1])*0.6 - mean(y[d==0 & S==1]) * 0.6 # P(S=1|D=0)
mean(y[d==1 & S==2])*0.2 - mean(y[d==0 & S==2]) * 0.2 # P(S=2|D=0)
mean(y[d==1 & S==3])*0.2 - mean(y[d==0 & S==3]) * 0.2 # P(S=3|D=0)
1.2 + 0.4 + 0.8

## ATE
mean(y[d==1 & S==1])*0.44 - mean(y[d==0 & S==1]) * 0.44 # P(S=1)
mean(y[d==1 & S==2])*0.24 - mean(y[d==0 & S==2]) * 0.24 # P(S=2)
mean(y[d==1 & S==3])*0.32 - mean(y[d==0 & S==3]) * 0.32 # P(S=3)
0.88 + 0.48 + 1.28
## or weight the ATT and ATC by P(D)
3*0.4 + 2.4*0.6 

## try regression technique
dfData = data.frame(y, d=factor(d), S=factor(S), y1, y0, p=ps)
summary(dfData)

f1 = lm(y ~ d, data=dfData)
summary(f1)

f2 = lm(y ~ d + S, data=dfData)
summary(f2)

f3 = lm(y ~ d + S, data=dfData, weights=1-p)
summary(f3)

f4 = lm(y ~ d + S, data=dfData, weights=p)
summary(f4)

library(arm)
matches = matching(z = d, score = dfData$p) 
matches$match.ind
matches$cnts
matches$pairs
matched = as.logical(matches$cnts)
#rownames(dfData) = 1:nrow(dfData)
dfData.matched = dfData[matched,]
dim(dfData.matched)

f5 = lm(y ~ d, data=dfData.matched)
summary(f5)

f6 = lm(y ~ d + S, data=dfData.matched)
summary(f6)

f7 = lm(y ~ d, weights = p, data=dfData.matched)
summary(f7)

f8 = lm(y ~ d + S, weights=p, data=dfData.matched)
summary(f8)

f8 = lm(y ~ d + S, weights=1-p, data=dfData.matched)
summary(f8)
