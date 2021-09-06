# Name: weightedRegressionMorganWinship.R
# Auth: uhkniazi
# Date: 25/08/2021
# Desc: calculate the ATE, ATC and ATT using the data from Morgan & Winship 
#       book on causal models, chapters 5 and 7

# data downloaded from https://www.cambridge.org/gb/academic/subjects/sociology/sociology-general-interest/counterfactuals-and-causal-inference-methods-and-principles-social-research-2nd-edition?format=HB&isbn=9781107065079

# import and explore the data
library(foreign)
dfData = read.dta(file.choose())

# simulate the data for matching demonstration 3 page 153
# and page 229

# propensity score function
# score = function(A, B) {
#   r = -2 +3*(A) -3*(A-0.1) +2*(A-.3) -2*(A-.5) +4*(A-.7)
#   -4*(A-.9) +1*(B) -1*(B-.1) +2*(B-.7) -2*(B -.9) +3*(A-.5)*(B -.5)
#   -3*(A-.7)*(B -.7)
#   return(plogis(r))
# }

score = function(A, B){
  return(rbeta(length(A), 10, 10))
}


# create grid of points
A = seq(0.1, 1, length.out=30)
B = seq(0.1, 1, length.out = 30)

z = outer(A, B, score)
persp(A, B, z, theta = -30, phi=30, col='green', 
      ticktype = 'detailed', zlim = c(0, 1), cex.axis=0.8)
# generate 2 random variables
set.seed(123)
A = round(runif(100, 0.1, 1), 1)
B = round(runif(100, 0.1, 1), 1)

y1 = 102 + 6*A + 4*B  + rnorm(100, 0, 5)
y0 = 100 + 3*A + 2*B  + rnorm(100, 0, 5)

ps = score(A, B)
d = rbinom(100, 1, ps)

y = y1*(d) + (1-d)*y0

dfData = data.frame(y, d=factor(d), A, B, y1, y0, p=ps)
mean(y1 - y0)

f1 = lm(y ~ d, data=dfData)
summary(f1)

f2 = lm(y ~ d + A + B, data=dfData)
summary(f2)

f3 = lm(y ~ d, weights = p, data=dfData)
summary(f3)

f4 = lm(y ~ d + A + B, weights=p, data=dfData)
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

f6 = lm(y ~ d + A + B, data=dfData.matched)
summary(f6)

f7 = lm(y ~ d, weights = p, data=dfData.matched)
summary(f7)

f8 = lm(y ~ d + A + B, weights=p, data=dfData.matched)
summary(f8)
