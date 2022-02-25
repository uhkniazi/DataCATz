# Name: matchingAndWeightedRegressionStanExample.R
# Auth: uhkniazi
# Date: 18/2/2022
# Desc: An example of matching and weighted regression including Stan example

# simulate the data for matching demonstration 3 page 153
# and page 229 - see morgan and winship books

# propensity score function
score = function(A, B) {
  r = -2 +3*(A) -3*(A-0.1) +2*(A-.3) -2*(A-.5) +4*(A-.7)
  -4*(A-.9) +1*(B) -1*(B-.1) +2*(B-.7) -2*(B -.9) +3*(A-.5)*(B -.5)
  -3*(A-.7)*(B -.7)
  return(plogis(r))
}

# score = function(A, B){
#   return(rbeta(length(A), 10, 10))
# }

# generate 2 random variables
set.seed(123)
A = round(runif(1000, 0.1, 1), 1)
B = round(runif(1000, 0.1, 1), 1)

y1 = 102 + 6*A + 4*B  + rnorm(1000, 0, 5)
y0 = 100 + 3*A + 2*B  + rnorm(1000, 0, 5)

## treatment assignment propensity
ps = score(A, B)
d = rbinom(1000, 1, ps)

## treatment assignment
y = y1*(d) + (1-d)*y0

dfData = data.frame(y, d=factor(d), A, B, y1, y0, p=ps,
                    p.ate=ifelse(d==1, 1/ps, 1/(1-ps)),
                    p.att=ifelse(d==1, 1, ps/(1-ps)),
                    p.atc=ifelse(d==1, (1-ps)/ps, 1))
mean(y1 - y0)
mean(y1) - mean(y0)

# try different models
f1 = lm(y ~ d, data=dfData)
summary(f1)

f2 = lm(y ~ d + A + B, data=dfData)
summary(f2)

f3 = lm(y ~ d, weights = p.ate, data=dfData)
summary(f3)

f4 = lm(y ~ d + A + B, weights=p.ate, data=dfData)
summary(f4)

library(arm)
matches = matching(z = d, score = dfData$p) 
matches$match.ind
matches$cnts
matches$pairs
matched = as.logical(matches$cnts)
dfData.matched = dfData[matched,]
dim(dfData.matched)

f5 = lm(y ~ d, data=dfData.matched)
summary(f5)

f6 = lm(y ~ d + A + B, data=dfData.matched)
summary(f6)

f7 = lm(y ~ d, weights = p.ate, data=dfData.matched)
summary(f7)

f8 = lm(y ~ d + A + B, weights=p.ate, data=dfData.matched)
summary(f8)

### try weighted regression with stan
## fit model with stan with various model sizes
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rethinking)

stanDso = rstan::stan_model(file='normResponseWeighted.stan')

str(dfData)
m1 = model.matrix(y ~ d + A + B, data=dfData)
head(m1)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m1), X=m1,
                 y=dfData$y, wts=dfData$p.ate)

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('betas', 'sigmaPop', 'mu', 'log_lik'),
                      cores=2)
print(fit.stan, c('sigmaPop', 'betas'), digits=3)

traceplot(fit.stan.4, 'populationMean')
traceplot(fit.stan.2, 'sigmaPop')
traceplot(fit.stan.2, 'betas')
