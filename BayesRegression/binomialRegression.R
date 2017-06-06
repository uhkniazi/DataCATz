# Name: binomialRegression.R
# Auth: uhkniazi
# Date: 06/06/2017
# Desc: Binomial regression, using fixed and random effects (2 level model)


## using the contraception data set
## see page 120 Chapter 6 of Lme4 book.
library(lattice)
library(lme4)
library(mlmRev)
library(car)
library(LearnBayes)

data("Contraception")

str(Contraception)

## we plot the data as the author does in this book
xyplot(ifelse(use == 'Y', 1, 0) ~ age|urban, data = Contraception, groups = livch, lty=1:4, col=1,
       type=c('g', 'smooth'), key = list(text=list(levels(Contraception$livch)), lines=list(lty=1:4),
                                         columns=4, cex=0.8),
       xlab='Age', ylab='Contraception Use', main=list(label='Contraception use Given Urban', cex=0.8))

## Figure suggests following things:
## 1- Contraception use and age have a quadradic trend.
## 2- middle age range women are more likely to use contraception
## 3- Urban women are more likely to use contraception
## 4- women with 0 childeren are less likely to use contraception.
## 5- women who have children are less different from each other

## using point 5, merge the groups 1, 2, and 3 in the variable livch
Contraception$children = factor(Contraception$livch != 0, labels=c('N', 'Y'))

str(Contraception)

xyplot(ifelse(use == 'Y', 1, 0) ~ age|urban, data = Contraception, groups = children, lty=1:2, col=1,
       type=c('g', 'smooth'), key = list(text=list(levels(Contraception$children)), lines=list(lty=1:2),
                                         columns=2, cex=0.8),
       xlab='Age', ylab='Contraception Use', main=list(label='Contraception use Given Urban', cex=0.8))

## we should perhaps also use an interaction term between age and number of children, as the slopes and interceptrs
## in cases with children are different to those without children.

## lets fit a model using the binomial glm without random effects
fit.1 = glm(use ~ age + I(age^2) + urban + children, data=Contraception, family = binomial(link='logit'))
summary(fit.1)


logit.inv = function(p) {exp(p)/(exp(p)+1) }

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
  lp = dcauchy(betas[1], 0, 10, log=T) + sum(dcauchy(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

lData = list(resp=ifelse(Contraception$use == 'N', 0, 1), mModMatrix=model.matrix(use ~ age + I(age^2) + urban + children, data=Contraception))
start = c(rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.2 = laplace(mylogpost, start, lData)
se = sqrt(diag(fit.2$var))

