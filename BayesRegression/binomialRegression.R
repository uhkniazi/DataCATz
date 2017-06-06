# Name: binomialRegression.R
# Auth: uhkniazi
# Date: 06/06/2017
# Desc: Binomial regression, using fixed and random effects (2 level model)


## using the contraception data set
## see page 120 Chapter 6 of Lme4 book.
library(lattice)
library(lme4)
library(mlmRev)

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

