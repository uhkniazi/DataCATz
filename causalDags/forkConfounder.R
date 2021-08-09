# Name: forkConfounder.R
# Auth: uhkniazi
# Date: 31/7/2021
# Desc: Showing a fork of mutual causation

library(dagitty)

dag.1 <- dagitty('dag {
  bb="0,0,1,1"
  C [pos="0.448,0.331"]
  D [exposure,pos="0.305,0.566"]
  Y [outcome,pos="0.628,0.567"]
  C -> D
  C -> Y
  D -> Y
}'
)

# simulate the data
set.seed(123)
C = rnorm(30, 5)
D = 5 * C + rnorm(30)
Y = 5 * C + 0 * D + rnorm(30)

pairs(dfData)

dfData = data.frame(C, D, Y)

fit.1 = lm(Y ~ D + C, data=dfData)
summary(fit.1)

fit.2 = lm(Y ~ D, data=dfData)
summary(fit.2)

fit.3 = lm(Y ~ 1, data=dfData)
summary(fit.3)

anova(fit.1, fit.2, fit.3)

