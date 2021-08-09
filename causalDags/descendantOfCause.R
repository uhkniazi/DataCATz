# Name: descendantOfCause.R
# Auth: uhkniazi
# Date: 9/8/2021
# Desc: adjusting for descendant of a cause

library(dagitty)

dag.1 <- dagitty('dag {
bb="0,0,1,1"
C [pos="0.620,0.401"]
D [exposure,pos="0.178,0.483"]
Y [outcome,pos="0.411,0.446"]
D -> Y
Y -> C
}
'
)
plot(dag.1)
# simulate the data
set.seed(123)
D = rnorm(30)
Y = D * 6 + rnorm(30)
C = 3 * Y + rnorm(30)

dfData = data.frame(C, D, Y)
pairs(dfData)
cor(Y, D)
cor(dfData)

fit.1 = lm(Y ~ D + C, data=dfData)
summary(fit.1)

fit.2 = lm(Y ~ D, data=dfData)
summary(fit.2)

fit.3 = lm(Y ~ 1, data=dfData)
summary(fit.3)

anova(fit.1, fit.2, fit.3)

fit.4 = lm(C ~ D + Y, data=dfData)
summary(fit.4)
