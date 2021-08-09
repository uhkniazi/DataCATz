# Name: InvertedForkCollider.R
# Auth: uhkniazi
# Date: 31/7/2021
# Desc: Showing a inverted fork with collider

library(dagitty)

dag.1 <- dagitty('dag {
bb="0,0,1,1"
C [pos="0.448,0.331"]
D [exposure,pos="0.305,0.566"]
Y [outcome,pos="0.628,0.567"]
D -> C
D -> Y
Y -> C
}
'
)
plot(dag.1)
# simulate the data
set.seed(123)
D = rnorm(30, 10, 2)
Y = rnorm(30, 6)
C = 5 * D + 3 * Y + rnorm(30)

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
