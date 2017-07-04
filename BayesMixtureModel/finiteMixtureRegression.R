# File: finiteMixtureRegression.R
# Auth: Umar Niazi,
# Date: 4/7/2017
# Desc: Fitting a finite mixture regression model

library(flexmix)

data("NPreg")
str(NPreg)

fit.flex = flexmix(yn ~ x, data=NPreg, k=2)
## fitted coefficients
parameters(fit.flex)
## predicted values, both return the same values
p = predict(fit.flex, newdata=NPreg)
p = do.call(cbind, p)
f = fitted(fit.flex)

## they will however return 2 components, 
head(f)

## if we want to aggregate these values, we need a weighted average, based on
## the weights assigned to the mixture components
## following the source code from the flexmix package
## https://rdrr.io/cran/flexmix/src/R/flexmix.R
# if (aggregate) {
#   prior_weights <- prior(object, newdata)
#   z <- lapply(x, function(z) matrix(rowSums(do.call("cbind", z) * prior_weights), nrow = nrow(z[[1]])))
# }
## lets check
p2 = predict(fit.flex, newdata=NPreg, aggregate=T)
p2 = do.call(cbind, p2)
head(p2)

## take a simple averga
head(rowMeans(p)) ## not the same
pr = prior(fit.flex)
p.agg = sweep(p, 2, pr, '*')
p.agg = rowSums(p.agg)
head(p.agg)
identical(as.numeric(p.agg), as.numeric(p2))



