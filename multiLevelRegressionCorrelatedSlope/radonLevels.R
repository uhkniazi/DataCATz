# Name: radonLevels.R
# Auth: umar niazi 
# Date: 13/12/2018
# Desc: multilevel regression example with correlated intercept and slope from Gelman Multilevel modelling book 2006



# load and format the data
# see instructions for data download here
# http://wiki.math.yorku.ca/index.php/Book:_Gelman_%26_Hill_%282007%29#Data_sets
# data formatting
# http://www.stat.columbia.edu/~gelman/arm/examples/radon/radon_setup.R
# Set up the radon data
setwd('multiLevelRegressionCorrelatedSlope/')
## this section is from http://www.stat.columbia.edu/~gelman/arm/examples/radon/radon_setup.R
# read in and clean the data

srrs2 <- read.table ("srrs2.dat", header=T, sep=",")
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
x <- floor

# get county index variable

county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

# get the county-level predictor

srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)
## end section for data formatting


# Complete-pooling and no-pooling estimates of county radon levels
# Consider the goal of estimating the distribution of
# radon levels of the houses within each of the 85 counties in Minnesota.
# page 252 onwards Gelman 2006 book