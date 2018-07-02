# Name: incumbent.R
# Auth: uhkniazi
# Date: 26/06/2018
# Desc: fitting various regression models, following the example from Gelman 2013

######################### fit a simple regression model first after loading the 1988 dataset
dfData = dfData = read.table('incumbentGelman/dataExternal/1988.asc', header=F)
colnames(dfData) = c('state', 'district', 'incumbency', 'democratVotes', 'republicanVotes')
str(dfData)

## previous election
dfData.prev = read.table('incumbentGelman/dataExternal/1986.asc', header=F)
colnames(dfData.prev) = c('state', 'district', 'incumbency', 'democratVotes', 'republicanVotes')


getWinner = function(d, r){
  w = -1
  if ((d == -9 | r == -9) & (d > 0 | r > 0)) {
    w = ifelse(d > r, 'd', 'r')
  } else if ((d > 0 & r > 0)) {
    w = ifelse(d > r, 'd', 'r')
  }
  return(w)
}

getTreatment = function(inc) {
  ifelse(inc == 0, 'open', 'incumbent')
}

# getResponse = function(inc, d, r, winner){
#   y = -1
#   if (inc == 0) {
#     # open seat
#     if (winner == 'd') y = d/(d+r) else y = r/(d+r)
#   } else if (inc == 1){
#     # not open seat
#     y = d / (d+r)
#   } else if (inc == -1) {
#     y = r / (d + r)
#   }
#   return(y)
# }

getResponse = function(inc, d, r){
  y = -1
  if (inc == 'r') {
    y = r / (d + r)
  } else if (inc == 'd'){
    y = d / (d+r)
  }
  return(y)
}

# check if the samples i.e. units of analysis are concordant in the 2 datasets
identical(dfData$state, dfData.prev$state)
identical(dfData$district, dfData.prev$district)
## get the winner for previous year
w = sapply(1:nrow(dfData.prev), function(x) getWinner(dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x]))
table(w)
## drop -1s as these are not concordant results
i = which(w == -1)
dfData = dfData[-i,]
dfData.prev = dfData.prev[-i,]
identical(dfData$state, dfData.prev$state)
identical(dfData$district, dfData.prev$district)

## perform same check for the current year
w = sapply(1:nrow(dfData), function(x) getWinner(dfData$democratVotes[x], dfData$republicanVotes[x]))
table(w)
## drop -1s as these are not concordant results
i = which(w == -1)
dfData = dfData[-i,]
dfData.prev = dfData.prev[-i,]
identical(dfData$state, dfData.prev$state)
identical(dfData$district, dfData.prev$district)

## get the winner of previous year as that is the current incumbent party
w = sapply(1:nrow(dfData.prev), function(x) getWinner(dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x]))
table(w)
dfData$incumbentParty = factor(w)

## get the treatment variable i.e. if open seat or incumbent candidate running
t = getTreatment(dfData$incumbency)
table(t)
dfData$treatment = factor(t)

## get previous election results for current incumbent party, i.e. winner of previous elections
w = sapply(1:nrow(dfData.prev), function(x) getWinner(dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x]))
table(w)
p = sapply(1:nrow(dfData.prev), function(x){
  getResponse(w[x], dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x])
}) 
summary(p)
dfData$previousProportion = p

## get the current year proportion for the incumbent party
p = sapply(1:nrow(dfData), function(x){
  getResponse(dfData$incumbentParty[x], dfData$democratVotes[x], dfData$republicanVotes[x])
}) 
summary(p)
dfData$response = p

# get the democratic vote for current and previous years
prev = dfData.prev$democratVotes / (dfData.prev$democratVotes + dfData.prev$republicanVotes)
cur = dfData$democratVotes / (dfData$democratVotes + dfData$republicanVotes)
summary(prev)
summary(cur)

## produce figure 14.1
plot(prev, cur, pch=c(20,2)[as.numeric(dfData$treatment)], xlab='1986', ylab='1988',
     xlim=c(0,1), ylim=c(0,1), main='Democratic Vote in 1988 vs 1986')


