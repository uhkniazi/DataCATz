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

getResponse = function(inc, d, r, winner){
  y = -1
  if (inc == 0) {
    # open seat
    if (winner == 'd') y = d/(d+r) else y = r/(d+r)
  } else if (inc == 1){
    # not open seat
    y = d / (d+r)
  } else if (inc == -1) {
    y = r / (d + r)
  }
  return(y)
}


## get the winner for current year
w = sapply(1:nrow(dfData), function(x) getWinner(dfData$democratVotes[x], dfData$republicanVotes[x]))
table(w)
i = which(w == -1)
dfData = dfData[-i,]
dfData.prev = dfData.prev[-i,]

dfData$winner = factor(w[-i])
dfData$treatment = factor(getTreatment(dfData$incumbency))
dfData$response = sapply(1:nrow(dfData), function(x){
  getResponse(dfData$incumbency[x], dfData$democratVotes[x], dfData$republicanVotes[x], dfData$winner[x])
}) 

dfData$oldResult = sapply(1:nrow(dfData), function(x){
  getResponse(dfData$incumbency[x], dfData.prev$democratVotes[x], dfData.prev$republicanVotes[x], dfData$winner[x])
}) 

