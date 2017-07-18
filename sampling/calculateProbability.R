# Name: calculateProbability.R
# Auth: uhkniazi
# Date: 18/7/2017
# Desc: calculate probabilities of events by simulation


######################################
sim.game = function() {
  r = sample(1:6, 2, replace=T)
  if (sum(r) == 7 || sum(r) == 11) return(1) else return(0)
}
win = replicate(1000, sim.game())
sum(win)
sum(win)/1000












## urn with N=4, M=2
urn = c('R', 'R', 'W', 'W')

mSims = sapply(1:1000, function(x) sample(urn, size = 3, replace = F))
mSims = t(mSims)

## how many times R occurs on the 2nd or 3rd draws but not on both
R.2 = mSims[,2] == 'R'
sum(R.2)/1000

R.3 = mSims[,3] == 'R'
sum(R.3)/1000

R.1 = mSims[,1] == 'R'
sum(R.1)/1000

sum(R.1 & R.2)/1000

R.2.3 = xor(R.2, R.3)
sum(R.2.3)/1000
sum(R.2 & R.3)/1000
sum(R.2 | R.3)/1000

## how many times R occurs on 1st draw while on also 2nd or 3rd draw
R.1 = mSims[,1] == 'R'
R.1And2Or3 = (R.1 & R.2.3)  
sum(R.1And2Or3)/10000

## how many times R occurs on the 1st and second draws
sum(R.1 & R.2)/10000

sum(mSims[,2] == 'R' & mSims[,3]=='W')/10000
sum(mSims[,2] == 'W' & mSims[,3]=='R')/10000
sum(mSims[,2] == 'R' & mSims[,3]=='R')/10000
sum(mSims[,2] == 'W' & mSims[,3]=='W')/10000

(1 * sum(mSims[,2] == 'R' & mSims[,3]=='W')/10000) + (1 * sum(mSims[,2] == 'W' & mSims[,3]=='R')/10000) +
  (0 * sum(mSims[,2] == 'R' & mSims[,3]=='R')/10000) + (2 * sum(mSims[,2] == 'W' & mSims[,3]=='W')/10000)
mean(R.1And2Or3)
mean(R.1)