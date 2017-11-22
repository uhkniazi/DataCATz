data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> NgroupsLvl1; // number of levels for the level 1 of hierarchy - population level
  int<lower=1, upper=NgroupsLvl1> NgroupsLvl2Map[Ntotal]; // mapping variable to map each value at lvl 2 of hierarchy to lvl 1 variables
  int y[Ntotal]; // number of success
  int<lower=1> N[Ntotal]; // total trials, at least one trial
}

parameters { // the parameters to track
  real<lower=0, upper=1> theta[Ntotal];
  real<lower=0, upper=1> alpha1[NgroupsLvl1];
  real<lower=0> beta1[NgroupsLvl1];
}
// transformed parameters {
// 
// }
model {
  // distribution for level paremeter - used in the likelihood
  theta ~ beta(alpha1[NgroupsLvl2Map], beta1[NgroupsLvl2Map]);
  ////// likelihood
  y ~ binomial(N, theta);
}
