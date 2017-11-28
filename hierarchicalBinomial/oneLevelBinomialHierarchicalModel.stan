data {
  int<lower=1> Ntotal; // number of observations
  int<lower=Ntotal> Ngroups1; // number of levels for the level 1 of hierarchy - group level
  int<lower=1, upper=Ngroups1> NgroupsMap[Ntotal]; // mapping variable to map data points to estimated parameters
  int y[Ntotal]; // number of success
  int<lower=1> N[Ntotal]; // total trials, at least one trial required
}

parameters { // the parameters to track
  real<lower=0, upper=1> theta[Ntotal];   // parameter of interest for level 1 of hierarachy
  real<lower=0, upper=1> alpha; // hyper-parameter or population level parameter  
  real<lower=0> beta; // hyper-parameter 
}
// transformed parameters {
// 
// }
model {
  // distribution for level paremeter - used in the likelihood
  theta ~ beta(alpha, beta);
  ////// likelihood
  y ~ binomial(N, theta[NgroupsMap]);
}

