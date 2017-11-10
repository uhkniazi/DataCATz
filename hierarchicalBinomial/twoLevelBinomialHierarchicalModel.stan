data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> NgroupsLvl1; // number of levels for the level 1 of hierarchy - population level
  int<lower=1> NgroupsLvl2; // number of levels for level 2
  int<lower=1, upper=NgroupsLvl1> NgroupsLvl2Map[NgroupsLvl2]; // mapping variable to map each value at lvl 1 of hierarchy to lvl 2 variables
  int<lower=1, upper=NgroupsLvl2> NdataMap[Ntotal]; // mapping variable to map each level 2 prior to the data
  int y[Ntotal]; // number of success
  int N[Ntotal]; // total trials
}

parameters { // the parameters to track
  real<lower=0, upper=1> theta[NgroupsLvl2];
  real<lower=0, upper=1> omega1[NgroupsLvl1];
  real<lower=2> kappa1[NgroupsLvl1];
  real<lower=2> kappa0;
  real<lower=0, upper=1> omega0;
}
// transformed parameters {
// 
// }
model {
  /////// level 0 - hyper-hyperparameter distributions
  // using Jeffreys non informative prior
  omega0 ~ beta(0.5, 0.5);
  kappa0 ~ gamma(0.5, 1e-4);
  ////// level 1 - hyperparameter distributions
  omega1 ~ beta(omega0*(kappa0-2)+1, (1-omega0)*(kappa0-2)+1);
  // using Jeffreys non informative prior for gamma distribution
  kappa1 ~ gamma(0.5, 1e-4);
  ///// level 2
  // distribution for level 2 paremeter - used in the likelihood
  for (i in 1:NgroupsLvl2){
    theta[i] ~ beta(omega1[NgroupsLvl2Map[i]]*(kappa1[NgroupsLvl2Map[i]]-2)+1, (1-omega1[NgroupsLvl2Map[i]])*(kappa1[NgroupsLvl2Map[i]]-2)+1);
  }
  
  ////// likelihood
  y ~ binomial(N, theta[NdataMap]);
}