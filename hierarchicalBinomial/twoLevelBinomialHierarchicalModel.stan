data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nlevels2; // number of levels for the level 2 of hierarchy
  int<lower=1> Nlevels1; // number of levels for population i.e. level 1
  int<lower=1, upper=Nlevels2> NgroupMap2[Ntotal]; // mapping variable to map each observation to level 2
  int<lower=1, upper=Nlevels1> NgroupMap1[Ntotal]; // mapping variable to map each observation to level 1
  int y[Ntotal]; // number of success
  int N[Ntotal]; // total trials
}

parameters { // the parameters to track
  real theta[NgroupsLvl2];
  real omega1[NgroupsLvl1];
  real<lower=2> kappa1[NgroupsLvl1];
  real<lower=2> kappa0;
  real omega0;
}
// transformed parameters {
// 
// }
model {
  /////// hyper-hyperparameter distributions
  // using Jeffreys non informative prior
  omega0 ~ beta(0.5, 0.5);
  kappa0 ~ gamma(0.5, 1e-4);
  ////// hyperparameter distributions
  omega1 ~ beta(omega0*(kappa0-2)+1, (1-omega0)*(kappa0-2)+1);
  // using Jeffreys non informative prior for gamma distribution
  kappa1 ~ gamma(0.5, 1e-4);
  // 1st level prior for theta
  theta ~ beta(omega1*(kappa1-2)+1, (1-omega1)*(kappa1-2)+1);
  ////// likelihood
  y ~ binomial(Ntotal, theta);
}