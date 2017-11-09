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
  // omega1 ~ beta(omega0, kappa0)
  real<lower=2> kappa0; // concentration parameter for top level beta distribution
  real omega0; // mean parameter for top level beta distribution
  // theta ~ beta(omega1, kappa1)
  real<lower=2> kappa1[Nlevels1]; // concentration parameters for level 1 beta distributions
  real omega1[Nlevels1]; // concentration parameters for level 1 beta distributions
  real<lower=0, upper=1> theta[Nlevels2]; // theta parameters 
    vector[Ncol] betas; // regression parameters
    real<lower=-1> sigma; // scale parameter for normal distribution  
  }
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X*betas; // fitted value using identity link
}
model {
  y ~ normal( mu , exp(sigma));
}