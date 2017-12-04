data {
  int<lower=1> Ntotal; // number of observations
  int<lower=Ntotal> Ngroups1; // number of levels for the level 1 of hierarchy - group level
  int<lower=1, upper=Ngroups1> NgroupsMap[Ntotal]; // mapping variable to map data points to estimated parameters
  int y[Ntotal]; // number of success
  int<lower=1> N[Ntotal]; // total trials, at least one trial required
}

parameters { // the parameters to track
  real<lower=0, upper=1> theta[Ntotal];   // parameter of interest for level 1 of hierarachy
  real<lower=0, upper=1> AlDivBe; // transformed hyper-parameter or population level parameter  
  real<lower=1> AlPlusBe; // transformed hyper-parameter 
}
transformed parameters {
  real alpha;
  real beta;
  alpha = AlPlusBe / (1 + (1/AlDivBe));
  beta = AlPlusBe - alpha;
}
model {
  AlDivBe ~ beta(0.5, 0.5);
  AlPlusBe ~ gamma(0.5, 1e-4);
  // distribution for population paremeter with hyperparameters alpha and beta
  theta ~ beta(alpha, beta);
  ////// likelihood
  y ~ binomial(N, theta[NgroupsMap]);
}

