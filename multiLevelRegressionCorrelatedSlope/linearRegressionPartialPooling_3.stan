data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ngroup1; // number of levels in group 1
  int<lower=1> Ngroup2; // number of levels in group 2
  vector[Ntotal] X; // level predictor i.e. continuous variable
  int<lower=1, upper=Ngroup1> Ngroup1Map[Ntotal]; // mapping variable from group to unit levels
  int<lower=1, upper=Ngroup2> Ngroup2Map[Ntotal]; // mapping variable from group to unit levels
  real y[Ntotal]; // response variable normally distributed
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ngroup1] coefGroup1;
  vector[Ngroup2] coefGroup2;
  real MuPopGrp1;
  real MuPopGrp2;
  real<lower=0> sigmaPop; // data standard deviation
  real<lower=0> sigmaRan1; // group level error for group 1
  real<lower=0> sigmaRan2; // group level error for group 2
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ntotal] coefGroup2_expanded;
  coefGroup2_expanded = coefGroup2[Ngroup2Map];
  for (i in 1:Ntotal){
    mu[i] = X[i] * coefGroup2_expanded[i];
  }
  mu = mu + coefGroup1[Ngroup1Map];
}
model {
  // using diffuse prior
  sigmaRan1 ~ cauchy(0, 2);
  sigmaRan2 ~ cauchy(0, 2);
  sigmaPop ~ cauchy(0, 2);
  coefGroup1 ~ normal(MuPopGrp1, sigmaRan1); //distribution for the intercepts 
  coefGroup2 ~ normal(MuPopGrp2, sigmaRan2); //distribution for the slopes
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
