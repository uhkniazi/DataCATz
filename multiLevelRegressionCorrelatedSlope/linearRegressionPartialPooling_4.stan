data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ngroup1; // number of levels in group 1
  vector[Ntotal] X; // level predictor i.e. continuous variable
  int<lower=1, upper=Ngroup1> Ngroup1Map[Ntotal]; // mapping variable from group to unit levels
  real y[Ntotal]; // response variable normally distributed
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[2] coefGroup1[Ngroup1]; // correlated coefficients
  vector[2] MuPopGrp1; // population distribution
  real<lower=-1, upper=1> rho; // correlation 
  real<lower=0> sigmaPop; // data standard deviation
  real<lower=0> sigmaRan1; // group level error for group 1
  real<lower=0> sigmaRan2; // group level error for group 2
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  for (i in 1:Ntotal){
    mu[i] = coefGroup1[Ngroup1Map[i]][1] + X[i] * coefGroup1[Ngroup1Map[i]][2];
  }
}
model {
  matrix[2, 2] mCov; // covariance matrix
  // using diffuse prior
  sigmaRan1 ~ cauchy(0, 2);
  sigmaRan2 ~ cauchy(0, 2);
  sigmaPop ~ cauchy(0, 2);
  mCov[1, 2] = rho*sigmaRan1*sigmaRan2;
  mCov[2, 1] = mCov[1, 2];
  mCov[1, 1] = sigmaRan1^2;
  mCov[2, 2] = sigmaRan2^2;
  // sample from the multi normal distribution
  for (i in 1:Ngroup1){
    coefGroup1[i] ~ multi_normal(MuPopGrp1, mCov);
  }
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
