data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ngroup1; // number of levels in group 1
 // int<lower=1> Ngroup2; // number of levels in group 2
  vector[Ntotal] X; // level predictor i.e. continuous variable
  int<lower=1, upper=Ngroup1> Ngroup1Map[Ntotal]; // mapping variable from group to unit levels
  //int<lower=1, upper=Ngroup2> Ngroup2Map[Ntotal]; // mapping variable from group to unit levels
  real y[Ntotal]; // response variable normally distributed
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  // vector[Ngroup1] coefGroup1;
  // vector[Ngroup2] coefGroup2;
  vector[2] coefGroup1[Ngroup1]; // correlated coefficients
  vector[2] MuPopGrp1; // population distribution
  // real MuPopGrp1;
  // real MuPopGrp2;
  real<lower=-1, upper=1> rho; // correlation 
  real<lower=0> sigmaPop; // data standard deviation
  real<lower=0> sigmaRan1; // group level error for group 1
  real<lower=0> sigmaRan2; // group level error for group 2
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  //vector[Ntotal] coefGroup1_slope;
  //coefGroup1_slope = coefGroup1[2][Ngroup1Map];
  for (i in 1:Ntotal){
    mu[i] = coefGroup1[Ngroup1Map[i]][1] + X[i] * coefGroup1[Ngroup1Map[i]][2];
  }
  //mu = mu + coefGroup1[Ngroup1Map][1];
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
  //coefGroup1 ~ normal(MuPopGrp1, sigmaRan1); //distribution for the intercepts 
  //coefGroup2 ~ normal(MuPopGrp2, sigmaRan2); //distribution for the slopes
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
