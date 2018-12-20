data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ngroup1; // number of levels in group 1
  vector[Ntotal] X; // unit level predictor i.e. continuous variable
  vector[Ntotal] Xg; // group level predictor, continuous variable
  int<lower=1, upper=Ngroup1> Ngroup1Map[Ntotal]; // mapping variable from group to unit levels
  real y[Ntotal]; // response variable normally distributed
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  matrix[2, Ngroup1] coefGroup1; // correlated coefficients
  vector[4] populationCoeff;
  cholesky_factor_corr[2] cholGroup1;
  real<lower=0> sigmaPop; // data standard deviation
  vector<lower=0>[2] sigmaRan1; // group level error for group 1
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  //matrix[2, Ngroup1] MuGroup1Means; // distribution centers for the deflections
  matrix[2, Ngroup1] coefGroup1_adjusted;
  coefGroup1_adjusted = diag_pre_multiply(sigmaRan1, cholGroup1) * coefGroup1;
  // // group level model
  // for (i in 1:Ngroup1){
  //   MuGroup1Means[1, i] = populationCoeff[1] + populationCoeff[2] * Xg[i];
  //   MuGroup1Means[2, i] = populationCoeff[3] + populationCoeff[4] * Xg[i];
  // }
  for (i in 1:Ntotal){
    mu[i] = populationCoeff[1] + coefGroup1_adjusted[1][Ngroup1Map[i]] + // intercept + adjustment
    populationCoeff[2]*Xg[i]  +  // group level predictor slope
    X[i] * populationCoeff[3] +  // unit level slope
    X[i] * coefGroup1_adjusted[2][Ngroup1Map[i]] + // slope adjustment  
    populationCoeff[4] * X[i] * Xg[i];
  }
}
model {
  // using diffuse prior
  sigmaRan1 ~ cauchy(0, 2);
  sigmaPop ~ cauchy(0, 2);
  cholGroup1 ~ lkj_corr_cholesky(2.0);
  to_vector(populationCoeff) ~ cauchy(0, 2);
  to_vector(coefGroup1) ~ normal(0, 1);
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
