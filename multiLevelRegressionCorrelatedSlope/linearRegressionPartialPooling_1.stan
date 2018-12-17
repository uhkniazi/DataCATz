data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  matrix[Ntotal, Ncol] X; // model matrix
  real y[Ntotal]; // response variable normally distributed
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real populationMean;
  real<lower=0> sigmaPop; // data standard deviation
  real<lower=0> sigmaRan; // group level error
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * betas; 
}
model {
  // using diffuse prior
  sigmaRan ~ cauchy(0, 2);
  sigmaPop ~ cauchy(0, 2);
  betas ~ normal(populationMean, sigmaRan); //prior for the betas
  // likelihood function
  y ~ normal(mu, sigmaPop);
}
