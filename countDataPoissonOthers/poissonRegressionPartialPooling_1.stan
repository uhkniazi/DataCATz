data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  matrix[Ntotal, Ncol] X; // model matrix
  vector[Ntotal] offset; // log of exposure
  int<lower=0> y[Ntotal]; // response count variable 
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real populationMean;
  real<lower=0.01> sigmaRan; // group level error
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * ((betas * sigmaRan)+populationMean); 
  mu = mu + offset;
}
model {
  sigmaRan ~ cauchy(0, 2);
  betas ~ normal(0, 1); //prior for the betas
  populationMean ~ cauchy(0, 2);
  // likelihood function
  y ~ poisson(exp(mu));
}
