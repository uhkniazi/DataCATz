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
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * betas; 
  mu = mu + offset;
}
model {
  betas ~ cauchy(0, 2); //prior for the betas
  // likelihood function
  y ~ poisson(exp(mu));
}
