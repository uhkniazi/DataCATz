data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    int y[Ntotal]; // response variable binomial distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X * betas; 
  mu = inv_logit(mu);
}
model {
  betas ~ cauchy(0, 10); //prior for the betas
  
  // likelihood function
  y ~ bernoulli(mu);
}