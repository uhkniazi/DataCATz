data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> NgroupMap[Ntotal]; // mapping variable to map each observation to a group 
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable normally distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
    real<lower=0> sigmaPop[2]; // 2 population standard deviations
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * betas; 
}
model {
  // using diffuse prior
  sigmaPop ~ cauchy(0, 2);
  betas ~ cauchy(0, 2); //prior for the betas
  // likelihood function
  y ~ normal(mu, sigmaPop[NgroupMap]);
}