data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable - normally distributed
}

parameters { // the parameters to track
    vector[Ncol] betas; // regression parameters
    real<lower=-1> sigma; // scale parameter for normal distribution  
  }
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X*betas; // fitted value using identity link
}
model {
  y ~ normal( mu , exp(sigma));
}