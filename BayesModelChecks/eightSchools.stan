data {
    int<lower=1> Ntotal; // number of observations
    real y[Ntotal]; // mean for each school
    real<lower=0> sig[Ntotal]; // the standard error for each school
  }
parameters {
    real mu; // the hypermarameter mean
    real<lower=0> tau; // hyperparameter standard deviation
    vector[Ntotal] eta; // deflections
  }
transformed parameters {
  vector[Ntotal] theta; // school level effects or deflections
  theta = mu + tau*eta;
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sig);
}