data {
    int<lower=1> Ntotal; // number of observations
    real y[Ntotal]; // observed data
}
parameters {
    real mu; // posterior mean
    real<lower=0> sig; // posterior sd
}
model {
  y ~ normal(mu, sig);
}