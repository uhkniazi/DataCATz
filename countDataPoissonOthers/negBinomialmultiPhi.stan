data {
  int<lower=1> Ntotal; // number of observations
  //vector<lower=0>[Ntotal] mu; // fitted value
  int<lower=0> y[Ntotal]; // response count variable 
}
parameters {
  // parameters to estimate in the model
  //real<lower=0, upper=1> phi_inv[Ntotal]; // over dispersion term
  real<lower=0, upper=1> phi_inv; // over dispersion term
  real mu;
}
// transformed parameters {
//   real phi[Ntotal];
//   phi = 1/phi_inv;
// }
model {
  phi_inv ~ uniform(0, 1);
  // likelihood function
  // for (i in 1:Ntotal){
  //   y[i] ~ neg_binomial_2(mu[i], 1/phi_inv[i]);  
  // }
  y ~ neg_binomial_2(mu, 1/phi_inv);
}
