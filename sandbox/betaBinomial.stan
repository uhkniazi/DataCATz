data {
  int<lower=1> iTotal; // number of observations
  int iResponse[iTotal]; // response variable - success per observation
  int iTrials[iTotal]; // trials per observation
}
parameters {
  // see here for choice of prior
  // https://discourse.mc-stan.org/t/prior-choice-for-beta-binomial-dispersion/24562/3
  // https://groups.google.com/g/stan-users/c/4iUO6TGopbs
  // https://github.com/paul-buerkner/brms/issues/377
  real<lower=1e-6, upper=0.99> rRho; // dispersion parameter = (1-rho)/rho
  real<lower=-15, upper=0> rFitted; // transform to proportion later
}
transformed parameters {
  real alpha;
  real beta;
  real rProportion; // probability of success  
  real rPhi; // dispersion parameter
  rPhi = (1-rRho)/rRho;
  rProportion = inv_logit(rFitted); 
  alpha = rProportion * rPhi;
  beta = (1-rProportion) * rPhi;
}
model {
  //rFitted ~ norm(, 1);
  rRho ~ beta(0.5, 0.5);
  // likelihood
  iResponse ~ beta_binomial(iTrials, alpha, beta);
}
// generated quantities{
//   vector[iTotal] log_lik;
//   for (i in 1:iTotal) {
//     log_lik[i] = beta_binomial_lpmf(iResponse[i] | iTrials[i], alpha[i], beta[i]);
//   }
// }

