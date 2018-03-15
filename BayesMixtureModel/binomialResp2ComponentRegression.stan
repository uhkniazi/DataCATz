data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol1; // total number of columns in model matrix 1
    int<lower=1> Ncol2; // total number of columns in model matrix 2
    matrix[Ntotal, Ncol1] X1; // model matrix for 1st component
    matrix[Ntotal, Ncol2] X2; // model matrix for 2nd component
    int<lower=2> iMixtures; // number of mixture distributions
    int y[Ntotal]; // response variable binomial distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
  
  vector[Ncol1] betasMix1; // regression parameters for each mixture component
  vector[Ncol2] betasMix2; // regression parameters for each mixture component
  vector[iMixtures] mu; // ordered intercept
}
transformed parameters {
  vector[Ntotal] muMix1; // number of fitted values
  vector[Ntotal] muMix2; // number of fitted values
  // calculate fitted values without intercept
  muMix1 = X1 * betasMix1;
  muMix2 = X2 * betasMix2;
  //mu = inv_logit(mu);
}
model {
  real ps[iMixtures]; // temporary variable for log components
  vector[iMixtures] iMixWeights; // weights for the number of mixtures (should sum to one)
  iMixWeights = rep_vector(0.5, 2);
  // any priors go here, get some information from ones we fit earlier
  mu[1] ~ cauchy(0, 2.5);
  mu[2] ~ cauchy(0, 2.5);
  betasMix1 ~ cauchy(0, 2.5);
  betasMix2 ~ cauchy(0, 2.5);
  iMixWeights ~ dirichlet(rep_vector(2.0, iMixtures));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // second loop for number of mixture components
    ps[1] = log(iMixWeights[1]) + bernoulli_lpmf(y[n] | inv_logit(mu[1]+muMix1[n]));
    ps[2] = log(iMixWeights[2]) + bernoulli_lpmf(y[n] | inv_logit(mu[2]+muMix2[n]));
    target += log_sum_exp(ps);
  }
}
