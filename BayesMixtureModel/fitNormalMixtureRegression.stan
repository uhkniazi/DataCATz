data {
    int<lower=1> Ntotal; // number of observations
    real y[Ntotal]; // response variable - normally distributed
    int<lower=2> iMixtures; // number of mixture distributions
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
}

parameters { // the parameters to track
    real<lower=0> sigma[iMixtures]; // scale parameters for normal distribution  
    simplex[iMixtures] iMixWeights; // weights for the number of mixtures (should sum to one)
    vector[(Ncol-1)] betasMix1; // regression parameters for each mixture component
    vector[(Ncol-1)] betasMix2; // regression parameters for each mixture component
    ordered[iMixtures] mu; // ordered intercept
  }
transformed parameters { // calculated parameters
    vector[Ntotal] muMix1; // number of fitted values
    vector[Ntotal] muMix2; // number of fitted values
    matrix[Ntotal, (Ncol-1)] mX2; // new model matrix without intercept
    mX2 = X[,2:Ncol]; 
    // calculate fitted values without intercept
    muMix1 = mX2 * betasMix1;
    muMix2 = mX2 * betasMix2;
}
model {
  // see stan manual page 187 for an example
  real ps[iMixtures]; // temporary variable for log components
  // any priors go here, get some information from ones we fit earlier
  mu[1] ~ normal(10, 10);
  mu[2] ~ normal(40, 10);
  betasMix1 ~ cauchy(0, 10);
  betasMix2 ~ cauchy(0, 10);
  sigma ~ cauchy(0, 2.5); // weak prior
  iMixWeights ~ dirichlet(rep_vector(2.0, iMixtures));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // second loop for number of mixture components
    ps[1] = log(iMixWeights[1]) + normal_lpdf(y[n] | mu[1]+muMix1[n], sigma[1]);
    ps[2] = log(iMixWeights[2]) + normal_lpdf(y[n] | mu[2]+muMix2[n], sigma[2]);
    target += log_sum_exp(ps);
  }
}
