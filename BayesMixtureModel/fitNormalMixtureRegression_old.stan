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
    matrix[Ncol, iMixtures] betas; // regression parameters for each mixture component
  }
transformed parameters { // calculated parameters 
    matrix[Ntotal, Ncol] mu; // number of fitted values
    mu = X * betas;
}
model {
  // see stan manual page 187 for an example
  real ps[iMixtures]; // temporary variable for log components
  // any priors go here 
  //mu ~ normal(100, 12);
  sigma ~ cauchy(0, 2.5); // weak prior
  iMixWeights ~ dirichlet(rep_vector(2.0, iMixtures));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // second loop for number of mixture components
    for (k in 1:iMixtures){
      ps[k] = log(iMixWeights[k]) + normal_lpdf(y[n] | mu[,k], sigma[k]);
    }
    target += log_sum_exp(ps);
  }
}
