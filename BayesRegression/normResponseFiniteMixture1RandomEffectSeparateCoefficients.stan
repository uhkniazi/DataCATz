data {
  int<lower=1> Ntotal; // number of observations
  real y[Ntotal]; // response variable - normally distributed
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=2> iMixtures; // number of mixture distributions
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1 
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
  real iIntercepts[iMixtures];
}

parameters { // the parameters to track
  ordered[iMixtures] mu; // number of means to track Breaking the Labeling Degeneracy by Enforcing an Ordering (population Intercepts)
  real<lower=0> sigma[iMixtures]; // scale parameters for normal distribution (population sigmas) 
  simplex[iMixtures] iMixWeights; // weights for the number of mixtures (should sum to one)
  // regression coefficients and other related parameters
  real<lower=0> sigmaRan1_1; // random effect standard deviation for group 1
  real<lower=0> sigmaRan1_2; // random effect standard deviation for group 1 component 2
  vector[Nclusters1] rGroupsJitter1_1; // number of random jitters for each level of cluster/group 1 - first set of regression coefficients
  vector[Nclusters1] rGroupsJitter1_2; // number of random jitters for each level of cluster/group 1 - second set 
}
transformed parameters {
  vector[Ntotal] muFitted_1; // fitted value from linear predictor for first component of model
  vector[Ntotal] muFitted_2; // fitted value from linear predictor for second component of model
  muFitted_1 =  rGroupsJitter1_1[NgroupMap1]; // we use separate regression coefficients for each component
  muFitted_2 =  rGroupsJitter1_2[NgroupMap1];
}
model {
  // see stan manual page 187 for an example
  real ps[iMixtures]; // temporary variable for log components
  sigmaRan1_1 ~ gamma(gammaShape, gammaRate); // scale parameter for the random jitter
  sigmaRan1_2 ~ gamma(gammaShape, gammaRate); // scale parameter for the random jitter
  rGroupsJitter1_1 ~ normal(0, sigmaRan1_1);
  rGroupsJitter1_2 ~ normal(0, sigmaRan1_2);
  // any priors for mixture components go here 
  mu[1] ~ normal(iIntercepts[1], 1);
  mu[2] ~ normal(iIntercepts[2], 1);
  sigma ~ gamma(gammaShape, gammaRate);
  iMixWeights ~ dirichlet(rep_vector(2.0, iMixtures));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // second loop for number of mixture components
    ps[1] = log(iMixWeights[1]) + normal_lpdf(y[n] | mu[1] + muFitted_1[n], sigma[1]);
    ps[2] = log(iMixWeights[2]) + normal_lpdf(y[n] | mu[2] + muFitted_2[n], sigma[2]);
    target += log_sum_exp(ps);
  }
}
