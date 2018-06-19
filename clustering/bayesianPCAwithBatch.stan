data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nvars; // number of variables/vectors in data matrix
  int<lower=1> Neigens; // number of eigen vectors to use
  vector[Nvars] y[Ntotal]; // data matrix where each vector is in column format
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each observation to group 1
  // additional parameters
  real gammaShape; // hyperparameters for the gamma distribution 
  real gammaRate;
}

transformed data {
  // create an identity covariance matrix for rotated components
  // fixed hyperparameter for the covariance matrix 
  matrix[Neigens, Neigens] mIdentity;
  // fixed hyperparameter for mean
  vector[Neigens] mu0;
  mIdentity = diag_matrix(rep_vector(1.0, Neigens));
  mu0 = rep_vector(0.0, Neigens);
}

parameters { // the parameters to track
  // eigen vectors of the covariance matrix of the data matrix
  matrix[Nvars, Neigens] mEigens; // note it is the transformation matrix 
  vector[Nvars] mu; // vector of means added to the centered rotated components 
  real<lower=0> sigma; // scale parameter for likelihood function
  vector<lower=0>[Neigens] sigma2; // scale parameter for the eigen vectors
  vector[Neigens] mComponents[Ntotal]; // vector with Ntotal components for the rotated data with dimensions NTotal, Neigens
  // variables for random effect section
  vector[1] beta; // regression intercept 
  real<lower=0> sigmaRan1; // random effect standard deviation for group 1
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
}
transformed parameters {
  row_vector[Ntotal] muDeflection; // fitted values from linear predictor
  muDeflection = (beta[1] + rGroupsJitter1[NgroupMap1])';
}

model {
  // prior for the sigma i.e. standard deviation of Eigens and likelihood
  sigma ~ cauchy(0, 2.5); // weak prior
  sigma2 ~ cauchy(0, 2.5); // weak prior
  // standard deviation of Deflections
  sigmaRan1 ~ gamma(gammaShape, gammaRate);
  beta ~ cauchy(0, 2.5);
  // random effects sample
  rGroupsJitter1 ~ normal(0, sigmaRan1);
  
  for (i in 1:Neigens){
    mEigens[,i] ~ normal(0.0, sigma2[i]);  
  }
  
  // go through each sample/observation and generate
  // a sample for each rotated components
  for (i in 1:Ntotal){
    mComponents[i] ~ multi_normal(mu0, mIdentity);
  }
  
  // now take a sample for the data as a function of
  // Eigens X Rotations + Centering
  // i.e. Operations X Inputs
  for (i in 1:Ntotal){
    y[i] ~ normal( mEigens * mComponents[i] + mu + muDeflection[i], sigma);
  }
}

