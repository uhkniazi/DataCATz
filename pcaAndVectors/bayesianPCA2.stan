data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nvars; // number of variables/vectors in data matrix
  int<lower=1> Neigens; // number of eigen vectors to use
  vector[Nvars] y[Ntotal]; // data matrix where each vector is in column format
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
}

model {
  // prior for the sigma
  sigma ~ cauchy(0, 2.5); // weak prior
  sigma2 ~ cauchy(0, 2.5); // weak prior
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
    y[i] ~ normal( mEigens * mComponents[i] + mu , sigma);
  }
}

