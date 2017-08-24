data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nvars; // number of variables/vectors in data matrix
  int<lower=1> Neigens; // number of eigen vectors to use
  matrix[Nvars, Ntotal] y; // data matrix where each vector is in row format
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
  vector[Neigens] mComponents[Ntotal]; // matrix for the rotated data in row vector format
}

model {
  // prior for the sigma
  sigma ~ cauchy(0, 2.5); // weak prior
  // go through each sample/observation and generate
  // a sample for each rotated components
  for (i in 1:Ntotal){
    mComponents[i] ~ multi_normal(mu0, mIdentity);
  }
  
  // now take a sample for the data as a function of
  // Eigens X Rotations + Centering
  // i.e. Operations X Inputs
  for (i in 1:Ntotal){
    y[,i] ~ normal( mEigens * mComponents[i] + mu , sigma);
  }
}